#!/usr/bin/env python3
"""
Group near-duplicate food names and aggregate their scientific names.

Steps:
- Clean names (drop "raw", normalise separators, remove '?', strip accents).
- Exact match by cleaned common/scientific name.
- Fuzzy match with character-bigram Dice similarity (length-aware thresholds).
- Optionally ask ChatGPT about borderline pairs.
- Aggregate all common + scientific names per group and write a 4-column CSV.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import unicodedata
from collections import defaultdict, OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import requests


RAW_WORD_RE = re.compile(r"\braw\b", flags=re.IGNORECASE)
NON_ALNUM_SPACE = re.compile(r"[^a-z0-9 ]+")
MULTISPACE = re.compile(r"\s+")
LOG_PREFIX = "[match]"


def log(message: str) -> None:
    print(f"{LOG_PREFIX} {message}", flush=True)


def load_env_file(path: Path) -> None:
    """Lightweight .env loader (KEY=VALUE). Existing env vars win."""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        key = key.strip()
        if key and key not in os.environ:
            os.environ[key] = value.strip()


def strip_accents(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    return "".join(ch for ch in normalized if not unicodedata.combining(ch))


def clean_display_name(value: str) -> str:
    """
    Lightly clean for output: drop 'raw', normalise separators, trim spacing, remove '?'.
    """
    value = (value or "").replace("_", " ").replace("-", " ")
    value = RAW_WORD_RE.sub("", value)
    value = value.replace("?", " ")
    value = MULTISPACE.sub(" ", value).strip()
    return value


def normalise_for_similarity(value: str) -> str:
    value = clean_display_name(value)
    value = strip_accents(value).lower()
    value = NON_ALNUM_SPACE.sub(" ", value)
    value = MULTISPACE.sub(" ", value).strip()
    return value


def collapse_no_space(value: str) -> str:
    return value.replace(" ", "")


def char_bigrams(value: str) -> Set[str]:
    if len(value) < 2:
        return {value} if value else set()
    return {value[i : i + 2] for i in range(len(value) - 1)}


def dice_similarity(a: Set[str], b: Set[str]) -> float:
    if not a or not b:
        return 0.0
    overlap = len(a & b)
    return (2.0 * overlap) / (len(a) + len(b))


def tokens(value: str) -> Tuple[str, ...]:
    if not value:
        return tuple()
    return tuple(value.split())


def pick_threshold(max_len: int) -> Tuple[float, float]:
    """
    Decide similarity thresholds based on length.
    Returns (merge_threshold, ai_threshold_lower_bound).
    """
    if max_len >= 8:
        return 0.9, 0.8
    return 0.6, 0.5


def pepper_conflict(tokens_a: Tuple[str, ...], tokens_b: Tuple[str, ...], norm_a: str, norm_b: str) -> bool:
    """
    Prevent merging different pepper types unless the non-pepper portion matches.
    """
    has_pepper = "pepper" in tokens_a or "pepper" in tokens_b
    if not has_pepper:
        return False
    if norm_a == norm_b:
        return False
    base_a = tuple(tok for tok in tokens_a if tok != "pepper")
    base_b = tuple(tok for tok in tokens_b if tok != "pepper")
    return base_a != base_b


@dataclass
class Record:
    row_index: int
    source_index: str
    raw_common: str
    raw_scientific: str
    row_data: Dict[str, str]
    clean_common: str
    clean_scientific: str
    norm_common: str
    norm_chars: str
    tokens: Tuple[str, ...]
    bigrams: Set[str]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Match similar food names and aggregate scientific names.")
    parser.add_argument(
        "-i",
        "--input",
        dest="input_csv",
        default="ndm_foods.csv",
        help="Input CSV to process (default: ndm_foods.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_csv",
        help="Output CSV path (default: <input>_matched.csv)",
    )
    parser.add_argument(
        "--use-openai",
        action="store_true",
        help="Ask the ChatGPT API about borderline pairs (requires OPENAI_API_KEY).",
    )
    parser.add_argument(
        "--openai-model",
        default="gpt-4o-mini",
        help="ChatGPT model to use when --use-openai is set (default: gpt-4o-mini).",
    )
    parser.add_argument(
        "--openai-batch-size",
        type=int,
        default=8,
        help="Number of borderline pairs to send per ChatGPT request (default: 8).",
    )
    parser.add_argument(
        "--max-openai-pairs",
        type=int,
        default=250,
        help="Safety cap on how many borderline pairs to send to ChatGPT (default: 250).",
    )
    parser.add_argument(
        "--collapse-groups",
        action="store_true",
        help="Write one row per merged group (row count shrinks) instead of one row per original index.",
    )
    return parser.parse_args(argv)


def load_records(path: Path) -> Tuple[str, List[str], List[Record]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV missing headers")
        index_field = reader.fieldnames[0]
        data_fields = [field for field in reader.fieldnames if field != index_field]

        records: List[Record] = []
        seen_rows: set[Tuple[str, ...]] = set()
        total_rows = 0
        duplicates_skipped = 0
        syn_all_field = "synonyms_all_sources"

        for row in reader:
            total_rows += 1
            dedupe_key = tuple((row.get(field, "") or "").strip() for field in data_fields)
            if dedupe_key in seen_rows:
                duplicates_skipped += 1
                continue
            seen_rows.add(dedupe_key)

            source_index = row.get(index_field, "")
            raw_common = row.get("food_com", "") or ""
            raw_scientific = row.get("food_sci", "") or ""

            clean_common = clean_display_name(raw_common or raw_scientific)
            clean_scientific = clean_display_name(raw_scientific)
            norm_common = normalise_for_similarity(raw_common or raw_scientific)
            norm_chars = collapse_no_space(norm_common)

            row_data = {field: row.get(field, "") or "" for field in data_fields}
            # Normalize synonyms_all_sources spacing if present
            if syn_all_field in row_data:
                row_data[syn_all_field] = "; ".join(
                    [part.strip() for part in row_data[syn_all_field].split(";") if part.strip()]
                )

            rec = Record(
                row_index=len(records),
                source_index=source_index,
                raw_common=raw_common,
                raw_scientific=raw_scientific,
                row_data=row_data,
                clean_common=clean_common,
                clean_scientific=clean_scientific,
                norm_common=norm_common,
                norm_chars=norm_chars,
                tokens=tokens(norm_common),
                bigrams=char_bigrams(norm_chars),
            )
            records.append(rec)

    # For exact common-name duplicates, keep only rows with the richest synonyms.
    # Priority: highest synonyms_all_sources count; tie-break by total synonyms across all sources; then first seen.
    grouped: Dict[str, List[Record]] = defaultdict(list)
    for rec in records:
        grouped[rec.clean_common].append(rec)

    filtered_records: List[Record] = []
    synonym_pref_dropped = 0
    for recs in grouped.values():
        if len(recs) == 1:
            filtered_records.append(recs[0])
            continue
        scores = []
        for r in recs:
            syn_all_count = count_synonyms(r.row_data.get(syn_all_field, ""))
            total_syn = (
                syn_all_count
                + count_synonyms(r.row_data.get("synonyms_open_tree_of_life", ""))
                + count_synonyms(r.row_data.get("synonyms_wiki_search", ""))
                + count_synonyms(r.row_data.get("synonyms_ncbi", ""))
            )
            has_sci = 1 if r.raw_scientific else 0
            scores.append((syn_all_count, total_syn, has_sci))

        max_score = max(scores)
        keep_idx = next(i for i, s in enumerate(scores) if s == max_score)
        keep_rec = recs[keep_idx]
        synonym_pref_dropped += len(recs) - 1
        filtered_records.append(keep_rec)

    # Reindex row_index after filtering
    for new_idx, rec in enumerate(filtered_records):
        rec.row_index = new_idx

    records = filtered_records

    log(
        f"Loaded {len(records)} unique rows from {path} "
        f"(input rows: {total_rows}, skipped {duplicates_skipped} exact duplicates, "
        f"dropped {synonym_pref_dropped} by synonym richness)"
    )
    return index_field, data_fields, records


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, idx: int) -> int:
        if self.parent[idx] != idx:
            self.parent[idx] = self.find(self.parent[idx])
        return self.parent[idx]

    def union(self, a: int, b: int) -> None:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a == root_b:
            return
        if root_a < root_b:
            self.parent[root_b] = root_a
        else:
            self.parent[root_a] = root_b


def exact_pass(records: List[Record], uf: UnionFind) -> int:
    buckets: Dict[str, List[int]] = defaultdict(list)
    for idx, rec in enumerate(records):
        if rec.norm_common:
            buckets[rec.norm_common].append(idx)
    merges = 0
    for indices in buckets.values():
        if len(indices) < 2:
            continue
        first = indices[0]
        for other in indices[1:]:
            uf.union(first, other)
            merges += 1
    log(f"Exact pass buckets: {len(buckets)}, unions applied: {merges}")
    return merges


def fuzzy_pass(records: List[Record], uf: UnionFind) -> Tuple[List[Tuple[int, int, float, float]], int]:
    """
    Compare within prefix buckets to avoid O(n^2) matching.
    Returns borderline pairs that may be sent to the model.
    """
    borderline: List[Tuple[int, int, float, float]] = []
    buckets: Dict[str, List[int]] = defaultdict(list)
    merges = 0

    for idx, rec in enumerate(records):
        key = rec.norm_chars[:2] if len(rec.norm_chars) >= 2 else rec.norm_chars
        buckets[key].append(idx)

    for candidates in buckets.values():
        size = len(candidates)
        if size < 2:
            continue
        for i in range(size):
            a_idx = candidates[i]
            rec_a = records[a_idx]
            for j in range(i + 1, size):
                b_idx = candidates[j]
                if uf.find(a_idx) == uf.find(b_idx):
                    continue
                rec_b = records[b_idx]
                if not rec_a.norm_chars or not rec_b.norm_chars:
                    continue
                if abs(len(rec_a.norm_chars) - len(rec_b.norm_chars)) > 3:
                    continue

                if pepper_conflict(rec_a.tokens, rec_b.tokens, rec_a.norm_common, rec_b.norm_common):
                    continue

                sim = dice_similarity(rec_a.bigrams, rec_b.bigrams)
                max_len = max(len(rec_a.norm_chars), len(rec_b.norm_chars))
                threshold, ai_lower = pick_threshold(max_len)

                if sim >= threshold:
                    uf.union(a_idx, b_idx)
                    merges += 1
                elif ai_lower <= sim < threshold:
                    borderline.append((a_idx, b_idx, sim, threshold))
    log(f"Fuzzy pass buckets: {len(buckets)}, unions applied: {merges}, borderline queued: {len(borderline)}")
    return borderline, merges


def chunked(sequence: Sequence[Tuple[int, int, float, float]], size: int) -> Iterable[List[Tuple[int, int, float, float]]]:
    for start in range(0, len(sequence), size):
        yield list(sequence[start : start + size])


def ask_chatgpt_about_pairs(
    pairs: List[Tuple[int, int, float, float]],
    records: List[Record],
    model: str,
    api_key: str,
    batch_size: int,
    max_pairs: int,
) -> Set[Tuple[int, int]]:
    """Send borderline pairs to ChatGPT and return those it said are the same."""
    approved: Set[Tuple[int, int]] = set()
    endpoint = "https://api.openai.com/v1/chat/completions"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}

    limited_pairs = pairs[:max_pairs]

    for batch in chunked(limited_pairs, batch_size):
        lines = []
        for idx, (a_idx, b_idx, sim, threshold) in enumerate(batch, start=1):
            a_name = records[a_idx].clean_common or records[a_idx].clean_scientific
            b_name = records[b_idx].clean_common or records[b_idx].clean_scientific
            lines.append(f"{idx}. {a_name!r} vs {b_name!r}")

        prompt = (
            "Decide if each pair of food names refers to the same food item. "
            'Return ONLY JSON in the shape {"answers": [true|false,...]} matching the order below '
            "(true = same food, false = different). Be conservative about merging peppers unless the type matches.\n"
            + "\n".join(lines)
        )
        payload = {
            "model": model,
            "messages": [
                {
                    "role": "system",
                    "content": (
                        "You resolve whether two food names refer to the same food item. "
                        "Be conservative about merging peppers unless they are clearly the exact type."
                    ),
                },
                {"role": "user", "content": prompt},
            ],
            "response_format": {"type": "json_object"},
        }

        try:
            response = requests.post(endpoint, headers=headers, json=payload, timeout=40)
            response.raise_for_status()
            data = response.json()
            content = data["choices"][0]["message"]["content"]
            parsed = json.loads(content)
            if isinstance(parsed, dict):
                answers = parsed.get("answers", [])
            elif isinstance(parsed, list):
                answers = parsed
            else:
                answers = []
        except Exception as exc:  # noqa: BLE001
            print(f"WARNING: ChatGPT call failed, skipping batch: {exc}")
            continue

        if not isinstance(answers, list):
            print("WARNING: Unexpected ChatGPT response shape; skipping batch.")
            continue

        for decision, (a_idx, b_idx, _, _) in zip(answers, batch):
            if isinstance(decision, bool) and decision:
                approved.add((a_idx, b_idx))

    return approved


def build_groups(records: List[Record], uf: UnionFind) -> Tuple[Dict[int, List[int]], Dict[int, int]]:
    groups: Dict[int, List[int]] = defaultdict(list)
    index_to_root: Dict[int, int] = {}
    for idx, _ in enumerate(records):
        root = uf.find(idx)
        groups[root].append(idx)
        index_to_root[idx] = root
    return groups, index_to_root


def dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    ordered = OrderedDict()
    for val in values:
        if not val:
            continue
        if val not in ordered:
            ordered[val] = None
    return list(ordered.keys())


def count_synonyms(value: str) -> int:
    return sum(1 for part in value.split(";") if part.strip())


def row_score(row: Dict[str, str]) -> Tuple[int, int]:
    syn_all = count_synonyms(row.get("synonyms_all_sources", ""))
    total_syn = (
        syn_all
        + count_synonyms(row.get("synonyms_open_tree_of_life", ""))
        + count_synonyms(row.get("synonyms_wiki_search", ""))
        + count_synonyms(row.get("synonyms_ncbi", ""))
    )
    has_sci = 1 if (row.get("food_sci") or "").strip() else 0
    return syn_all, total_syn + has_sci


def filter_by_common_best(rows: List[Dict[str, str]]) -> Tuple[List[Dict[str, str]], int]:
    buckets: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("food_com") or "").strip()
        buckets[key].append(row)

    filtered: List[Dict[str, str]] = []
    dropped = 0
    for items in buckets.values():
        if len(items) == 1:
            filtered.append(items[0])
            continue
        scores = [row_score(r) for r in items]
        max_score = max(scores)
        keep_idx = next(i for i, s in enumerate(scores) if s == max_score)
        filtered.append(items[keep_idx])
        dropped += len(items) - 1
    return filtered, dropped


def write_output(
    output_path: Path,
    index_field: str,
    data_fields: Sequence[str],
    records: List[Record],
    groups: Dict[int, List[int]],
    index_to_root: Dict[int, int],
) -> None:
    fieldnames = [index_field] + list(data_fields)

    # For each cleaned common name, keep the richest synonyms_all_sources; clear others.
    max_synonyms_by_common: Dict[str, int] = defaultdict(int)
    syn_field = "synonyms_all_sources"
    for rec in records:
        if not rec.clean_common:
            continue
        count = count_synonyms(rec.row_data.get(syn_field, ""))
        if count > max_synonyms_by_common[rec.clean_common]:
            max_synonyms_by_common[rec.clean_common] = count

    rows: List[Dict[str, str]] = []
    for rec in records:
        rec_idx = rec.row_index
        root = index_to_root.get(rec_idx, rec_idx)
        members = groups.get(root, [rec_idx])

        # Aggregate names for this group
        common_names = dedupe_preserve_order(
            [records[m].clean_common for m in members if records[m].clean_common]
        )
        resolved = rec.clean_common or rec.clean_scientific

        row_out = {index_field: rec.source_index, **rec.row_data}
        if "food_com" in row_out:
            row_out["food_com"] = "; ".join(common_names) if common_names else resolved
        if "food_sci" in row_out:
            row_out["food_sci"] = rec.clean_scientific
        if syn_field in row_out:
            max_syns = max_synonyms_by_common.get(rec.clean_common, 0)
            if count_synonyms(row_out[syn_field]) < max_syns:
                row_out[syn_field] = ""

        rows.append(row_out)

    rows, dropped = filter_by_common_best(rows)
    if dropped:
        log(f"Dropped {dropped} rows after output common-name filtering")

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)


def write_collapsed_output(
    output_path: Path,
    index_field: str,
    data_fields: Sequence[str],
    records: List[Record],
    groups: Dict[int, List[int]],
) -> None:
    """
    Emit one row per merged group.
    """
    fieldnames = [index_field] + list(data_fields)
    rows: List[Dict[str, str]] = []

    syn_field = "synonyms_all_sources"

    for root, members in groups.items():
        # Use the smallest source index as the representative
        representative = min((records[m].source_index for m in members if records[m].source_index), default=str(root))
        common_names = dedupe_preserve_order(
            [records[m].clean_common for m in members if records[m].clean_common]
        )
        resolved = common_names[0] if common_names else records[members[0]].clean_common

        # Pick member with the richest synonyms_all_sources (fallback to first)
        richest_member = members[0]
        best_count = -1
        for m in members:
            count = count_synonyms(records[m].row_data.get(syn_field, ""))
            if count > best_count:
                best_count = count
                richest_member = m

        representative_row = {index_field: representative, **records[richest_member].row_data}
        if "food_com" in representative_row:
            representative_row["food_com"] = "; ".join(common_names) if common_names else resolved
        if "food_sci" in representative_row:
            representative_row["food_sci"] = records[richest_member].clean_scientific
        rows.append(representative_row)

    rows, dropped = filter_by_common_best(rows)
    if dropped:
        log(f"Dropped {dropped} rows after collapsed common-name filtering")

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> None:
    load_env_file(Path(".env"))

    args = parse_args(argv)
    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise SystemExit(f"Input CSV not found: {input_path}")

    default_output_name = f"{input_path.stem}_matched_collapsed.csv" if args.collapse_groups else f"{input_path.stem}_matched.csv"
    output_path = Path(args.output_csv) if args.output_csv else input_path.with_name(default_output_name)

    index_field, data_fields, records = load_records(input_path)

    uf = UnionFind(len(records))
    exact_merges = exact_pass(records, uf)
    borderline_pairs, fuzzy_merges = fuzzy_pass(records, uf)
    log(
        f"After matching: exact merges={exact_merges}, fuzzy merges={fuzzy_merges}, "
        f"borderline queued={len(borderline_pairs)}"
    )

    if args.use_openai:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise SystemExit("OPENAI_API_KEY not set; required for --use-openai.")
        approved = ask_chatgpt_about_pairs(
            borderline_pairs,
            records,
            model=args.openai_model,
            api_key=api_key,
            batch_size=max(1, args.openai_batch_size),
            max_pairs=max(1, args.max_openai_pairs),
        )
        log(f"ChatGPT approved {len(approved)} borderline merges (of {len(borderline_pairs)})")
        for a_idx, b_idx in approved:
            uf.union(a_idx, b_idx)
    else:
        log("Skipping ChatGPT step (no --use-openai)")

    groups, index_to_root = build_groups(records, uf)
    log(f"Built {len(groups)} groups out of {len(records)} rows")
    if args.collapse_groups:
        write_collapsed_output(output_path, index_field, data_fields, records, groups)
        log(f"Wrote collapsed (deduped) CSV to {output_path}")
    else:
        write_output(output_path, index_field, data_fields, records, groups, index_to_root)
        log(f"Wrote resolved CSV to {output_path}")
    log(f"Columns: {[index_field] + data_fields}")


if __name__ == "__main__":
    main()
