from __future__ import annotations

import csv
import sys
import time
from collections import OrderedDict
from pathlib import Path

import requests

INPUT_PATH = Path("ndm_foods.csv")
OUTPUT_PATH = INPUT_PATH  # overwrite in place once synonyms are populated

USER_AGENT = "gFoodsScraper/0.1 (+https://example.com/contact)"
MATCH_URL = "https://api.opentreeoflife.org/v3/tnrs/match_names"

REQUEST_HEADERS = {"User-Agent": USER_AGENT}
REQUEST_TIMEOUT = 30
THROTTLE_SECONDS = 0.15
CHUNK_SIZE = 40

OPEN_TREE_FIELD = "synonyms_open_tree_of_life"
WIKI_FIELD = "synonyms_wiki_search"
NCBI_FIELD = "synonyms_ncbi"


def format_scientific_name(name: str) -> str:
    return name.replace("_", " ").strip()


def select_best_match(matches: list[dict]) -> dict | None:
    if not matches:
        return None

    sorted_matches = sorted(matches, key=lambda m: m.get("score", 0), reverse=True)

    for match in sorted_matches:
        if not match.get("is_synonym") and not match.get("is_approximate_match"):
            return match

    for match in sorted_matches:
        if not match.get("is_synonym"):
            return match

    return sorted_matches[0]


def fetch_synonyms_batch(names: list[str]) -> dict[str, list[str]]:
    payload = {"names": names, "include_synonyms": True}
    response = requests.post(
        MATCH_URL, json=payload, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
    )
    response.raise_for_status()

    data = response.json()
    results = data.get("results", [])
    synonyms_map: dict[str, list[str]] = {}

    for idx, name in enumerate(names):
        result = results[idx] if idx < len(results) else None
        if not result:
            synonyms_map[name] = []
            continue

        best_match = select_best_match(result.get("matches", []))
        if best_match:
            taxon = best_match.get("taxon", {})
            synonyms = taxon.get("synonyms") or []
            synonyms_map[name] = list(OrderedDict.fromkeys(synonyms))
        else:
            synonyms_map[name] = []

    return synonyms_map


def main() -> None:
    if not INPUT_PATH.exists():
        raise SystemExit(f"Input CSV not found: {INPUT_PATH}")

    with INPUT_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV missing headers")

        index_field = reader.fieldnames[0]
        rows = list(reader)

    formatted_names = [
        format_scientific_name(row.get("food_sci", "")) for row in rows
    ]
    unique_names = sorted({name for name in formatted_names if name})

    synonyms_cache: dict[str, str] = {}

    for start in range(0, len(unique_names), CHUNK_SIZE):
        batch = unique_names[start : start + CHUNK_SIZE]
        try:
            batch_synonyms = fetch_synonyms_batch(batch)
        except requests.RequestException as exc:
            print(
                f"WARNING: Failed to fetch synonyms for batch starting with {batch[0]}: {exc}",
                file=sys.stderr,
            )
            batch_synonyms = {name: [] for name in batch}

        for name, synonyms in batch_synonyms.items():
            synonyms_cache[name] = "; ".join(synonyms)

        time.sleep(THROTTLE_SECONDS)

    updated_rows: list[dict[str, str]] = []

    for row, formatted_name in zip(rows, formatted_names):
        scientific_raw = row.get("food_sci", "")
        synonyms_otol = synonyms_cache.get(formatted_name, "") if formatted_name else ""

        updated_rows.append(
            {
                index_field: row.get(index_field, ""),
                "food_com": row.get("food_com", ""),
                "food_sci": scientific_raw,
                OPEN_TREE_FIELD: synonyms_otol,
                WIKI_FIELD: row.get(WIKI_FIELD, ""),
                NCBI_FIELD: row.get(NCBI_FIELD, ""),
            }
        )

    fieldnames = [
        index_field,
        "food_com",
        "food_sci",
        OPEN_TREE_FIELD,
        WIKI_FIELD,
        NCBI_FIELD,
    ]

    with OUTPUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(updated_rows)


if __name__ == "__main__":
    main()
