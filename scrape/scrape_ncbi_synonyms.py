from __future__ import annotations

import csv
import time
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
from xml.etree import ElementTree

import requests

DEFAULT_INPUT_PATH = Path("ndm_foods.csv")

USER_AGENT = "gFoodsScraper/0.1 (+https://example.com/contact)"
NCBI_TOOL = "gFoodsScraper"
NCBI_EMAIL = "example@example.com"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

REQUEST_HEADERS = {"User-Agent": USER_AGENT}
REQUEST_TIMEOUT = 25
THROTTLE_SECONDS = 0.34  # stay comfortably within NCBI rate limits

OPEN_TREE_FIELD = "synonyms_open_tree_of_life"
WIKI_FIELD = "synonyms_wiki_search"
NCBI_FIELD = "synonyms_ncbi"

PROGRESS_PREFIX = "[ncbi]"

SESSION = requests.Session()


def format_scientific_name(name: str) -> str:
    return name.replace("_", " ").strip()


def format_common_name(name: str) -> str:
    return name.replace("_", " ").strip()


def build_key(scientific: str, common: str, fallback_index: int) -> Tuple[str, str]:
    sci = format_scientific_name(scientific)
    if sci:
        return ("sci", sci)

    com = format_common_name(common)
    if com:
        return ("com", com)

    return ("row", str(fallback_index))


def candidate_queries(scientific: str, common: str) -> List[str]:
    candidates: List[str] = []
    sci = format_scientific_name(scientific)
    com = format_common_name(common)

    if sci:
        candidates.append(sci)
        if sci.endswith("."):
            candidates.append(sci[:-1])
        lower = sci.lower()
        if lower.endswith(" sp.") or lower.endswith(" sp") or lower.endswith(" spp.") or lower.endswith(" spp"):
            candidates.append(sci.rsplit(" ", 1)[0])

    if com:
        candidates.append(com)

    return list(OrderedDict.fromkeys([c for c in candidates if c]))


def search_taxid(query: str) -> Optional[str]:
    params = {
        "db": "taxonomy",
        "term": f"{query}[All Names]",
        "retmode": "xml",
        "tool": NCBI_TOOL,
        "email": NCBI_EMAIL,
    }
    try:
        response = SESSION.get(
            ESEARCH_URL, params=params, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"{PROGRESS_PREFIX} WARNING: esearch failed for '{query}': {exc}")
        return None

    root = ElementTree.fromstring(response.text)
    id_list = root.find("IdList")
    if id_list is None:
        return None

    for id_elem in id_list.findall("Id"):
        if id_elem.text:
            return id_elem.text.strip()
    return None


def extract_synonyms(xml_text: str) -> List[str]:
    root = ElementTree.fromstring(xml_text)
    taxon = root.find("Taxon")
    if taxon is None:
        return []

    other_names = taxon.find("OtherNames")
    if other_names is None:
        return []

    values: List[str] = []
    synonym_tags = [
        "Synonym",
        "EquivalentName",
        "GenbankSynonym",
        "GenbankCommonName",
        "CommonName",
    ]

    for tag in synonym_tags:
        for elem in other_names.findall(tag):
            text = elem.text.strip() if elem.text else ""
            if text:
                values.append(text)

    for name_elem in other_names.findall("Name"):
        disp = name_elem.findtext("DispName")
        if disp:
            values.append(disp.strip())

    # remove duplicates while preserving order
    ordered = OrderedDict()
    for value in values:
        cleaned = value.strip()
        if cleaned and cleaned not in ordered:
            ordered[cleaned] = None

    return list(ordered.keys())


def fetch_synonyms_for_taxid(taxid: str) -> List[str]:
    params = {
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml",
        "tool": NCBI_TOOL,
        "email": NCBI_EMAIL,
    }
    try:
        response = SESSION.get(
            EFETCH_URL, params=params, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"{PROGRESS_PREFIX} WARNING: efetch failed for taxid {taxid}: {exc}")
        return []

    return extract_synonyms(response.text)


def main(input_path: Path, output_path: Path, limit: Optional[int] = None) -> None:
    if not input_path.exists():
        raise SystemExit(f"Input CSV not found: {input_path}")

    with input_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV missing headers")

        index_field = reader.fieldnames[0]
        rows = list(reader)

    if limit is not None:
        rows = rows[:limit]

    print(
        f"{PROGRESS_PREFIX} loaded {len(rows)} rows from {input_path}",
        flush=True,
    )

    name_candidates: Dict[Tuple[str, str], List[str]] = {}
    for idx, row in enumerate(rows):
        scientific_raw = row.get("food_sci", "")
        common_raw = row.get("food_com", "")
        key = build_key(scientific_raw, common_raw, idx)
        queries = candidate_queries(scientific_raw, common_raw)
        if queries:
            name_candidates.setdefault(key, queries)

    print(
        f"{PROGRESS_PREFIX} prepared {len(name_candidates)} unique name keys",
        flush=True,
    )

    query_cache: Dict[str, Optional[str]] = {}
    taxid_cache: Dict[Tuple[str, str], Optional[str]] = {}
    synonym_cache: Dict[str, str] = {}

    updated_rows: List[dict] = []

    total_keys = len(rows)
    for idx, row in enumerate(rows, start=1):
        if idx % 500 == 0 or idx == total_keys:
            print(
                f"{PROGRESS_PREFIX} processing row {idx}/{total_keys}",
                flush=True,
            )

        scientific_raw = row.get("food_sci", "")
        common_raw = row.get("food_com", "")
        key = build_key(scientific_raw, common_raw, idx - 1)
        queries = name_candidates.get(key, [])

        synonyms_ncbi = ""
        taxid: Optional[str] = None

        if key in taxid_cache:
            taxid = taxid_cache[key]
        else:
            for query in queries:
                if query not in query_cache:
                    query_cache[query] = search_taxid(query)
                    time.sleep(THROTTLE_SECONDS)
                taxid = query_cache[query]
                if taxid:
                    break
            taxid_cache[key] = taxid

        if taxid:
            if taxid not in synonym_cache:
                synonyms = fetch_synonyms_for_taxid(taxid)
                synonym_cache[taxid] = "; ".join(synonyms)
                time.sleep(THROTTLE_SECONDS)
            synonyms_ncbi = synonym_cache.get(taxid, "")

        updated_rows.append(
            {
                index_field: row.get(index_field, ""),
                "food_com": row.get("food_com", ""),
                "food_sci": row.get("food_sci", ""),
                OPEN_TREE_FIELD: row.get(OPEN_TREE_FIELD, ""),
                WIKI_FIELD: row.get(WIKI_FIELD, ""),
                NCBI_FIELD: synonyms_ncbi,
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

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(updated_rows)

    print(f"{PROGRESS_PREFIX} wrote results to {output_path}", flush=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Populate NCBI synonyms for foods CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help="Path to input CSV (default: ndm_foods.csv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path to output CSV (default: overwrite input)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of rows to process (for testing).",
    )

    args = parser.parse_args()

    output_path = args.output or args.input

    main(args.input, output_path, args.limit)
