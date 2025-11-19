from __future__ import annotations

import csv
import time
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import requests

INPUT_PATH = Path("ndm_foods.csv")
OUTPUT_PATH = INPUT_PATH

USER_AGENT = "gFoodsScraper/0.1 (+https://example.com/contact)"
SEARCH_URL = "https://www.wikidata.org/w/api.php"
ENTITY_URL = "https://www.wikidata.org/w/api.php"

REQUEST_HEADERS = {"User-Agent": USER_AGENT}
REQUEST_TIMEOUT = 25
SEARCH_LIMIT = 5
ENTITY_CHUNK_SIZE = 40
TITLE_CHUNK_SIZE = 40
THROTTLE_SECONDS = 0.02
PROGRESS_PREFIX = "[wiki]"

OPEN_TREE_FIELD = "synonyms_open_tree_of_life"
WIKI_FIELD = "synonyms_wiki_search"
NCBI_FIELD = "synonyms_ncbi"


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

    # ensure uniqueness for records lacking both fields
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

    # Deduplicate while preserving order
    return list(OrderedDict.fromkeys([c for c in candidates if c]))


SESSION = requests.Session()


def search_entity_id(query: str) -> Optional[str]:
    params = {
        "action": "wbsearchentities",
        "format": "json",
        "language": "en",
        "search": query,
        "limit": SEARCH_LIMIT,
    }

    try:
        response = SESSION.get(
            SEARCH_URL, params=params, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"WARNING: Wikidata search failed for '{query}': {exc}")
        return None

    data = response.json()
    results = data.get("search", [])
    if not results:
        return None

    query_lower = query.lower()

    # Prefer exact label matches
    for result in results:
        label = (result.get("label") or "").lower()
        if label == query_lower:
            return result["id"]

    # Prefer entries where the match text equals the query
    for result in results:
        match = result.get("match", {})
        match_text = (match.get("text") or "").lower()
        if match_text == query_lower:
            return result["id"]

    # Fallback to the highest-scoring result
    return results[0]["id"]


def fetch_entities(qids: Iterable[str]) -> Dict[str, dict]:
    qids_list = [qid for qid in qids if qid]
    entities: Dict[str, dict] = {}

    for start in range(0, len(qids_list), ENTITY_CHUNK_SIZE):
        batch = qids_list[start : start + ENTITY_CHUNK_SIZE]
        params = {
            "action": "wbgetentities",
            "format": "json",
            "ids": "|".join(batch),
            "props": "labels|aliases|claims",
            "languages": "en",
        }

        try:
            response = SESSION.get(
                ENTITY_URL, params=params, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"WARNING: Failed fetching entities {batch[0]}...: {exc}")
            continue

        data = response.json()
        entities.update(data.get("entities", {}))
        time.sleep(THROTTLE_SECONDS)

    return entities


def parse_normalized_map(normalized: Optional[dict]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if normalized is None:
        return mapping
    if isinstance(normalized, list):
        for item in normalized:
            if isinstance(item, dict) and "from" in item and "to" in item:
                mapping[item["from"]] = item["to"]
    elif isinstance(normalized, dict):
        for item in normalized.values():
            if isinstance(item, dict) and "from" in item and "to" in item:
                mapping[item["from"]] = item["to"]
    return mapping


def fetch_entities_for_titles(
    title_to_keys: Dict[str, List[Tuple[str, str]]]
) -> Tuple[
    Dict[Tuple[str, str], dict],
    Dict[Tuple[str, str], str],
    Dict[str, dict],
    Set[Tuple[str, str]],
]:
    titles = list(title_to_keys.keys())
    key_to_entity: Dict[Tuple[str, str], dict] = {}
    key_to_qid: Dict[Tuple[str, str], str] = {}
    qid_to_entity: Dict[str, dict] = {}
    unmatched: Set[Tuple[str, str]] = set()
    unresolved_titles: Set[str] = set()

    if not titles:
        return key_to_entity, key_to_qid, qid_to_entity, unmatched

    total_batches = (len(titles) + TITLE_CHUNK_SIZE - 1) // TITLE_CHUNK_SIZE

    for start in range(0, len(titles), TITLE_CHUNK_SIZE):
        batch = titles[start : start + TITLE_CHUNK_SIZE]
        batch_index = start // TITLE_CHUNK_SIZE + 1
        print(
            f"{PROGRESS_PREFIX} fetching title batch {batch_index}/{total_batches} "
            f"({len(batch)} titles)",
            flush=True,
        )
        params = {
            "action": "wbgetentities",
            "format": "json",
            "sites": "enwiki",
            "titles": "|".join(batch),
            "props": "aliases|claims|sitelinks",
            "languages": "en",
        }

        try:
            response = SESSION.get(
                ENTITY_URL, params=params, headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"WARNING: Failed fetching titles batch starting {batch[0]}: {exc}")
            unmatched.update([key for title in batch for key in title_to_keys[title]])
            continue

        data = response.json()
        title_entity_map: Dict[str, Tuple[str, dict]] = {}
        for entity_id, entity in data.get("entities", {}).items():
            if entity.get("missing") == "":
                title = entity.get("title")
                if title:
                    unresolved_titles.add(title)
                continue
            enwiki = entity.get("sitelinks", {}).get("enwiki", {})
            title = enwiki.get("title")
            if not title:
                continue
            underscore_title = title.replace(" ", "_")
            title_entity_map[underscore_title] = (entity_id, entity)
            qid_to_entity[entity_id] = entity

        for requested_title in batch:
            keys = title_to_keys.get(requested_title, [])
            if requested_title in title_entity_map:
                qid, entity = title_entity_map[requested_title]
                for key in keys:
                    key_to_entity[key] = entity
                    key_to_qid[key] = qid
            else:
                unresolved_titles.add(requested_title)

        time.sleep(THROTTLE_SECONDS)

    single_title_cache: Dict[str, Tuple[str, dict]] = {}

    unresolved_titles = sorted(unresolved_titles)
    for idx, title in enumerate(unresolved_titles, start=1):
        print(
            f"{PROGRESS_PREFIX} resolving leftover title {idx}/{len(unresolved_titles)}: {title}",
            flush=True,
        )
        if title in single_title_cache:
            result = single_title_cache[title]
        else:
            params = {
                "action": "wbgetentities",
                "format": "json",
                "sites": "enwiki",
                "titles": title,
                "props": "aliases|claims|sitelinks",
                "languages": "en",
                "normalize": "true",
            }
            try:
                response = SESSION.get(
                    ENTITY_URL,
                    params=params,
                    headers=REQUEST_HEADERS,
                    timeout=REQUEST_TIMEOUT,
                )
                response.raise_for_status()
                data = response.json()
            except requests.RequestException as exc:
                print(f"WARNING: Failed to resolve title {title}: {exc}")
                data = {}

            entity_data = data.get("entities", {})
            entity = None
            qid = None
            for entity_id, item in entity_data.items():
                if item.get("missing") == "":
                    continue
                enwiki = item.get("sitelinks", {}).get("enwiki", {})
                if enwiki.get("title"):
                    entity = item
                    qid = entity_id
                    break

            if entity and qid:
                single_title_cache[title] = (qid, entity)
            else:
                single_title_cache[title] = (None, None)
            time.sleep(THROTTLE_SECONDS)
            result = single_title_cache[title]

        qid, entity = result
        keys = title_to_keys.get(title, [])
        if qid and entity:
            qid_to_entity[qid] = entity
            for key in keys:
                key_to_entity[key] = entity
                key_to_qid[key] = qid
        else:
            unmatched.update(keys)

    return key_to_entity, key_to_qid, qid_to_entity, unmatched


def extract_aliases(entity: dict) -> List[str]:
    aliases = []
    for alias in entity.get("aliases", {}).get("en", []):
        value = alias.get("value", "").strip()
        if value:
            aliases.append(value)
    return aliases


def extract_p1420_ids(entity: dict) -> List[str]:
    p1420_claims = entity.get("claims", {}).get("P1420", [])
    qids: List[str] = []
    for claim in p1420_claims:
        mainsnak = claim.get("mainsnak", {})
        if mainsnak.get("snaktype") != "value":
            continue
        datavalue = mainsnak.get("datavalue", {})
        value = datavalue.get("value")
        if isinstance(value, dict) and value.get("entity-type") == "item":
            qid = value.get("id")
            if qid:
                qids.append(qid)
    return qids


def english_label(entity: Optional[dict]) -> Optional[str]:
    if not entity:
        return None
    label = entity.get("labels", {}).get("en")
    if label:
        return label.get("value")
    return None


def dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    ordered = OrderedDict()
    for value in values:
        cleaned = value.strip()
        if not cleaned:
            continue
        if cleaned not in ordered:
            ordered[cleaned] = None
    return list(ordered.keys())


def main() -> None:
    if not INPUT_PATH.exists():
        raise SystemExit(f"Input CSV not found: {INPUT_PATH}")

    with INPUT_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV missing headers")

        index_field = reader.fieldnames[0]
        rows = list(reader)

    print(
        f"{PROGRESS_PREFIX} loaded {len(rows)} rows from {INPUT_PATH}",
        flush=True,
    )

    # Map names to candidate queries
    name_candidates: Dict[Tuple[str, str], List[str]] = {}
    for idx, row in enumerate(rows):
        scientific_raw = row.get("food_sci", "")
        common_raw = row.get("food_com", "")
        key = build_key(scientific_raw, common_raw, idx)
        queries = candidate_queries(scientific_raw, common_raw)
        if queries:
            name_candidates.setdefault(key, queries)

    title_to_keys: Dict[str, List[Tuple[str, str]]] = {}
    for key, queries in name_candidates.items():
        primary = queries[0]
        title = primary.replace(" ", "_")
        if title:
            title_to_keys.setdefault(title, []).append(key)

    print(
        f"{PROGRESS_PREFIX} prepared {len(name_candidates)} unique name keys "
        f"across {len(title_to_keys)} primary titles",
        flush=True,
    )

    key_to_entity, key_to_qid, qid_to_entity, unmatched_keys = fetch_entities_for_titles(
        title_to_keys
    )

    if unmatched_keys:
        print(
            f"{PROGRESS_PREFIX} falling back to search for {len(unmatched_keys)} keys",
            flush=True,
        )
        query_cache: Dict[str, Optional[str]] = {}
        fallback_qids: Dict[Tuple[str, str], str] = {}
        for key in unmatched_keys:
            queries = name_candidates.get(key, [])
            qid: Optional[str] = None
            for query in queries:
                if query not in query_cache:
                    query_cache[query] = search_entity_id(query)
                    time.sleep(THROTTLE_SECONDS)
                qid = query_cache[query]
                if qid:
                    break
            if qid:
                fallback_qids[key] = qid

        additional_entities = fetch_entities(fallback_qids.values())
        for key, qid in fallback_qids.items():
            entity = additional_entities.get(qid)
            if entity:
                key_to_entity[key] = entity
                key_to_qid[key] = qid
                qid_to_entity[qid] = entity
        unresolved_after_search = unmatched_keys - set(fallback_qids.keys())
        if unresolved_after_search:
            print(
                f"{PROGRESS_PREFIX} no Wikidata entity found for {len(unresolved_after_search)} keys",
                flush=True,
            )

    entity_data = qid_to_entity.copy()

    # Collect aliases and referenced synonym QIDs
    key_to_aliases: Dict[Tuple[str, str], List[str]] = {}
    key_to_synonym_qids: Dict[Tuple[str, str], List[str]] = {}
    referenced_qids: Set[str] = set()

    for key, entity in key_to_entity.items():
        aliases = extract_aliases(entity) if entity else []
        synonym_qids = extract_p1420_ids(entity) if entity else []

        key_to_aliases[key] = aliases
        key_to_synonym_qids[key] = synonym_qids
        referenced_qids.update(synonym_qids)

    # Fetch referenced entities for P1420 labels
    missing_qids = referenced_qids - set(entity_data.keys())
    if missing_qids:
        print(
            f"{PROGRESS_PREFIX} fetching {len(missing_qids)} synonym-linked entities",
            flush=True,
        )
        entity_data.update(fetch_entities(missing_qids))

    key_to_synonyms: Dict[Tuple[str, str], str] = {}
    for key in name_candidates.keys():
        aliases = key_to_aliases.get(key, [])
        synonym_qids = key_to_synonym_qids.get(key, [])

        synonym_labels = [
            english_label(entity_data.get(qid)) for qid in synonym_qids
        ]
        combined = aliases + [label for label in synonym_labels if label]

        # remove case-insensitive duplicates and the canonical name itself
        canonical_lower = key[1].lower()
        filtered: List[str] = []
        seen_lower: Set[str] = set()
        for item in combined:
            item_clean = item.strip()
            if not item_clean:
                continue
            item_lower = item_clean.lower()
            if item_lower == canonical_lower:
                continue
            if item_lower not in seen_lower:
                filtered.append(item_clean)
                seen_lower.add(item_lower)

        key_to_synonyms[key] = "; ".join(filtered)

    updated_rows: List[dict] = []
    for idx, row in enumerate(rows):
        scientific_raw = row.get("food_sci", "")
        common_raw = row.get("food_com", "")
        key = build_key(scientific_raw, common_raw, idx)
        synonyms_wiki = key_to_synonyms.get(key, "")

        updated_rows.append(
            {
                index_field: row.get(index_field, ""),
                "food_com": row.get("food_com", ""),
                "food_sci": row.get("food_sci", ""),
                OPEN_TREE_FIELD: row.get(OPEN_TREE_FIELD, ""),
                WIKI_FIELD: synonyms_wiki,
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
    print(f"{PROGRESS_PREFIX} wrote results to {OUTPUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
