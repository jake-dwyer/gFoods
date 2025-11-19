# gFoods Synonym Scrapers

This workspace holds a mock food taxonomy dataset (`ndm_foods.csv`) and utilities for enriching and cleaning it:

- `scrape/scrape_opentree_synonyms.py` — Open Tree of Life
- `scrape/scrape_wikisearch_synonyms.py` — Wikidata (wiki search)
- `scrape/scrape_ncbi_synonyms.py` — NCBI Taxonomy
- `dedupe_ndm_foods.py` — removes duplicate rows while preserving column order
- `build_synonym_overlap.py` — writes an all-sources overlap column
- `match_foods.py` — groups near-duplicate common names, aggregates scientific names, and can compress rows

Each script reads the CSV, fills its corresponding synonym column, and rewrites the file (unless you target a different output path).

---

## 1. Environment Setup

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` pins the non-stdlib dependencies (`requests` and `beautifulsoup4`) needed across the scrapers. Keep it updated if you add more packages later.

> **Tip:** All commands below assume you are inside the virtual environment you just created.

---

## 2. Data File

The canonical CSV is `ndm_foods.csv`. Columns:

1. anonymous row index (from the source dataset)
2. `food_com` — common name
3. `food_sci` — scientific/taxon name
4. `synonyms_open_tree_of_life`
5. `synonyms_wiki_search`
6. `synonyms_ncbi`

Run each scraper to populate its respective column. They are idempotent: re-running will overwrite the target column with freshly fetched data.

Before bulk runs consider backing up the CSV:

```bash
cp ndm_foods.csv ndm_foods.backup.csv
```

---

## 3. Synonym Overlap

`build_synonym_overlap.py` scans the three synonym columns and captures values shared across all of them in `synonyms_all_sources`.

```bash
.venv/bin/python build_synonym_overlap.py             # writes ndm_foods_with_overlap.csv
.venv/bin/python build_synonym_overlap.py --in-place  # overwrite the input file
```

Flags worth noting:

- `-o/--output` — choose a different target CSV.
- `--output-column` — change the destination column name.
- `--case-sensitive` — require exact casing instead of the default case-insensitive match.

The script reports how many rows contained overlaps and the total strings written so you can gauge coverage.

---

## 4. Open Tree of Life Synonyms

```bash
.venv/bin/python scrape/scrape_opentree_synonyms.py
```

- Batches scientific names (40 per request) against Open Tree of Life’s TNRS API.
- Populates `synonyms_open_tree_of_life` with a semicolon-delimited list.
- Prints progress for batches and throttles politely.

Approximate runtime: ~5 minutes for ~20k records.

---

## 5. Wikidata Synonyms

```bash
.venv/bin/python scrape/scrape_wikisearch_synonyms.py
```

- Resolves each entry to a Wikidata item via `wbgetentities` (falling back to `wbsearchentities`).
- Collects English aliases and taxon-synonym links (P1420 items) and stores them in `synonyms_wiki_search`.
- Emits detailed progress logs (batches, leftover titles, fallbacks) so you can track long runs.

This job is the slowest; expect 30–40 minutes. For smoke tests:

```bash
.venv/bin/python scrape/scrape_wikisearch_synonyms.py --limit 500 --output ndm_wiki_sample.csv
```

Use the `--limit` flag to process a subset and inspect results before committing to the full dataset. You can also direct output elsewhere while testing.

---

## 6. NCBI Synonyms

```bash
.venv/bin/python scrape/scrape_ncbi_synonyms.py
```

- Uses NCBI E-utilities (`esearch` + `efetch`) to find Taxonomy IDs and pull synonym/common-name lists.
- Populates `synonyms_ncbi` with unique values, semicolon-separated.
- Runs with a conservative delay (~0.34 s/request) to respect NCBI rate limits. With ~20k rows expect a multi-hour run.

Helpful flags:

```bash
.venv/bin/python scrape/scrape_ncbi_synonyms.py --limit 500 --output ndm_ncbi_sample.csv
```

`--limit` lets you test on a small slice; `--output` sets a different target CSV (default overwrites the input).

> **Note:** NCBI asks for a real contact; update `NCBI_EMAIL` inside the script before long runs.

---

## 7. Matching and Grouping Foods (`match_foods.py`)

This script cleans names (drops “raw”, normalises separators, strips accents/`?`), groups similar common names using exact + bigram fuzzy matching (length-aware thresholds: long names need ~0.9, short ~0.6; borderline pairs can go to ChatGPT), and writes a CSV with the same headers as the input.

Key behaviours (current branch):
- Drops exact duplicate rows.
- For exact common-name duplicates, keeps the row with the richest synonyms (`synonyms_all_sources` preferred, then total synonyms across sources, then presence of `food_sci`).
- Exact match first, then fuzzy match. Pepper guard avoids cross-pepper merges. Borderline pairs optionally sent to ChatGPT (`--use-openai`, needs `OPENAI_API_KEY` in `.env`).
- Aggregates common names per group but preserves the row’s scientific name; after grouping, prunes duplicate common names again to the richest synonym row.
- Output headers mirror the input (e.g., index, `food_com`, `food_sci`, `synonyms_*`, `synonyms_all_sources`).

Examples:
```bash
# Default (no OpenAI)
.venv/bin/python match_foods.py -i ndm_foods_with_overlap.csv

# With ChatGPT review of borderline pairs
.venv/bin/python match_foods.py -i ndm_foods_with_overlap.csv --use-openai --max-openai-pairs 50 --openai-batch-size 10

# Collapse to one row per merged group
.venv/bin/python match_foods.py -i ndm_foods_with_overlap.csv --collapse-groups
```

## 8. Full Workflow (this branch)

1. **Back up** `ndm_foods.csv`.
2. **Scrape synonyms** (any order, avoid parallel runs):
   - `.venv/bin/python scrape/scrape_opentree_synonyms.py`
   - `.venv/bin/python scrape/scrape_wikisearch_synonyms.py`
   - `.venv/bin/python scrape/scrape_ncbi_synonyms.py`
   Each overwrites only its target column in `ndm_foods.csv` by default; use `-o` to write elsewhere.
3. **Build overlap**: `.venv/bin/python build_synonym_overlap.py` (default output `ndm_foods_with_overlap.csv`; `--in-place` to overwrite).
4. **Match/clean**: `.venv/bin/python match_foods.py -i ndm_foods_with_overlap.csv [--use-openai --max-openai-pairs N --openai-batch-size M]`
   - Keeps the original column schema.
   - Prefers rows with richer synonyms for duplicate common names; blank or poorer rows are dropped.
   - Fuzzy thresholds as above; peppers guarded.
   - Outputs `<input>_matched.csv`.

---

## 9. Requirements Summary

- Python 3.10+
- `requests`
- `beautifulsoup4` (only for the wiki scraper)

No other runtime dependencies are required.

---

## 10. Troubleshooting

| Issue | Likely Cause | Fix |
|-------|--------------|-----|
| `requests.exceptions.HTTPError` from OpenTree/NCBI | API hiccup or throttle | rerun after a pause; ensure logs show polite spacing |
| Empty synonym column for an entry | Source lacks synonyms or name mismatch | inspect the source site manually; adjust `food_sci` if needed |
| Wikidata script slow | Large dataset; API throttling | use `--limit` for testing; run overnight for full dataset |
| NCBI API 429 errors | Too many requests too quickly | increase `THROTTLE_SECONDS` constant |

---
