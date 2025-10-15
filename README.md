# gFoods Synonym Scrapers

This workspace holds a mock food taxonomy dataset (`ndm_foods.csv`) and three standalone Python scripts that populate synonym columns from different sources:

- `scrape_opentree_synonyms.py` — Open Tree of Life
- `scrape_wikisearch_synonyms.py` — Wikidata (wiki search)
- `scrape_ncbi_synonyms.py` — NCBI Taxonomy

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

## 3. Open Tree of Life Synonyms

```bash
.venv/bin/python scrape_opentree_synonyms.py
```

- Batches scientific names (40 per request) against Open Tree of Life’s TNRS API.
- Populates `synonyms_open_tree_of_life` with a semicolon-delimited list.
- Prints progress for batches and throttles politely.

Approximate runtime: ~5 minutes for ~20k records.

---

## 4. Wikidata Synonyms

```bash
.venv/bin/python scrape_wikisearch_synonyms.py
```

- Resolves each entry to a Wikidata item via `wbgetentities` (falling back to `wbsearchentities`).
- Collects English aliases and taxon-synonym links (P1420 items) and stores them in `synonyms_wiki_search`.
- Emits detailed progress logs (batches, leftover titles, fallbacks) so you can track long runs.

This job is the slowest; expect 30–40 minutes. For smoke tests:

```bash
.venv/bin/python scrape_wikisearch_synonyms.py --limit 500 --output ndm_wiki_sample.csv
```

Use the `--limit` flag to process a subset and inspect results before committing to the full dataset. You can also direct output elsewhere while testing.

---

## 5. NCBI Synonyms

```bash
.venv/bin/python scrape_ncbi_synonyms.py
```

- Uses NCBI E-utilities (`esearch` + `efetch`) to find Taxonomy IDs and pull synonym/common-name lists.
- Populates `synonyms_ncbi` with unique values, semicolon-separated.
- Runs with a conservative delay (~0.34 s/request) to respect NCBI rate limits. With ~20k rows expect a multi-hour run.

Helpful flags:

```bash
.venv/bin/python scrape_ncbi_synonyms.py --limit 500 --output ndm_ncbi_sample.csv
```

`--limit` lets you test on a small slice; `--output` sets a different target CSV (default overwrites the input).

> **Note:** NCBI asks for a real contact; update `NCBI_EMAIL` inside the script before long runs.

---

## 6. Workflow Suggestions

1. **Back up** `ndm_foods.csv`.
2. **Run the scrapers sequentially** (OpenTree → Wikidata → NCBI). Each script only touches its own column, so order is flexible but avoid parallel runs to prevent throttling issues.
3. **Version control**: commit after each successful script to isolate changes per source.
4. **Error handling**: scripts log warnings when a lookup fails; rerun those specific names or inspect the CSV to decide on manual fixes.

---

## 7. Requirements Summary

- Python 3.10+
- `requests`
- `beautifulsoup4` (only for the wiki scraper)

No other runtime dependencies are required.

---

## 8. Troubleshooting

| Issue | Likely Cause | Fix |
|-------|--------------|-----|
| `requests.exceptions.HTTPError` from OpenTree/NCBI | API hiccup or throttle | rerun after a pause; ensure logs show polite spacing |
| Empty synonym column for an entry | Source lacks synonyms or name mismatch | inspect the source site manually; adjust `food_sci` if needed |
| Wikidata script slow | Large dataset; API throttling | use `--limit` for testing; run overnight for full dataset |
| NCBI API 429 errors | Too many requests too quickly | increase `THROTTLE_SECONDS` constant |

---

## 9. TODO

- Scrape for missing taxon names
- Scrape for missing common names
- Clean data, remove duplicates
- Add N/A fields when no such scientific name exists / synonyms etc.
- Add another script to path over the enitre CSV, returning only synonyms that all 3 sources recognize.
