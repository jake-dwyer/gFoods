#!/usr/bin/env python3
"""
Generate an overlap column for ndm_foods.csv by intersecting synonym sources.

The script reads the CSV, finds synonym strings shared across all specified
columns, and writes them to a new semicolon-delimited column named
`synonyms_all_sources` by default. Results are written to `<input>_with_overlap.csv`
unless `--in-place` is provided.
"""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

DEFAULT_COLUMNS = [
    "synonyms_open_tree_of_life",
    "synonyms_wiki_search",
    "synonyms_ncbi",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Intersect synonym columns and store the shared values in a new column."
    )
    parser.add_argument(
        "input_csv",
        nargs="?",
        default="ndm_foods.csv",
        help="Path to the source CSV (default: ndm_foods.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path for the updated CSV (default: <input>_with_overlap.csv unless --in-place)",
    )
    parser.add_argument(
        "-c",
        "--columns",
        nargs="+",
        default=DEFAULT_COLUMNS,
        help="Synonym columns to intersect (default: OpenTree, Wikidata, NCBI)",
    )
    parser.add_argument(
        "--output-column",
        default="synonyms_all_sources",
        help="Name of the column to store the shared synonyms (default: synonyms_all_sources)",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the input CSV instead of writing to <input>_with_overlap.csv.",
    )
    parser.add_argument(
        "--case-sensitive",
        action="store_true",
        help="Treat synonyms as case-sensitive when computing the overlap (default: case-insensitive).",
    )
    return parser.parse_args(argv)


def normalise(value: str, *, case_sensitive: bool) -> str:
    value = value.strip()
    return value if case_sensitive else value.casefold()


def split_synonyms(value: str, *, case_sensitive: bool) -> Tuple[Dict[str, str], set[str]]:
    mapping: Dict[str, str] = {}
    synonyms = set()
    for raw in value.split(";"):
        item = raw.strip()
        if not item:
            continue
        key = normalise(item, case_sensitive=case_sensitive)
        if key not in mapping:  # preserve first spelling we encounter
            mapping[key] = item
        synonyms.add(key)
    return mapping, synonyms


def compute_overlap(
    row: Mapping[str, str],
    columns: Iterable[str],
    *,
    case_sensitive: bool,
) -> Tuple[list[str], int]:
    mappings: list[Dict[str, str]] = []
    synonym_sets: list[set[str]] = []

    for column in columns:
        values = row.get(column, "") or ""
        mapping, synonyms = split_synonyms(values, case_sensitive=case_sensitive)
        if not synonyms:
            return [], 0  # One column missing values means no intersection
        mappings.append(mapping)
        synonym_sets.append(synonyms)

    shared = set.intersection(*synonym_sets) if synonym_sets else set()
    if not shared:
        return [], 0

    ordered = OrderedDict()
    for synonym_key in sorted(shared):
        for mapping in mappings:
            if synonym_key in mapping:
                ordered[synonym_key] = mapping[synonym_key]
                break
    return list(ordered.values()), len(shared)


def process_rows(
    rows: Iterable[OrderedDict[str, str]],
    columns: Sequence[str],
    output_column: str,
    *,
    case_sensitive: bool,
) -> Tuple[list[OrderedDict[str, str]], int, int]:
    updated_rows: list[OrderedDict[str, str]] = []
    rows_with_overlap = 0
    total_synonyms = 0

    for row in rows:
        overlap_values, shared_size = compute_overlap(
            row, columns, case_sensitive=case_sensitive
        )
        row[output_column] = "; ".join(overlap_values)
        if overlap_values:
            rows_with_overlap += 1
        total_synonyms += shared_size
        updated_rows.append(row)

    return updated_rows, rows_with_overlap, total_synonyms


def load_csv(path: Path) -> Tuple[Sequence[str], list[OrderedDict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"No header row found in {path}")
        rows = [OrderedDict(row) for row in reader]
    return fieldnames, rows


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    input_path = Path(args.input_csv)

    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    fieldnames, rows = load_csv(input_path)

    missing_columns = [column for column in args.columns if column not in fieldnames]
    if missing_columns:
        raise ValueError(f"Column(s) not found: {', '.join(missing_columns)}")

    if args.output_column not in fieldnames:
        fieldnames = list(fieldnames) + [args.output_column]

    updated_rows, rows_with_overlap, total_synonyms = process_rows(
        rows,
        args.columns,
        args.output_column,
        case_sensitive=args.case_sensitive,
    )

    if args.in_place:
        output_path = Path(args.output) if args.output else input_path
    else:
        default_output = input_path.with_name(f"{input_path.stem}_with_overlap.csv")
        output_path = Path(args.output) if args.output else default_output

    write_csv(output_path, fieldnames, updated_rows)

    print(f"Rows processed: {len(updated_rows)}")
    print(f"Rows with overlap: {rows_with_overlap}")
    print(f"Total synonyms written: {total_synonyms}")
    print(f"Output file: {output_path}")


if __name__ == "__main__":
    main()
