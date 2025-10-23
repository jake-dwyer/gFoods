#!/usr/bin/env python3
"""
Utility script for removing duplicate rows from ndm_foods.csv.

By default the script keeps the first occurrence for each unique `food_com`
value and writes the deduplicated dataset to `ndm_foods_deduped.csv`.
Use `--columns` to change the subset of columns used as the deduplication key,
or `--in-place` to overwrite the input file.
"""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, Mapping, Sequence, Tuple


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove duplicate rows from ndm_foods.csv")
    parser.add_argument(
        "input_csv",
        nargs="?",
        default="ndm_foods.csv",
        help="Path to the source CSV (default: ndm_foods.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path for the deduplicated CSV (default: <input>_deduped.csv, ignored with --in-place)",
    )
    parser.add_argument(
        "-c",
        "--columns",
        default=["food_com"],
        nargs="+",
        help="Column names that form the deduplication key (default: food_com)",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite the input CSV with the deduplicated rows",
    )
    return parser.parse_args(argv)


def make_key(row: Mapping[str, str], columns: Iterable[str]) -> Tuple[str, ...]:
    return tuple(row.get(column, "") for column in columns)


def load_rows(path: Path) -> Tuple[Sequence[str], list[OrderedDict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"No header row found in {path}")
        rows = [OrderedDict(row) for row in reader]
    return fieldnames, rows


def dedupe_rows(
    rows: Iterable[OrderedDict[str, str]], columns: Sequence[str]
) -> Tuple[list[OrderedDict[str, str]], int]:
    seen: set[Tuple[str, ...]] = set()
    deduped: list[OrderedDict[str, str]] = []
    dropped = 0
    for row in rows:
        key = make_key(row, columns)
        if key in seen:
            dropped += 1
            continue
        seen.add(key)
        deduped.append(row)
    return deduped, dropped


def write_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[OrderedDict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    fieldnames, rows = load_rows(input_path)

    missing_columns = [col for col in args.columns if col not in fieldnames]
    if missing_columns:
        missing = ", ".join(missing_columns)
        available = ", ".join(fieldnames)
        raise ValueError(f"Column(s) not found: {missing} (available: {available})")

    deduped_rows, removed_count = dedupe_rows(rows, args.columns)

    if args.in_place:
        output_path = input_path
    else:
        default_output = input_path.with_name(f"{input_path.stem}_deduped.csv")
        output_path = Path(args.output) if args.output else default_output

    write_rows(output_path, fieldnames, deduped_rows)

    print(f"Input rows: {len(rows)}")
    print(f"Removed duplicates: {removed_count}")
    print(f"Output rows: {len(deduped_rows)}")
    print(f"Wrote deduplicated CSV to: {output_path}")


if __name__ == "__main__":
    main()
