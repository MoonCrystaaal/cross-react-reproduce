#!/usr/bin/env python3
"""Rewrite existing allergen ranking CSVs with competition ranks for score ties."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True, help="Ranking CSV or directory")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV or directory")
    return parser.parse_args(argv)


def discover(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    files = sorted(path.glob("*.allergen_rankings.csv"))
    if not files:
        raise RuntimeError(f"No per-query ranking CSVs found: {path}")
    return files


def score_value(row: dict[str, str]) -> float:
    column = "score" if "score" in row else "best_pcc"
    value = row.get(column, "").strip()
    return float("-inf") if not value else float(value)


def rerank_file(source: Path, destination: Path) -> tuple[int, int]:
    with source.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or ())
        rows = list(reader)
    required = {"rank", "query_allergen", "target_allergen"}
    missing = sorted(required - set(fieldnames))
    if missing:
        raise RuntimeError(f"{source}: missing columns {missing}")
    if "score" not in fieldnames and "best_pcc" not in fieldnames:
        raise RuntimeError(f"{source}: missing score or best_pcc column")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["query_allergen"]].append(row)
    output_rows: list[dict[str, str]] = []
    tie_count = 0
    for query in sorted(grouped):
        query_rows = sorted(
            grouped[query],
            key=lambda row: (-score_value(row), row["target_allergen"]),
        )
        previous_score: float | None = None
        current_rank = 0
        for position, row in enumerate(query_rows, start=1):
            score = score_value(row)
            if position == 1 or score != previous_score:
                current_rank = position
            else:
                tie_count += 1
            row["rank"] = str(current_rank)
            output_rows.append(row)
            previous_score = score

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(destination) + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    os.replace(temporary, destination)
    return len(output_rows), tie_count


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    files = discover(args.input)
    total_rows = 0
    total_ties = 0
    for source in files:
        destination = args.output if args.input.is_file() else args.output / source.name
        rows, ties = rerank_file(source, destination)
        total_rows += rows
        total_ties += ties
    print(f"files={len(files)} rows={total_rows} tied_rows={total_ties}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

