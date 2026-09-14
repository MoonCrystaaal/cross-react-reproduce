#!/usr/bin/env python3
"""Evaluate allergen rankings using expected metrics for score ties."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Sequence
from xml.etree import ElementTree as ET


VERSION = "ranking-benchmark-v3-tie-expected"
DEFAULT_KS = (1, 3, 5, 10, 20)
REQUIRED_RANKING_BASE_COLUMNS = {
    "rank",
    "query_allergen",
    "target_allergen",
}
PER_QUERY_COLUMNS = (
    "tool",
    "score_name",
    "query_allergen",
    "k",
    "relevant_target_count",
    "retrieved_relevant_count",
    "hit_at_k",
    "recall_at_k",
    "dcg_at_k",
    "idcg_at_k",
    "ndcg_at_k",
    "weighted_hits_at_k",
    "weighted_recall_at_k",
)
SUMMARY_COLUMNS = (
    "tool",
    "score_name",
    "k",
    "evaluated_queries",
    "queries_with_hit",
    "hits_at_k",
    "recall_at_k",
    "ndcg_at_k",
    "weighted_hits_at_k",
    "weighted_recall_at_k",
)


class BenchmarkError(RuntimeError):
    """Raised for an input or validation error."""


@dataclass(frozen=True)
class GroundTruth:
    relevant: dict[str, frozenset[str]]
    allergens: frozenset[str]
    row_count: int
    positive_pair_count: int


@dataclass(frozen=True)
class Ranking:
    query: str
    ranks: tuple[int, ...]
    targets: tuple[str, ...]
    scores: tuple[float, ...]
    source: Path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate per-query allergen ranking CSV files against the ordered_pairs "
            "sheet in pair_label_dataset.xlsx."
        )
    )
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--rankings-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--tool",
        choices=("cross-react", "surface-id"),
        required=True,
        help=(
            "Score semantics: cross-react uses score_name=pcc and max(0, pcc); "
            "surface-id uses score_name=frac_geo without transformation."
        ),
    )
    parser.add_argument("--ks", type=int, nargs="+", default=list(DEFAULT_KS))
    parser.add_argument("--sheet", default="ordered_pairs")
    parser.add_argument(
        "--strict",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fail on missing/extra queries or incomplete target rankings (default: true).",
    )
    return parser.parse_args(argv)


def normalize_name(value: object, *, context: str) -> str:
    name = str(value).strip()
    if not name:
        raise BenchmarkError(f"Empty allergen name: {context}")
    return name


def cell_column(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference.upper())
    if not match:
        raise BenchmarkError(f"Invalid XLSX cell reference: {reference}")
    number = 0
    for char in match.group(1):
        number = number * 26 + ord(char) - ord("A") + 1
    return number - 1


def read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    path = "xl/sharedStrings.xml"
    if path not in archive.namelist():
        return []
    root = ET.fromstring(archive.read(path))
    values = []
    for item in root.findall("{*}si"):
        values.append("".join(node.text or "" for node in item.findall(".//{*}t")))
    return values


def resolve_sheet_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationship_id = None
    for sheet in workbook.findall(".//{*}sheet"):
        if sheet.attrib.get("name") == sheet_name:
            relationship_id = sheet.attrib.get(
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            )
            break
    if relationship_id is None:
        available = [sheet.attrib.get("name", "") for sheet in workbook.findall(".//{*}sheet")]
        raise BenchmarkError(f"Sheet {sheet_name!r} not found; available={available}")

    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    for relationship in relationships.findall("{*}Relationship"):
        if relationship.attrib.get("Id") == relationship_id:
            target = relationship.attrib["Target"].replace("\\", "/")
            if target.startswith("/"):
                return target.lstrip("/")
            if target.startswith("xl/"):
                return target
            return "xl/" + target.lstrip("./")
    raise BenchmarkError(f"Relationship for sheet {sheet_name!r} was not found")


def xlsx_rows(path: Path, sheet_name: str) -> Iterator[list[object]]:
    try:
        with zipfile.ZipFile(path) as archive:
            shared_strings = read_shared_strings(archive)
            sheet_path = resolve_sheet_path(archive, sheet_name)
            with archive.open(sheet_path) as stream:
                for _, row in ET.iterparse(stream, events=("end",)):
                    if not row.tag.endswith("}row"):
                        continue
                    cells: dict[int, object] = {}
                    for cell in row.findall("{*}c"):
                        index = cell_column(cell.attrib.get("r", ""))
                        cell_type = cell.attrib.get("t")
                        value_node = cell.find("{*}v")
                        if cell_type == "inlineStr":
                            value: object = "".join(
                                node.text or "" for node in cell.findall(".//{*}t")
                            )
                        elif value_node is None:
                            value = ""
                        elif cell_type == "s":
                            value = shared_strings[int(value_node.text or "0")]
                        elif cell_type == "b":
                            value = value_node.text == "1"
                        else:
                            raw = value_node.text or ""
                            try:
                                value = float(raw) if "." in raw else int(raw)
                            except ValueError:
                                value = raw
                        cells[index] = value
                    width = max(cells, default=-1) + 1
                    yield [cells.get(index, "") for index in range(width)]
                    row.clear()
    except (KeyError, zipfile.BadZipFile, ET.ParseError, OSError) as exc:
        raise BenchmarkError(f"Could not read XLSX {path}: {exc}") from exc


def load_ground_truth(path: Path, sheet_name: str = "ordered_pairs") -> GroundTruth:
    if not path.is_file():
        raise BenchmarkError(f"Ground-truth file not found: {path}")
    rows = xlsx_rows(path, sheet_name)
    try:
        header = [str(value).strip() for value in next(rows)]
    except StopIteration as exc:
        raise BenchmarkError(f"Ground-truth sheet is empty: {sheet_name}") from exc

    required = ("allergen1", "allergen2", "label")
    missing = [column for column in required if column not in header]
    if missing:
        raise BenchmarkError(f"Missing ground-truth columns: {missing}")
    positions = {column: header.index(column) for column in required}

    relevant: dict[str, set[str]] = defaultdict(set)
    allergens: set[str] = set()
    seen_pairs: set[tuple[str, str]] = set()
    row_count = 0
    positive_pair_count = 0
    for row_number, row in enumerate(rows, start=2):
        if not any(str(value).strip() for value in row):
            continue
        row_count += 1
        try:
            left = normalize_name(row[positions["allergen1"]], context=f"row {row_number}")
            right = normalize_name(row[positions["allergen2"]], context=f"row {row_number}")
            raw_label = row[positions["label"]]
        except IndexError as exc:
            raise BenchmarkError(f"Incomplete ground-truth row {row_number}: {row}") from exc
        if left == right:
            raise BenchmarkError(f"Self-pair in ground truth at row {row_number}: {left}")
        pair = tuple(sorted((left, right)))
        if pair in seen_pairs:
            raise BenchmarkError(f"Duplicate unordered pair at row {row_number}: {pair}")
        seen_pairs.add(pair)
        allergens.update((left, right))
        try:
            label = int(raw_label)
        except (TypeError, ValueError) as exc:
            raise BenchmarkError(f"Invalid label at row {row_number}: {raw_label!r}") from exc
        if label not in (0, 1):
            raise BenchmarkError(f"Non-binary label at row {row_number}: {label}")
        if label == 1:
            positive_pair_count += 1
            relevant[left].add(right)
            relevant[right].add(left)

    if not relevant:
        raise BenchmarkError("Ground truth contains no positive pairs")
    return GroundTruth(
        relevant={query: frozenset(targets) for query, targets in relevant.items()},
        allergens=frozenset(allergens),
        row_count=row_count,
        positive_pair_count=positive_pair_count,
    )


def score_name_for_tool(tool: str) -> str:
    if tool == "cross-react":
        return "pcc"
    if tool == "surface-id":
        return "frac_geo"
    raise BenchmarkError(f"Unsupported tool: {tool}")


def prepare_score(score: float, tool: str, *, context: str) -> float:
    if not math.isfinite(score):
        raise BenchmarkError(f"{context}: score must be finite")
    if tool == "cross-react":
        if score < -1.000001 or score > 1.000001:
            raise BenchmarkError(f"{context}: PCC outside [-1, 1]: {score}")
        return max(0.0, min(1.0, score))
    if tool == "surface-id":
        if score < -0.000001 or score > 1.000001:
            raise BenchmarkError(f"{context}: frac_geo outside [0, 1]: {score}")
        return min(1.0, max(0.0, score))
    raise BenchmarkError(f"Unsupported tool: {tool}")


def load_ranking(path: Path, tool: str) -> Ranking:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            columns = set(reader.fieldnames or ())
            missing = sorted(REQUIRED_RANKING_BASE_COLUMNS - columns)
            if missing:
                raise BenchmarkError(f"{path.name}: missing columns {missing}")
            if "score" in columns:
                score_column = "score"
            elif "best_pcc" in columns:
                score_column = "best_pcc"
            else:
                raise BenchmarkError(
                    f"{path.name}: missing score column; expected 'score' "
                    "or legacy 'best_pcc'"
                )
            parsed: list[tuple[int, str, str, float, float]] = []
            for row_number, row in enumerate(reader, start=2):
                query = normalize_name(row["query_allergen"], context=f"{path.name}:{row_number}")
                target = normalize_name(row["target_allergen"], context=f"{path.name}:{row_number}")
                try:
                    rank = int(row["rank"])
                except (TypeError, ValueError) as exc:
                    raise BenchmarkError(
                        f"{path.name}:{row_number}: invalid rank {row['rank']!r}"
                    ) from exc
                if rank < 1:
                    raise BenchmarkError(f"{path.name}:{row_number}: rank must be positive")
                try:
                    score = float(row[score_column])
                except (TypeError, ValueError) as exc:
                    raise BenchmarkError(
                        f"{path.name}:{row_number}: invalid {score_column} "
                        f"{row[score_column]!r}"
                    ) from exc
                prepared_score = prepare_score(
                    score,
                    tool,
                    context=f"{path.name}:{row_number}",
                )
                if row.get("tool") and row["tool"].strip() != tool:
                    raise BenchmarkError(
                        f"{path.name}:{row_number}: tool={row['tool']!r} "
                        f"does not match --tool {tool!r}"
                    )
                expected_score_name = score_name_for_tool(tool)
                if row.get("score_name") and row["score_name"].strip() != expected_score_name:
                    raise BenchmarkError(
                        f"{path.name}:{row_number}: score_name={row['score_name']!r} "
                        f"does not match {expected_score_name!r}"
                    )
                parsed.append((rank, query, target, prepared_score, score))
    except OSError as exc:
        raise BenchmarkError(f"Could not read ranking {path}: {exc}") from exc

    if not parsed:
        raise BenchmarkError(f"Ranking file has no rows: {path}")
    queries = {row[1] for row in parsed}
    if len(queries) != 1:
        raise BenchmarkError(f"{path.name}: multiple query names {sorted(queries)}")
    parsed.sort(key=lambda row: (row[0], row[2]))
    ranks = [row[0] for row in parsed]
    expected_rank = 1
    offset = 0
    previous_raw_score: float | None = None
    while offset < len(parsed):
        rank = parsed[offset][0]
        if rank != expected_rank:
            raise BenchmarkError(
                f"{path.name}: rank {rank} at position {offset + 1}; "
                f"expected competition rank {expected_rank}"
            )
        end = offset + 1
        while end < len(parsed) and parsed[end][0] == rank:
            end += 1
        raw_scores = [row[4] for row in parsed[offset:end]]
        if any(score != raw_scores[0] for score in raw_scores[1:]):
            raise BenchmarkError(
                f"{path.name}: rank {rank} contains different scores"
            )
        if previous_raw_score is not None and raw_scores[0] >= previous_raw_score:
            raise BenchmarkError(
                f"{path.name}: scores are not strictly descending between tie groups"
            )
        previous_raw_score = raw_scores[0]
        expected_rank += end - offset
        offset = end
    targets = [row[2] for row in parsed]
    if len(targets) != len(set(targets)):
        duplicates = sorted({target for target in targets if targets.count(target) > 1})
        raise BenchmarkError(f"{path.name}: duplicate targets {duplicates}")
    scores = [row[3] for row in parsed]
    return Ranking(
        query=next(iter(queries)),
        ranks=tuple(ranks),
        targets=tuple(targets),
        scores=tuple(scores),
        source=path,
    )


def discover_rankings(rankings_dir: Path) -> list[Path]:
    if not rankings_dir.is_dir():
        raise BenchmarkError(f"Rankings directory not found: {rankings_dir}")
    files = sorted(rankings_dir.glob("*.allergen_rankings.csv"))
    if not files:
        raise BenchmarkError(f"No per-query ranking files found in {rankings_dir}")
    return files


def collect_rankings(files: Iterable[Path], tool: str) -> dict[str, Ranking]:
    rankings: dict[str, Ranking] = {}
    for path in files:
        ranking = load_ranking(path, tool)
        if ranking.query in rankings:
            raise BenchmarkError(
                f"Multiple ranking files for {ranking.query!r}: "
                f"{rankings[ranking.query].source.name}, {path.name}"
            )
        rankings[ranking.query] = ranking
    return rankings


def validate_rankings(
    truth: GroundTruth, rankings: dict[str, Ranking]
) -> tuple[dict[str, object], list[str]]:
    expected_queries = set(truth.relevant)
    found_queries = set(rankings)
    errors: list[str] = []
    details: dict[str, object] = {
        "expected_queries": len(expected_queries),
        "ranking_files_found": len(rankings),
        "missing_queries": sorted(expected_queries - found_queries),
        "unexpected_queries": sorted(found_queries - expected_queries),
        "query_errors": {},
    }
    if details["missing_queries"]:
        errors.append("Missing expected query ranking files")
    if details["unexpected_queries"]:
        errors.append("Found ranking files for queries without positive ground truth")

    query_errors: dict[str, list[str]] = {}
    for query in sorted(expected_queries & found_queries):
        ranking = rankings[query]
        problems: list[str] = []
        target_set = set(ranking.targets)
        expected_targets = set(truth.allergens) - {query}
        if query in target_set:
            problems.append("ranking contains the query itself")
        missing_candidates = sorted(expected_targets - target_set)
        unexpected_candidates = sorted(target_set - expected_targets)
        if missing_candidates:
            problems.append(f"missing candidate targets: {missing_candidates}")
        if unexpected_candidates:
            problems.append(f"unexpected candidate targets: {unexpected_candidates}")
        missing_relevant = sorted(set(truth.relevant[query]) - target_set)
        if missing_relevant:
            problems.append(f"missing relevant targets: {missing_relevant}")
        if problems:
            query_errors[query] = problems
            errors.append(f"Invalid ranking for {query}")
    details["query_errors"] = query_errors
    details["validation_errors"] = errors
    return details, errors


def compute_query_rows(
    tool: str,
    score_name: str,
    query: str,
    ranks: Sequence[int],
    targets: Sequence[str],
    scores: Sequence[float],
    relevant: frozenset[str],
    ks: Sequence[int],
) -> list[dict[str, object]]:
    if not (len(ranks) == len(targets) == len(scores)):
        raise BenchmarkError(f"Rank/target/score length mismatch for {query}")
    rows: list[dict[str, object]] = []
    for k in ks:
        retrieved = 0.0
        weighted_relevant = 0.0
        dcg = 0.0
        probability_no_hit = 1.0
        weighted_hits = 0.0
        offset = 0
        while offset < len(targets):
            rank = ranks[offset]
            end = offset + 1
            while end < len(targets) and ranks[end] == rank:
                end += 1
            group_size = end - offset
            slots = min(group_size, max(0, k - rank + 1))
            if slots <= 0:
                break
            group_relevance = [
                1 if target in relevant else 0 for target in targets[offset:end]
            ]
            relevant_in_group = sum(group_relevance)
            expected_group_relevant = slots * relevant_in_group / group_size
            retrieved += expected_group_relevant
            weighted_relevant += (slots / group_size) * sum(
                score
                for is_relevant, score in zip(
                    group_relevance, scores[offset:end]
                )
                if is_relevant
            )
            expected_relevance_per_position = relevant_in_group / group_size
            dcg += expected_relevance_per_position * sum(
                1.0 / math.log2(position + 1)
                for position in range(rank, rank + slots)
            )
            if relevant_in_group > 0:
                if slots == group_size:
                    group_hit_probability = 1.0
                else:
                    group_hit_probability = 1.0 - (
                        math.comb(group_size - relevant_in_group, slots)
                        / math.comb(group_size, slots)
                    )
                group_score = max(
                    score
                    for is_relevant, score in zip(
                        group_relevance, scores[offset:end]
                    )
                    if is_relevant
                )
                weighted_hits += probability_no_hit * group_hit_probability * group_score
                probability_no_hit *= 1.0 - group_hit_probability
            offset = end
        hit = 1.0 - probability_no_hit
        recall = retrieved / len(relevant)
        idcg = sum(
            1.0 / math.log2(rank + 1)
            for rank in range(1, min(k, len(relevant)) + 1)
        )
        ndcg = dcg / idcg
        weighted_recall = weighted_relevant / len(relevant)
        rows.append(
            {
                "tool": tool,
                "score_name": score_name,
                "query_allergen": query,
                "k": k,
                "relevant_target_count": len(relevant),
                "retrieved_relevant_count": retrieved,
                "hit_at_k": hit,
                "recall_at_k": recall,
                "dcg_at_k": dcg,
                "idcg_at_k": idcg,
                "ndcg_at_k": ndcg,
                "weighted_hits_at_k": weighted_hits,
                "weighted_recall_at_k": weighted_recall,
            }
        )
    return rows


def aggregate_rows(
    per_query: Sequence[dict[str, object]],
    ks: Sequence[int],
    tool: str,
    score_name: str,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for k in ks:
        rows = [row for row in per_query if row["k"] == k]
        count = len(rows)
        if count == 0:
            continue
        hits = sum(float(row["hit_at_k"]) for row in rows)
        output.append(
            {
                "tool": tool,
                "score_name": score_name,
                "k": k,
                "evaluated_queries": count,
                "queries_with_hit": hits,
                "hits_at_k": hits / count,
                "recall_at_k": sum(float(row["recall_at_k"]) for row in rows) / count,
                "ndcg_at_k": sum(float(row["ndcg_at_k"]) for row in rows) / count,
                "weighted_hits_at_k": sum(
                    float(row["weighted_hits_at_k"]) for row in rows
                )
                / count,
                "weighted_recall_at_k": sum(
                    float(row["weighted_recall_at_k"]) for row in rows
                )
                / count,
            }
        )
    return output


def format_csv_value(value: object) -> object:
    return f"{value:.10f}" if isinstance(value, float) else value


def write_csv_atomic(path: Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    temporary = Path(str(path) + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_csv_value(value) for key, value in row.items()})
    os.replace(temporary, path)


def write_json_atomic(path: Path, payload: object) -> None:
    temporary = Path(str(path) + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_ks(values: Sequence[int]) -> tuple[int, ...]:
    if any(value < 1 for value in values):
        raise BenchmarkError("All K values must be positive integers")
    return tuple(sorted(set(values)))


def run(args: argparse.Namespace) -> int:
    ks = normalize_ks(args.ks)
    score_name = score_name_for_tool(args.tool)
    truth = load_ground_truth(args.ground_truth, args.sheet)
    files = discover_rankings(args.rankings_dir)
    rankings = collect_rankings(files, args.tool)
    validation, errors = validate_rankings(truth, rankings)
    validation.update(
        {
            "benchmark_version": VERSION,
            "tool": args.tool,
            "score_name": score_name,
            "ground_truth_rows": truth.row_count,
            "positive_pairs": truth.positive_pair_count,
            "allergen_count": len(truth.allergens),
        }
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json_atomic(args.output_dir / "validation_report.json", validation)
    if errors and args.strict:
        raise BenchmarkError(
            f"Strict validation failed with {len(errors)} error(s); "
            f"see {args.output_dir / 'validation_report.json'}"
        )

    evaluable = sorted(set(truth.relevant) & set(rankings))
    if not args.strict:
        invalid = set(validation["query_errors"])
        evaluable = [query for query in evaluable if query not in invalid]
    if not evaluable:
        raise BenchmarkError("No valid queries are available for evaluation")

    per_query: list[dict[str, object]] = []
    for query in evaluable:
        ranking = rankings[query]
        per_query.extend(
            compute_query_rows(
                args.tool,
                score_name,
                query,
                ranking.ranks,
                ranking.targets,
                ranking.scores,
                truth.relevant[query],
                ks,
            )
        )
    summary = aggregate_rows(per_query, ks, args.tool, score_name)
    write_csv_atomic(args.output_dir / "per_query_metrics.csv", PER_QUERY_COLUMNS, per_query)
    write_csv_atomic(args.output_dir / "summary_metrics.csv", SUMMARY_COLUMNS, summary)
    configuration = {
        "benchmark_version": VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool": args.tool,
        "score_name": score_name,
        "ground_truth": str(args.ground_truth.resolve()),
        "ground_truth_sha256": sha256_file(args.ground_truth),
        "rankings_dir": str(args.rankings_dir.resolve()),
        "ranking_file_count": len(files),
        "output_dir": str(args.output_dir.resolve()),
        "sheet": args.sheet,
        "ks": list(ks),
        "strict": args.strict,
        "tie_policy": "expected metrics over all permutations within equal-score groups",
        "evaluated_queries": len(evaluable),
        "metric_definitions": {
            "hits_at_k": "macro mean expected hit probability at K over permutations within score ties",
            "recall_at_k": "macro mean expected recall at K over permutations within score ties",
            "ndcg_at_k": "macro mean expected binary-relevance NDCG at K over permutations within score ties",
            "weighted_hits_at_k": (
                "macro mean of the maximum tool-specific score among relevant targets in top K; "
                "0 when no relevant target is retrieved"
            ),
            "weighted_recall_at_k": (
                "macro mean of the sum of tool-specific scores for relevant targets in top K "
                "divided by the query's total relevant-target count"
            ),
        },
    }
    write_json_atomic(args.output_dir / "benchmark_configuration.json", configuration)

    print(f"Evaluated queries: {len(evaluable)}")
    print(f"Summary: {args.output_dir / 'summary_metrics.csv'}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except BenchmarkError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
