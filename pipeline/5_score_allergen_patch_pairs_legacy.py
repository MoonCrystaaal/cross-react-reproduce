#!/usr/bin/env python3
"""Score selected query allergens by exhaustive patch-to-patch legacy PCC."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

from legacy_pcc import (
    MATRIX_ELEMENT_COUNT,
    MATRIX_SIZE,
    build_legacy_stats,
    legacy_crossreact_scores_from_stats,
)


AA_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]
MATRIX_COLUMNS = [f"m_{left}_{right}" for left in AA_ORDER for right in AA_ORDER]
METADATA_COLUMNS = [
    "patch_id",
    "source_pdb",
    "accession",
    "allergen",
    "patch_size",
    "connected_pair_count",
    "matrix_total_sum",
    "unknown_pair_count",
]
PATCH_HIT_COLUMNS = [
    "query_allergen",
    "query_accession",
    "query_source_pdb",
    "query_patch_id",
    "target_allergen",
    "target_accession",
    "target_source_pdb",
    "target_patch_id",
    "pcc",
]
PAIR_SCORE_COLUMNS = [
    "query_allergen",
    "target_allergen",
    "best_pcc",
    "query_patch_id",
    "target_patch_id",
    "query_source_pdb",
    "target_source_pdb",
    "query_accession",
    "target_accession",
    "query_patch_count",
    "target_patch_count",
    "patch_pairs_tested",
    "valid_pcc_count",
    "nan_pcc_count",
    "patch_hits",
]
RANKING_COLUMNS = ["rank"] + PAIR_SCORE_COLUMNS
SUMMARY_COLUMNS = [
    "query_allergen",
    "query_patch_count",
    "target_allergen_count",
    "patch_pairs_tested",
    "valid_pcc_count",
    "nan_pcc_count",
    "patch_hits",
]
ALGORITHM_VERSION = "legacy-patch-pairs-v2-competition-rank"


@dataclass(frozen=True)
class PatchDataset:
    patch_ids: np.ndarray
    source_pdbs: np.ndarray
    accessions: np.ndarray
    allergens: np.ndarray
    row_sums: np.ndarray
    sums: np.ndarray
    sumsq: np.ndarray
    indices_by_allergen: Dict[str, np.ndarray]

    @property
    def patch_count(self) -> int:
        return int(self.patch_ids.shape[0])


@dataclass(frozen=True)
class QueryPaths:
    hit_file: Path
    pair_file: Path
    ranking_file: Path
    checkpoint_file: Path


def parse_args() -> argparse.Namespace:
    run_root = Path.cwd().resolve()
    parser = argparse.ArgumentParser(
        description=(
            "Compare every patch of selected query allergens against every target "
            "patch with the legacy Cross-React PCC formula."
        )
    )
    parser.add_argument(
        "--flat-dir",
        type=Path,
        default=run_root / "3_connectivity" / "connectivity_matrix_flat",
        help="Directory containing stage 3 *.connectivity_flat.csv files",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=run_root / "5_allergen_pcc",
        help="Root directory for patch hits, pair scores, rankings, and logs",
    )
    parser.add_argument(
        "--query-list",
        type=Path,
        required=True,
        help="Text file containing one exact allergen name per line",
    )
    parser.add_argument(
        "--patch-hit-threshold",
        type=float,
        required=True,
        help="Store every finite patch pair whose legacy PCC is >= this value",
    )
    parser.add_argument(
        "--exclude-self",
        action="store_true",
        help="Exclude targets with the same allergen name as the query",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=256,
        help="Maximum number of query and target patches per score block (default: 256)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of query allergens processed concurrently with threads (default: 1)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse query outputs whose checkpoint matches the current inputs and options",
    )
    parser.add_argument(
        "--target-match",
        default="*.connectivity_flat.csv",
        help="Glob for stage 3 flat files (default: *.connectivity_flat.csv)",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "allergen_patch_pcc.log", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )


def sanitize_stem(name: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return sanitized or "unnamed"


def read_query_list(path: Path) -> List[str]:
    queries: List[str] = []
    seen = set()
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            if value in seen:
                logging.warning(
                    "Duplicate query allergen ignored at line %d: %s", line_number, value
                )
                continue
            seen.add(value)
            queries.append(value)
    if not queries:
        raise ValueError(f"Query list contains no allergen names: {path}")
    return queries


def iter_flat_files(flat_dir: Path, match: str) -> List[Path]:
    return sorted(flat_dir.glob(match), key=lambda item: item.name)


def load_patch_dataset(flat_dir: Path, match: str) -> tuple[PatchDataset, List[Path]]:
    files = iter_flat_files(flat_dir, match)
    if not files:
        raise ValueError(f"No target flat matrix files found in {flat_dir}")

    patch_ids: List[str] = []
    source_pdbs: List[str] = []
    accessions: List[str] = []
    allergens: List[str] = []
    legacy_vectors: List[Sequence[float]] = []
    sums: List[float] = []
    sumsq: List[float] = []

    for path in files:
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or [])
            missing = [
                column
                for column in METADATA_COLUMNS + MATRIX_COLUMNS
                if column not in fields
            ]
            if missing:
                preview = ", ".join(missing[:5])
                raise ValueError(f"Missing columns in {path.name}: {preview}")

            for line_number, row in enumerate(reader, start=2):
                allergen = row["allergen"].strip()
                patch_id = row["patch_id"].strip()
                if not allergen or not patch_id:
                    raise ValueError(
                        f"Blank allergen or patch_id in {path.name}, line {line_number}"
                    )
                try:
                    flat_values = np.fromiter(
                        (float(row[column]) for column in MATRIX_COLUMNS),
                        dtype=np.float64,
                        count=MATRIX_ELEMENT_COUNT,
                    )
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid matrix value in {path.name}, line {line_number}"
                    ) from exc
                matrix = flat_values.reshape(MATRIX_SIZE, MATRIX_SIZE)
                row_sums, matrix_sum, matrix_sumsq = build_legacy_stats(matrix)
                patch_ids.append(patch_id)
                source_pdbs.append(row["source_pdb"].strip())
                accessions.append(row["accession"].strip())
                allergens.append(allergen)
                legacy_vectors.append(row_sums)
                sums.append(matrix_sum)
                sumsq.append(matrix_sumsq)

    if not patch_ids:
        raise ValueError(f"No patch rows found in {flat_dir}")

    allergen_array = np.asarray(allergens, dtype=object)
    indices_by_allergen = {
        allergen: np.flatnonzero(allergen_array == allergen)
        for allergen in sorted(set(allergens))
    }
    dataset = PatchDataset(
        patch_ids=np.asarray(patch_ids, dtype=object),
        source_pdbs=np.asarray(source_pdbs, dtype=object),
        accessions=np.asarray(accessions, dtype=object),
        allergens=allergen_array,
        row_sums=np.asarray(legacy_vectors, dtype=np.float64),
        sums=np.asarray(sums, dtype=np.float64),
        sumsq=np.asarray(sumsq, dtype=np.float64),
        indices_by_allergen=indices_by_allergen,
    )
    return dataset, files


def build_input_signature(
    files: Sequence[Path], threshold: float, exclude_self: bool
) -> str:
    digest = hashlib.sha256()
    digest.update(ALGORITHM_VERSION.encode("ascii"))
    digest.update(f"|threshold={threshold:.17g}|exclude_self={exclude_self}".encode("ascii"))
    for path in files:
        stat = path.stat()
        digest.update(
            f"|{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}".encode(
                "utf-8", errors="surrogatepass"
            )
        )
    return digest.hexdigest()


def make_query_paths(output_root: Path, query: str) -> QueryPaths:
    key = sanitize_stem(query)
    return QueryPaths(
        hit_file=output_root / "patch_hits" / f"{key}.patch_hits.csv.gz",
        pair_file=(
            output_root
            / "pair_scores"
            / "by_query"
            / f"{key}.allergen_pair_best_scores.csv"
        ),
        ranking_file=(
            output_root
            / "rankings"
            / "by_query"
            / f"{key}.allergen_rankings.csv"
        ),
        checkpoint_file=output_root / "checkpoints" / f"{key}.done.json",
    )


def atomic_path(final_path: Path) -> Path:
    final_path.parent.mkdir(parents=True, exist_ok=True)
    return Path(str(final_path) + ".tmp")


def write_csv_atomic(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    temporary = atomic_path(path)
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_json_atomic(path: Path, payload: dict) -> None:
    temporary = atomic_path(path)
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def query_is_complete(paths: QueryPaths, query: str, signature: str) -> bool:
    required = [paths.hit_file, paths.pair_file, paths.ranking_file, paths.checkpoint_file]
    if not all(path.exists() for path in required):
        return False
    try:
        with paths.checkpoint_file.open("r", encoding="utf-8") as handle:
            checkpoint = json.load(handle)
    except (OSError, ValueError, TypeError):
        return False
    return checkpoint.get("query_allergen") == query and checkpoint.get("signature") == signature


def format_score(value: float) -> str:
    return "" if not math.isfinite(value) else f"{value:.12f}"


def best_candidate_from_block(
    scores: np.ndarray,
    query_indices: np.ndarray,
    target_indices: np.ndarray,
    dataset: PatchDataset,
) -> tuple[float, int, int] | None:
    finite = np.isfinite(scores)
    if not finite.any():
        return None
    best_score = float(np.max(scores[finite]))
    positions = np.argwhere(finite & (scores == best_score))
    candidates = [
        (
            str(dataset.patch_ids[int(query_indices[row])]),
            str(dataset.patch_ids[int(target_indices[column])]),
            int(query_indices[row]),
            int(target_indices[column]),
        )
        for row, column in positions
    ]
    _, _, query_index, target_index = min(candidates)
    return best_score, query_index, target_index


def score_one_query(
    query: str,
    dataset: PatchDataset,
    output_root: Path,
    threshold: float,
    exclude_self: bool,
    block_size: int,
    signature: str,
    resume: bool,
) -> QueryPaths:
    paths = make_query_paths(output_root, query)
    if resume and query_is_complete(paths, query, signature):
        logging.info("Resume: reusing completed query %s", query)
        return paths

    started = time.monotonic()
    query_indices = dataset.indices_by_allergen[query]
    target_allergens = sorted(dataset.indices_by_allergen)
    if exclude_self:
        target_allergens = [name for name in target_allergens if name != query]

    hit_temporary = atomic_path(paths.hit_file)
    pair_rows: List[dict] = []
    try:
        with gzip.open(hit_temporary, "wt", newline="", encoding="utf-8") as handle:
            hit_writer = csv.DictWriter(handle, fieldnames=PATCH_HIT_COLUMNS)
            hit_writer.writeheader()

            for target in target_allergens:
                target_indices = dataset.indices_by_allergen[target]
                patch_pairs_tested = int(query_indices.size * target_indices.size)
                valid_count = 0
                hit_count = 0
                best_score = float("nan")
                best_query_index: int | None = None
                best_target_index: int | None = None

                for query_start in range(0, query_indices.size, block_size):
                    query_block = query_indices[query_start : query_start + block_size]
                    for target_start in range(0, target_indices.size, block_size):
                        target_block = target_indices[
                            target_start : target_start + block_size
                        ]
                        scores = legacy_crossreact_scores_from_stats(
                            dataset.row_sums[query_block],
                            dataset.sums[query_block],
                            dataset.sumsq[query_block],
                            dataset.row_sums[target_block],
                            dataset.sums[target_block],
                            dataset.sumsq[target_block],
                        )
                        finite = np.isfinite(scores)
                        valid_count += int(np.count_nonzero(finite))

                        candidate = best_candidate_from_block(
                            scores, query_block, target_block, dataset
                        )
                        if candidate is not None:
                            candidate_score, candidate_query, candidate_target = candidate
                            candidate_ids = (
                                str(dataset.patch_ids[candidate_query]),
                                str(dataset.patch_ids[candidate_target]),
                            )
                            current_ids = (
                                str(dataset.patch_ids[best_query_index]),
                                str(dataset.patch_ids[best_target_index]),
                            ) if best_query_index is not None and best_target_index is not None else None
                            if (
                                not math.isfinite(best_score)
                                or candidate_score > best_score
                                or (candidate_score == best_score and candidate_ids < current_ids)
                            ):
                                best_score = candidate_score
                                best_query_index = candidate_query
                                best_target_index = candidate_target

                        hit_positions = np.argwhere(finite & (scores >= threshold))
                        hit_count += int(hit_positions.shape[0])
                        hit_writer.writerows(
                            {
                                "query_allergen": query,
                                "query_accession": str(dataset.accessions[int(query_block[row])]),
                                "query_source_pdb": str(dataset.source_pdbs[int(query_block[row])]),
                                "query_patch_id": str(dataset.patch_ids[int(query_block[row])]),
                                "target_allergen": target,
                                "target_accession": str(dataset.accessions[int(target_block[column])]),
                                "target_source_pdb": str(dataset.source_pdbs[int(target_block[column])]),
                                "target_patch_id": str(dataset.patch_ids[int(target_block[column])]),
                                "pcc": f"{float(scores[row, column]):.12f}",
                            }
                            for row, column in hit_positions
                        )

                pair_rows.append(
                    {
                        "query_allergen": query,
                        "target_allergen": target,
                        "best_pcc": format_score(best_score),
                        "query_patch_id": "" if best_query_index is None else str(dataset.patch_ids[best_query_index]),
                        "target_patch_id": "" if best_target_index is None else str(dataset.patch_ids[best_target_index]),
                        "query_source_pdb": "" if best_query_index is None else str(dataset.source_pdbs[best_query_index]),
                        "target_source_pdb": "" if best_target_index is None else str(dataset.source_pdbs[best_target_index]),
                        "query_accession": "" if best_query_index is None else str(dataset.accessions[best_query_index]),
                        "target_accession": "" if best_target_index is None else str(dataset.accessions[best_target_index]),
                        "query_patch_count": int(query_indices.size),
                        "target_patch_count": int(target_indices.size),
                        "patch_pairs_tested": patch_pairs_tested,
                        "valid_pcc_count": valid_count,
                        "nan_pcc_count": patch_pairs_tested - valid_count,
                        "patch_hits": hit_count,
                    }
                )

        os.replace(hit_temporary, paths.hit_file)
    finally:
        if hit_temporary.exists():
            hit_temporary.unlink()

    ranked_rows = sorted(
        pair_rows,
        key=lambda row: (
            row["best_pcc"] == "",
            -float(row["best_pcc"]) if row["best_pcc"] != "" else 0.0,
            row["target_allergen"],
            row["query_patch_id"],
            row["target_patch_id"],
        ),
    )
    ranking_rows = []
    previous_score: str | None = None
    current_rank = 0
    for position, row in enumerate(ranked_rows, start=1):
        score = row["best_pcc"]
        if position == 1 or score != previous_score:
            current_rank = position
        ranking_rows.append(dict(row, rank=current_rank))
        previous_score = score
    write_csv_atomic(paths.pair_file, PAIR_SCORE_COLUMNS, pair_rows)
    write_csv_atomic(paths.ranking_file, RANKING_COLUMNS, ranking_rows)
    write_json_atomic(
        paths.checkpoint_file,
        {
            "algorithm_version": ALGORITHM_VERSION,
            "query_allergen": query,
            "signature": signature,
            "query_patch_count": int(query_indices.size),
            "target_allergen_count": len(target_allergens),
        },
    )
    logging.info(
        "Completed query %s: query_patches=%d targets=%d elapsed=%.1fs",
        query,
        query_indices.size,
        len(target_allergens),
        time.monotonic() - started,
    )
    return paths


def combine_query_csvs(
    output_path: Path,
    fieldnames: Sequence[str],
    input_paths: Sequence[Path],
) -> None:
    def rows() -> Iterable[dict]:
        for input_path in input_paths:
            with input_path.open("r", newline="", encoding="utf-8") as handle:
                yield from csv.DictReader(handle)

    write_csv_atomic(output_path, fieldnames, rows())


def build_summary(query_paths: Sequence[QueryPaths]) -> List[dict]:
    summaries: List[dict] = []
    for paths in query_paths:
        with paths.pair_file.open("r", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        summaries.append(
            {
                "query_allergen": rows[0]["query_allergen"],
                "query_patch_count": rows[0]["query_patch_count"],
                "target_allergen_count": len(rows),
                "patch_pairs_tested": sum(int(row["patch_pairs_tested"]) for row in rows),
                "valid_pcc_count": sum(int(row["valid_pcc_count"]) for row in rows),
                "nan_pcc_count": sum(int(row["nan_pcc_count"]) for row in rows),
                "patch_hits": sum(int(row["patch_hits"]) for row in rows),
            }
        )
    return summaries


def validate_args(args: argparse.Namespace) -> None:
    if not args.flat_dir.is_dir():
        raise ValueError(f"Flat matrix directory not found: {args.flat_dir}")
    if not args.query_list.is_file():
        raise ValueError(f"Query list not found: {args.query_list}")
    if not math.isfinite(args.patch_hit_threshold):
        raise ValueError("--patch-hit-threshold must be finite")
    if not 0.0 <= args.patch_hit_threshold <= 1.0:
        raise ValueError("--patch-hit-threshold must be between 0 and 1")
    if args.block_size < 1:
        raise ValueError("--block-size must be at least 1")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")


def main() -> int:
    args = parse_args()
    args.flat_dir = args.flat_dir.resolve()
    args.output_root = args.output_root.resolve()
    args.query_list = args.query_list.resolve()
    configure_logging(args.output_root / "logs")

    try:
        validate_args(args)
        queries = read_query_list(args.query_list)
        dataset, flat_files = load_patch_dataset(args.flat_dir, args.target_match)
        missing = [query for query in queries if query not in dataset.indices_by_allergen]
        if missing:
            raise ValueError(
                "Query allergens absent from connectivity matrices: " + ", ".join(missing)
            )
        file_keys = [sanitize_stem(query) for query in queries]
        if len(file_keys) != len(set(file_keys)):
            raise ValueError("Query allergen names collide after filename sanitization")

        signature = build_input_signature(
            flat_files, args.patch_hit_threshold, args.exclude_self
        )
        logging.info(
            "Loaded %d patches from %d files across %d allergens",
            dataset.patch_count,
            len(flat_files),
            len(dataset.indices_by_allergen),
        )
        logging.info(
            "Selected queries=%d threshold=%.12g exclude_self=%s workers=%d block_size=%d",
            len(queries),
            args.patch_hit_threshold,
            args.exclude_self,
            args.workers,
            args.block_size,
        )

        completed: Dict[str, QueryPaths] = {}
        if args.workers == 1:
            for query in queries:
                completed[query] = score_one_query(
                    query,
                    dataset,
                    args.output_root,
                    args.patch_hit_threshold,
                    args.exclude_self,
                    args.block_size,
                    signature,
                    args.resume,
                )
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(
                        score_one_query,
                        query,
                        dataset,
                        args.output_root,
                        args.patch_hit_threshold,
                        args.exclude_self,
                        args.block_size,
                        signature,
                        args.resume,
                    ): query
                    for query in queries
                }
                for future in as_completed(futures):
                    query = futures[future]
                    completed[query] = future.result()

        ordered_paths = [completed[query] for query in queries]
        combine_query_csvs(
            args.output_root / "pair_scores" / "allergen_pair_best_scores.csv",
            PAIR_SCORE_COLUMNS,
            [paths.pair_file for paths in ordered_paths],
        )
        combine_query_csvs(
            args.output_root / "rankings" / "allergen_rankings.csv",
            RANKING_COLUMNS,
            [paths.ranking_file for paths in ordered_paths],
        )
        summaries = build_summary(ordered_paths)
        write_csv_atomic(
            args.output_root / "summary" / "allergen_scoring_summary.csv",
            SUMMARY_COLUMNS,
            summaries,
        )
        write_json_atomic(
            args.output_root / "run_configuration.json",
            {
                "algorithm_version": ALGORITHM_VERSION,
                "block_size": args.block_size,
                "exclude_self": args.exclude_self,
                "flat_dir": str(args.flat_dir),
                "input_signature": signature,
                "patch_hit_threshold": args.patch_hit_threshold,
                "queries": queries,
                "query_list": str(args.query_list),
                "target_match": args.target_match,
                "workers": args.workers,
            },
        )
        logging.info("All selected queries completed: %d", len(queries))
        return 0
    except Exception:
        logging.exception("Allergen patch-pair scoring failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
