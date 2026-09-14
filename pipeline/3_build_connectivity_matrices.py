#!/usr/bin/env python3
"""
Build 20x20 connectivity matrices from 10 A patch membership outputs.

Default behavior:
- Input patch member directory: ./2_patches_10A/patch_members
- Patch residue set output: ./3_connectivity/patch_residue_sets
- Square matrix output: ./3_connectivity/connectivity_matrices
- Flattened matrix output: ./3_connectivity/connectivity_matrix_flat
- Summary output: ./3_connectivity/connectivity_summary
- Logs: ./3_connectivity/connectivity_logs

Connectivity definition:
- Candidate residues inside each patch are taken from patch_members outputs
- Residues are deduplicated within each patch by (chain_id, residue_number, residue_name)
- Representative atom coordinates come from member atom coordinates
- Two residues are connected when representative-atom distance <= 8.0 A
- Each connected residue pair increments both [i,j] and [j,i]
- Same-type contacts therefore add 2 to the diagonal cell
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = Path.cwd().resolve()
STAGE_ROOT = RUN_ROOT / "3_connectivity"
DEFAULT_PATCH_MEMBERS_DIR = RUN_ROOT / "2_patches_10A" / "patch_members"
DEFAULT_RESIDUE_SET_DIR = STAGE_ROOT / "patch_residue_sets"
DEFAULT_MATRIX_DIR = STAGE_ROOT / "connectivity_matrices"
DEFAULT_FLAT_DIR = STAGE_ROOT / "connectivity_matrix_flat"
DEFAULT_SUMMARY_DIR = STAGE_ROOT / "connectivity_summary"
DEFAULT_LOG_DIR = STAGE_ROOT / "connectivity_logs"

DEFAULT_CUTOFF = 8.0
AA_ORDER = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]
AA_INDEX = {aa: idx for idx, aa in enumerate(AA_ORDER)}


@dataclass(frozen=True)
class PatchMemberRow:
    patch_id: str
    source_pdb: str
    accession: str
    allergen: str
    center_chain_id: str
    center_residue_number: str
    center_residue_name: str
    center_atom_name: str
    center_x: float
    center_y: float
    center_z: float
    member_chain_id: str
    member_residue_number: str
    member_residue_name: str
    member_atom_name: str
    member_x: float
    member_y: float
    member_z: float
    distance_to_center: float
    is_center: bool


@dataclass(frozen=True)
class PatchResidue:
    patch_id: str
    source_pdb: str
    accession: str
    allergen: str
    chain_id: str
    residue_number: str
    residue_name: str
    atom_name: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class PatchMatrixSummary:
    patch_id: str
    source_pdb: str
    accession: str
    allergen: str
    patch_size: int
    connected_pair_count: int
    matrix_total_sum: int
    unknown_pair_count: int


@dataclass(frozen=True)
class FileSummary:
    patch_members_file: str
    source_pdb: str
    accession: str
    allergen: str
    patch_count: int
    total_patch_residue_count: int
    total_connected_pairs: int
    total_unknown_pairs: int
    status: str
    error: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build connectivity matrices from patch membership CSV files."
        )
    )
    parser.add_argument(
        "--patch-members-dir",
        type=Path,
        default=DEFAULT_PATCH_MEMBERS_DIR,
        help=f"Input patch member directory (default: {DEFAULT_PATCH_MEMBERS_DIR})",
    )
    parser.add_argument(
        "--residue-set-dir",
        type=Path,
        default=DEFAULT_RESIDUE_SET_DIR,
        help=f"Patch residue set output directory (default: {DEFAULT_RESIDUE_SET_DIR})",
    )
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        default=DEFAULT_MATRIX_DIR,
        help=f"Square matrix output directory (default: {DEFAULT_MATRIX_DIR})",
    )
    parser.add_argument(
        "--flat-dir",
        type=Path,
        default=DEFAULT_FLAT_DIR,
        help=f"Flattened matrix output directory (default: {DEFAULT_FLAT_DIR})",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=DEFAULT_SUMMARY_DIR,
        help=f"Summary output directory (default: {DEFAULT_SUMMARY_DIR})",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Log output directory (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF,
        help=f"Connectivity cutoff in angstrom (default: {DEFAULT_CUTOFF})",
    )
    parser.add_argument(
        "--match",
        default="*.patch_members.csv",
        help="Glob pattern used to select patch member files (default: *.patch_members.csv)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N matching files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing per-file outputs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which files would be processed without building matrices",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "connectivity_run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def sanitize_stem(name: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return sanitized or "unnamed"


def base_stem_from_patch_members(path: Path) -> str:
    suffix = ".patch_members.csv"
    name = path.name
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return path.stem


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def iter_patch_member_files(
    patch_members_dir: Path, match: str, limit: int | None
) -> List[Path]:
    files = sorted(patch_members_dir.glob(match), key=lambda path: path.name)
    if limit is not None:
        files = files[:limit]
    return files


def load_patch_member_rows(path: Path) -> List[PatchMemberRow]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [
            PatchMemberRow(
                patch_id=row["patch_id"],
                source_pdb=row["source_pdb"],
                accession=row["accession"],
                allergen=row["allergen"],
                center_chain_id=row["center_chain_id"],
                center_residue_number=row["center_residue_number"],
                center_residue_name=row["center_residue_name"],
                center_atom_name=row["center_atom_name"],
                center_x=float(row["center_x"]),
                center_y=float(row["center_y"]),
                center_z=float(row["center_z"]),
                member_chain_id=row["member_chain_id"],
                member_residue_number=row["member_residue_number"],
                member_residue_name=row["member_residue_name"],
                member_atom_name=row["member_atom_name"],
                member_x=float(row["member_x"]),
                member_y=float(row["member_y"]),
                member_z=float(row["member_z"]),
                distance_to_center=float(row["distance_to_center"]),
                is_center=parse_bool(row["is_center"]),
            )
            for row in reader
        ]
    if not rows:
        raise ValueError(f"Patch member CSV is empty: {path.name}")
    return rows


def group_rows_by_patch_id(
    rows: Sequence[PatchMemberRow],
) -> Dict[str, List[PatchMemberRow]]:
    grouped: Dict[str, List[PatchMemberRow]] = {}
    for row in rows:
        grouped.setdefault(row.patch_id, []).append(row)
    return grouped


def deduplicate_patch_members(
    grouped_rows: Dict[str, List[PatchMemberRow]]
) -> Dict[str, List[PatchResidue]]:
    deduped: Dict[str, List[PatchResidue]] = {}

    for patch_id, rows in grouped_rows.items():
        residues: List[PatchResidue] = []
        seen = set()
        for row in rows:
            key = (
                row.member_chain_id,
                row.member_residue_number,
                row.member_residue_name,
            )
            if key in seen:
                continue
            seen.add(key)
            residues.append(
                PatchResidue(
                    patch_id=patch_id,
                    source_pdb=row.source_pdb,
                    accession=row.accession,
                    allergen=row.allergen,
                    chain_id=row.member_chain_id,
                    residue_number=row.member_residue_number,
                    residue_name=row.member_residue_name,
                    atom_name=row.member_atom_name,
                    x=row.member_x,
                    y=row.member_y,
                    z=row.member_z,
                )
            )
        deduped[patch_id] = residues

    return deduped


def euclidean_distance(a: PatchResidue, b: PatchResidue) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def build_empty_matrix() -> List[List[int]]:
    return [[0 for _ in AA_ORDER] for _ in AA_ORDER]


def build_connectivity_matrix(
    residues: Sequence[PatchResidue], cutoff: float
) -> Tuple[List[List[int]], int, int]:
    matrix = build_empty_matrix()
    connected_pairs = 0
    unknown_pairs = 0

    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            left = residues[i]
            right = residues[j]
            distance = euclidean_distance(left, right)
            if distance > cutoff:
                continue

            idx_left = AA_INDEX.get(left.residue_name)
            idx_right = AA_INDEX.get(right.residue_name)
            if idx_left is None or idx_right is None:
                unknown_pairs += 1
                continue

            matrix[idx_left][idx_right] += 1
            matrix[idx_right][idx_left] += 1
            connected_pairs += 1

    return matrix, connected_pairs, unknown_pairs


def flatten_matrix(matrix: Sequence[Sequence[int]]) -> List[int]:
    flat: List[int] = []
    for row in matrix:
        flat.extend(int(value) for value in row)
    return flat


def matrix_total_sum(matrix: Sequence[Sequence[int]]) -> int:
    return sum(sum(int(value) for value in row) for row in matrix)


def write_patch_residue_sets_csv(
    output_path: Path, patch_residue_sets: Dict[str, List[PatchResidue]]
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_id",
                "source_pdb",
                "accession",
                "allergen",
                "chain_id",
                "residue_number",
                "residue_name",
                "atom_name",
                "x",
                "y",
                "z",
            ]
        )
        for patch_id in sorted(patch_residue_sets):
            for residue in patch_residue_sets[patch_id]:
                writer.writerow(
                    [
                        residue.patch_id,
                        residue.source_pdb,
                        residue.accession,
                        residue.allergen,
                        residue.chain_id,
                        residue.residue_number,
                        residue.residue_name,
                        residue.atom_name,
                        f"{residue.x:.6f}",
                        f"{residue.y:.6f}",
                        f"{residue.z:.6f}",
                    ]
                )


def write_square_matrix_csv(
    output_path: Path,
    matrix: Sequence[Sequence[int]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["aa"] + AA_ORDER)
        for aa, row in zip(AA_ORDER, matrix):
            writer.writerow([aa] + [int(value) for value in row])


def write_flat_matrix_csv(
    output_path: Path,
    summaries: Sequence[PatchMatrixSummary],
    matrices: Dict[str, List[List[int]]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    matrix_columns = [f"m_{aa_row}_{aa_col}" for aa_row in AA_ORDER for aa_col in AA_ORDER]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_id",
                "source_pdb",
                "accession",
                "allergen",
                "patch_size",
                "connected_pair_count",
                "matrix_total_sum",
                "unknown_pair_count",
            ]
            + matrix_columns
        )
        for summary in summaries:
            flat = flatten_matrix(matrices[summary.patch_id])
            writer.writerow(
                [
                    summary.patch_id,
                    summary.source_pdb,
                    summary.accession,
                    summary.allergen,
                    summary.patch_size,
                    summary.connected_pair_count,
                    summary.matrix_total_sum,
                    summary.unknown_pair_count,
                ]
                + flat
            )


def write_patch_summary_csv(
    output_path: Path,
    summaries: Sequence[PatchMatrixSummary],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_id",
                "source_pdb",
                "accession",
                "allergen",
                "patch_size",
                "connected_pair_count",
                "matrix_total_sum",
                "unknown_pair_count",
            ]
        )
        for summary in summaries:
            writer.writerow(
                [
                    summary.patch_id,
                    summary.source_pdb,
                    summary.accession,
                    summary.allergen,
                    summary.patch_size,
                    summary.connected_pair_count,
                    summary.matrix_total_sum,
                    summary.unknown_pair_count,
                ]
            )


def write_file_summary_csv(
    output_path: Path,
    rows: Sequence[FileSummary],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_members_file",
                "source_pdb",
                "accession",
                "allergen",
                "patch_count",
                "total_patch_residue_count",
                "total_connected_pairs",
                "total_unknown_pairs",
                "status",
                "error",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.patch_members_file,
                    row.source_pdb,
                    row.accession,
                    row.allergen,
                    row.patch_count,
                    row.total_patch_residue_count,
                    row.total_connected_pairs,
                    row.total_unknown_pairs,
                    row.status,
                    row.error,
                ]
            )


def count_csv_rows(csv_path: Path) -> int:
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def main() -> int:
    args = parse_args()
    configure_logging(args.log_dir)

    patch_members_dir = args.patch_members_dir.resolve()
    residue_set_dir = args.residue_set_dir.resolve()
    matrix_dir = args.matrix_dir.resolve()
    flat_dir = args.flat_dir.resolve()
    summary_dir = args.summary_dir.resolve()

    if not patch_members_dir.exists():
        logging.error("Patch member directory not found: %s", patch_members_dir)
        return 1

    patch_member_files = iter_patch_member_files(
        patch_members_dir, args.match, args.limit
    )
    if not patch_member_files:
        logging.error(
            "No patch member files matched pattern '%s' in %s",
            args.match,
            patch_members_dir,
        )
        return 1

    logging.info("Script directory: %s", SCRIPT_DIR)
    logging.info("Run root: %s", RUN_ROOT)
    logging.info("Stage root: %s", STAGE_ROOT)
    logging.info("Patch member directory: %s", patch_members_dir)
    logging.info("Matched patch member files: %d", len(patch_member_files))
    logging.info("Connectivity cutoff: %.3f A", args.cutoff)

    if args.dry_run:
        for patch_member_path in patch_member_files:
            logging.info("[DRY RUN] %s", patch_member_path.name)
        return 0

    file_summaries: List[FileSummary] = []

    for index, patch_member_path in enumerate(patch_member_files, start=1):
        base_stem = sanitize_stem(base_stem_from_patch_members(patch_member_path))
        residue_set_output = residue_set_dir / f"{base_stem}.patch_residue_sets.csv"
        flat_output = flat_dir / f"{base_stem}.connectivity_flat.csv"
        patch_summary_output = summary_dir / f"{base_stem}.connectivity_summary.csv"
        matrix_output_dir = matrix_dir / base_stem

        if (
            residue_set_output.exists()
            and flat_output.exists()
            and patch_summary_output.exists()
            and matrix_output_dir.exists()
            and not args.overwrite
        ):
            patch_count = count_csv_rows(patch_summary_output)
            logging.info(
                "[%d/%d] Skipping existing outputs: %s",
                index,
                len(patch_member_files),
                base_stem,
            )
            file_summaries.append(
                FileSummary(
                    patch_members_file=patch_member_path.name,
                    source_pdb="",
                    accession="",
                    allergen="",
                    patch_count=patch_count,
                    total_patch_residue_count=0,
                    total_connected_pairs=0,
                    total_unknown_pairs=0,
                    status="skipped_existing",
                    error="",
                )
            )
            continue

        try:
            rows = load_patch_member_rows(patch_member_path)
            source_pdb = rows[0].source_pdb
            accession = rows[0].accession
            allergen = rows[0].allergen

            grouped = group_rows_by_patch_id(rows)
            patch_residue_sets = deduplicate_patch_members(grouped)
            write_patch_residue_sets_csv(residue_set_output, patch_residue_sets)

            patch_summaries: List[PatchMatrixSummary] = []
            matrices: Dict[str, List[List[int]]] = {}

            matrix_output_dir.mkdir(parents=True, exist_ok=True)
            for patch_id in sorted(patch_residue_sets):
                residues = patch_residue_sets[patch_id]
                matrix, connected_pairs, unknown_pairs = build_connectivity_matrix(
                    residues, args.cutoff
                )
                matrices[patch_id] = matrix

                matrix_path = matrix_output_dir / f"{sanitize_stem(patch_id)}.connectivity_matrix.csv"
                write_square_matrix_csv(matrix_path, matrix)

                patch_summaries.append(
                    PatchMatrixSummary(
                        patch_id=patch_id,
                        source_pdb=residues[0].source_pdb if residues else source_pdb,
                        accession=residues[0].accession if residues else accession,
                        allergen=residues[0].allergen if residues else allergen,
                        patch_size=len(residues),
                        connected_pair_count=connected_pairs,
                        matrix_total_sum=matrix_total_sum(matrix),
                        unknown_pair_count=unknown_pairs,
                    )
                )

            write_flat_matrix_csv(flat_output, patch_summaries, matrices)
            write_patch_summary_csv(patch_summary_output, patch_summaries)

            total_patch_residue_count = sum(
                summary.patch_size for summary in patch_summaries
            )
            total_connected_pairs = sum(
                summary.connected_pair_count for summary in patch_summaries
            )
            total_unknown_pairs = sum(
                summary.unknown_pair_count for summary in patch_summaries
            )

            file_summaries.append(
                FileSummary(
                    patch_members_file=patch_member_path.name,
                    source_pdb=source_pdb,
                    accession=accession,
                    allergen=allergen,
                    patch_count=len(patch_summaries),
                    total_patch_residue_count=total_patch_residue_count,
                    total_connected_pairs=total_connected_pairs,
                    total_unknown_pairs=total_unknown_pairs,
                    status="ok",
                    error="",
                )
            )
            logging.info(
                "[%d/%d] %s -> patches=%d connected_pairs=%d",
                index,
                len(patch_member_files),
                patch_member_path.name,
                len(patch_summaries),
                total_connected_pairs,
            )
        except Exception as exc:  # pragma: no cover - file-level failure path
            logging.exception("Failed to build connectivity matrices for %s", patch_member_path.name)
            file_summaries.append(
                FileSummary(
                    patch_members_file=patch_member_path.name,
                    source_pdb="",
                    accession="",
                    allergen="",
                    patch_count=0,
                    total_patch_residue_count=0,
                    total_connected_pairs=0,
                    total_unknown_pairs=0,
                    status="error",
                    error=str(exc),
                )
            )

    write_file_summary_csv(
        summary_dir / "connectivity_generation_counts.csv",
        file_summaries,
    )

    ok_count = sum(row.status == "ok" for row in file_summaries)
    error_count = sum(row.status == "error" for row in file_summaries)
    skipped_count = sum(row.status == "skipped_existing" for row in file_summaries)
    logging.info(
        "Done. ok=%d error=%d skipped=%d total=%d",
        ok_count,
        error_count,
        skipped_count,
        len(file_summaries),
    )
    return 0 if error_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
