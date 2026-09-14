#!/usr/bin/env python3
"""
Build a query connectivity matrix and score it against target patch matrices using
the legacy Cross-React correlation formula.

Default behavior:
- Target flat matrix directory: ./3_connectivity/connectivity_matrix_flat
- Query matrix output: ./4_pcc/query_matrix
- PCC score output: ./4_pcc/pcc_scores
- Ranked result output: ./4_pcc/pcc_ranked
- Logs: ./4_pcc/pcc_logs

Expected query epitope CSV columns:
- residue_number
- chain_id (optional; required only if the same residue number exists in multiple chains)
- residue_name (optional but recommended)

Representative atom rules:
- non-glycine: CB
- glycine: CA

Connectivity rule:
- representative-atom distance <= 8.0 A
- each connected residue pair increments both [i,j] and [j,i]
- same-type contacts therefore add 2 to the diagonal cell
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
from typing import Dict, List, Sequence, Tuple

from legacy_pcc import build_legacy_stats, legacy_crossreact_score_from_stats


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = Path.cwd().resolve()
STAGE_ROOT = RUN_ROOT / "4_pcc"
DEFAULT_TARGET_FLAT_DIR = RUN_ROOT / "3_connectivity" / "connectivity_matrix_flat"
DEFAULT_QUERY_MATRIX_DIR = STAGE_ROOT / "query_matrix"
DEFAULT_PCC_SCORES_DIR = STAGE_ROOT / "pcc_scores"
DEFAULT_PCC_RANKED_DIR = STAGE_ROOT / "pcc_ranked"
DEFAULT_LOG_DIR = STAGE_ROOT / "pcc_logs"

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
MATRIX_COLUMNS = [f"m_{aa_row}_{aa_col}" for aa_row in AA_ORDER for aa_col in AA_ORDER]
MATRIX_SIZE = len(AA_ORDER)
MATRIX_ELEMENT_COUNT = MATRIX_SIZE * MATRIX_SIZE


@dataclass(frozen=True)
class QueryEpitopeRow:
    chain_id: str | None
    residue_number: str
    residue_name: str | None


@dataclass(frozen=True)
class RepresentativeAtom:
    chain_id: str
    residue_number: str
    residue_name: str
    atom_name: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class QueryResidue:
    chain_id: str
    residue_number: str
    residue_name: str
    atom_name: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class TargetPatchRow:
    patch_id: str
    source_pdb: str
    accession: str
    allergen: str
    patch_size: int
    connected_pair_count: int
    matrix_total_sum: int
    unknown_pair_count: int
    legacy_row_sums: Tuple[float, ...]
    legacy_sum: float
    legacy_sumsq: float
    flat_source_file: str


@dataclass(frozen=True)
class PccScoreRow:
    query_name: str
    target_flat_file: str
    target_source_pdb: str
    target_accession: str
    target_allergen: str
    patch_id: str
    patch_size: int
    connected_pair_count: int
    matrix_total_sum: int
    unknown_pair_count: int
    pcc: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a query matrix and calculate PCC scores against target patches."
    )
    parser.add_argument(
        "--query-pdb",
        type=Path,
        required=True,
        help="Path to the query PDB file",
    )
    parser.add_argument(
        "--query-epitope-csv",
        type=Path,
        required=True,
        help="CSV file listing query epitope residues",
    )
    parser.add_argument(
        "--target-flat-dir",
        type=Path,
        default=DEFAULT_TARGET_FLAT_DIR,
        help=f"Target flat matrix directory (default: {DEFAULT_TARGET_FLAT_DIR})",
    )
    parser.add_argument(
        "--query-matrix-dir",
        type=Path,
        default=DEFAULT_QUERY_MATRIX_DIR,
        help=f"Query matrix output directory (default: {DEFAULT_QUERY_MATRIX_DIR})",
    )
    parser.add_argument(
        "--pcc-scores-dir",
        type=Path,
        default=DEFAULT_PCC_SCORES_DIR,
        help=f"PCC score output directory (default: {DEFAULT_PCC_SCORES_DIR})",
    )
    parser.add_argument(
        "--pcc-ranked-dir",
        type=Path,
        default=DEFAULT_PCC_RANKED_DIR,
        help=f"Ranked result output directory (default: {DEFAULT_PCC_RANKED_DIR})",
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
        "--target-match",
        default="*.connectivity_flat.csv",
        help="Glob pattern for target flat matrix files (default: *.connectivity_flat.csv)",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "pcc_run.log"
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


def parse_representative_atoms(pdb_path: Path) -> Dict[Tuple[str, str], RepresentativeAtom]:
    atoms: Dict[Tuple[str, str], RepresentativeAtom] = {}

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            if len(line) < 54:
                continue

            alt_loc = line[16].strip()
            if alt_loc not in {"", "A"}:
                continue

            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()
            if insertion_code:
                residue_number = f"{residue_number}{insertion_code}"

            target_atom = "CA" if residue_name == "GLY" else "CB"
            if atom_name != target_atom:
                continue

            key = (chain_id, residue_number)
            if key in atoms:
                continue

            atoms[key] = RepresentativeAtom(
                chain_id=chain_id,
                residue_number=residue_number,
                residue_name=residue_name,
                atom_name=atom_name,
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )

    return atoms


def load_query_epitope_csv(path: Path) -> List[QueryEpitopeRow]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"residue_number"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Query epitope CSV must contain columns {sorted(required)}; missing {sorted(missing)}"
            )

        rows = []
        for row in reader:
            residue_number = row.get("residue_number", "").strip()
            if not residue_number:
                continue
            residue_name = row.get("residue_name", "").strip() or None
            chain_id = row.get("chain_id", "").strip() or None
            rows.append(
                QueryEpitopeRow(
                    chain_id=chain_id,
                    residue_number=residue_number,
                    residue_name=residue_name,
                )
            )

    if not rows:
        raise ValueError(f"Query epitope CSV is empty: {path.name}")
    return rows


def build_query_residue_set(
    epitope_rows: Sequence[QueryEpitopeRow],
    representative_atoms: Dict[Tuple[str, str], RepresentativeAtom],
) -> List[QueryResidue]:
    residues: List[QueryResidue] = []
    seen = set()
    atoms_by_residue_number: Dict[str, List[RepresentativeAtom]] = {}

    for atom in representative_atoms.values():
        atoms_by_residue_number.setdefault(atom.residue_number, []).append(atom)

    for row in epitope_rows:
        if row.chain_id:
            key = (row.chain_id, row.residue_number)
            atom = representative_atoms.get(key)
            if atom is None:
                raise KeyError(
                    f"Representative atom not found for query residue {row.chain_id}:{row.residue_number}"
                )
        else:
            matches = atoms_by_residue_number.get(row.residue_number, [])
            if not matches:
                raise KeyError(
                    f"Representative atom not found for query residue {row.residue_number}"
                )
            if len(matches) > 1:
                chains = ", ".join(sorted(atom.chain_id or "(blank)" for atom in matches))
                raise KeyError(
                    "Multiple chains match query residue "
                    f"{row.residue_number} ({chains}). Add chain_id to the CSV."
                )
            atom = matches[0]

        if row.residue_name and row.residue_name != atom.residue_name:
            residue_label = (
                f"{atom.chain_id}:{row.residue_number}" if atom.chain_id else row.residue_number
            )
            raise ValueError(
                f"Residue name mismatch for {residue_label}: "
                f"CSV has {row.residue_name}, PDB has {atom.residue_name}"
            )
        unique_key = (atom.chain_id, atom.residue_number, atom.residue_name)
        if unique_key in seen:
            continue
        seen.add(unique_key)
        residues.append(
            QueryResidue(
                chain_id=atom.chain_id,
                residue_number=atom.residue_number,
                residue_name=atom.residue_name,
                atom_name=atom.atom_name,
                x=atom.x,
                y=atom.y,
                z=atom.z,
            )
        )

    return residues


def euclidean_distance(a: QueryResidue, b: QueryResidue) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def build_empty_matrix() -> List[List[int]]:
    return [[0 for _ in AA_ORDER] for _ in AA_ORDER]


def build_connectivity_matrix(
    residues: Sequence[QueryResidue], cutoff: float
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


def tuple_matrix_from_flat(flat_values: Sequence[int]) -> Tuple[Tuple[int, ...], ...]:
    if len(flat_values) != MATRIX_ELEMENT_COUNT:
        raise ValueError(
            f"Expected {MATRIX_ELEMENT_COUNT} matrix values, got {len(flat_values)}"
        )
    return tuple(
        tuple(int(flat_values[row * MATRIX_SIZE + col]) for col in range(MATRIX_SIZE))
        for row in range(MATRIX_SIZE)
    )


def iter_target_flat_files(target_flat_dir: Path, match: str) -> List[Path]:
    return sorted(target_flat_dir.glob(match), key=lambda path: path.name)


def load_target_flat_rows(target_flat_dir: Path, match: str) -> List[TargetPatchRow]:
    rows: List[TargetPatchRow] = []
    for path in iter_target_flat_files(target_flat_dir, match):
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            missing = [column for column in MATRIX_COLUMNS if column not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"Target flat file missing matrix columns: {path.name}")

            for row in reader:
                matrix = tuple_matrix_from_flat(
                    [int(float(row[column])) for column in MATRIX_COLUMNS]
                )
                legacy_row_sums, legacy_sum, legacy_sumsq = build_legacy_stats(matrix)
                rows.append(
                    TargetPatchRow(
                        patch_id=row["patch_id"],
                        source_pdb=row["source_pdb"],
                        accession=row["accession"],
                        allergen=row["allergen"],
                        patch_size=int(row["patch_size"]),
                        connected_pair_count=int(row["connected_pair_count"]),
                        matrix_total_sum=int(row["matrix_total_sum"]),
                        unknown_pair_count=int(row["unknown_pair_count"]),
                        legacy_row_sums=legacy_row_sums,
                        legacy_sum=legacy_sum,
                        legacy_sumsq=legacy_sumsq,
                        flat_source_file=path.name,
                    )
                )

    if not rows:
        raise ValueError(f"No target flat matrix rows found in {target_flat_dir}")
    return rows


def score_query_against_targets(
    query_name: str,
    query_matrix: Sequence[Sequence[int]],
    target_rows: Sequence[TargetPatchRow],
) -> List[PccScoreRow]:
    query_row_sums, query_sum, query_sumsq = build_legacy_stats(query_matrix)
    scores: List[PccScoreRow] = []
    for row in target_rows:
        scores.append(
            PccScoreRow(
                query_name=query_name,
                target_flat_file=row.flat_source_file,
                target_source_pdb=row.source_pdb,
                target_accession=row.accession,
                target_allergen=row.allergen,
                patch_id=row.patch_id,
                patch_size=row.patch_size,
                connected_pair_count=row.connected_pair_count,
                matrix_total_sum=row.matrix_total_sum,
                unknown_pair_count=row.unknown_pair_count,
                pcc=legacy_crossreact_score_from_stats(
                    query_row_sums=query_row_sums,
                    query_sum=query_sum,
                    query_sumsq=query_sumsq,
                    target_row_sums=row.legacy_row_sums,
                    target_sum=row.legacy_sum,
                    target_sumsq=row.legacy_sumsq,
                ),
            )
        )
    return scores


def pcc_sort_key(row: PccScoreRow) -> Tuple[int, float, str]:
    if math.isnan(row.pcc):
        return (1, 0.0, row.patch_id)
    return (0, -row.pcc, row.patch_id)


def select_best_patch_per_target(
    rows: Sequence[PccScoreRow],
) -> List[PccScoreRow]:
    best: Dict[Tuple[str, str, str], PccScoreRow] = {}
    for row in rows:
        key = (row.target_source_pdb, row.target_accession, row.target_allergen)
        current = best.get(key)
        if current is None or pcc_sort_key(row) < pcc_sort_key(current):
            best[key] = row
    return sorted(best.values(), key=pcc_sort_key)


def write_query_matrix_csv(
    output_path: Path,
    matrix: Sequence[Sequence[int]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["aa"] + AA_ORDER)
        for aa, row in zip(AA_ORDER, matrix):
            writer.writerow([aa] + [int(value) for value in row])


def write_query_flat_csv(
    output_path: Path,
    query_name: str,
    matrix: Sequence[Sequence[int]],
    residue_count: int,
    connected_pair_count: int,
    unknown_pair_count: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "query_name",
                "residue_count",
                "connected_pair_count",
                "matrix_total_sum",
                "unknown_pair_count",
            ]
            + MATRIX_COLUMNS
        )
        writer.writerow(
            [
                query_name,
                residue_count,
                connected_pair_count,
                matrix_total_sum(matrix),
                unknown_pair_count,
            ]
            + flatten_matrix(matrix)
        )


def write_query_residue_set_csv(
    output_path: Path,
    residues: Sequence[QueryResidue],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "chain_id",
                "residue_number",
                "residue_name",
                "atom_name",
                "x",
                "y",
                "z",
            ]
        )
        for residue in residues:
            writer.writerow(
                [
                    residue.chain_id,
                    residue.residue_number,
                    residue.residue_name,
                    residue.atom_name,
                    f"{residue.x:.6f}",
                    f"{residue.y:.6f}",
                    f"{residue.z:.6f}",
                ]
            )


def write_pcc_scores_csv(
    output_path: Path,
    rows: Sequence[PccScoreRow],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "query_name",
                "target_flat_file",
                "target_source_pdb",
                "target_accession",
                "target_allergen",
                "patch_id",
                "patch_size",
                "connected_pair_count",
                "matrix_total_sum",
                "unknown_pair_count",
                "pcc",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.query_name,
                    row.target_flat_file,
                    row.target_source_pdb,
                    row.target_accession,
                    row.target_allergen,
                    row.patch_id,
                    row.patch_size,
                    row.connected_pair_count,
                    row.matrix_total_sum,
                    row.unknown_pair_count,
                    "" if math.isnan(row.pcc) else f"{row.pcc:.12f}",
                ]
            )


def main() -> int:
    args = parse_args()
    configure_logging(args.log_dir)

    query_pdb = args.query_pdb.resolve()
    query_epitope_csv = args.query_epitope_csv.resolve()
    target_flat_dir = args.target_flat_dir.resolve()
    query_matrix_dir = args.query_matrix_dir.resolve()
    pcc_scores_dir = args.pcc_scores_dir.resolve()
    pcc_ranked_dir = args.pcc_ranked_dir.resolve()

    if not query_pdb.exists():
        logging.error("Query PDB not found: %s", query_pdb)
        return 1
    if not query_epitope_csv.exists():
        logging.error("Query epitope CSV not found: %s", query_epitope_csv)
        return 1
    if not target_flat_dir.exists():
        logging.error("Target flat matrix directory not found: %s", target_flat_dir)
        return 1

    query_name = sanitize_stem(query_epitope_csv.stem)
    logging.info("Script directory: %s", SCRIPT_DIR)
    logging.info("Run root: %s", RUN_ROOT)
    logging.info("Stage root: %s", STAGE_ROOT)
    logging.info("Query PDB: %s", query_pdb)
    logging.info("Query epitope CSV: %s", query_epitope_csv)
    logging.info("Target flat matrix directory: %s", target_flat_dir)
    logging.info("Connectivity cutoff: %.3f A", args.cutoff)

    epitope_rows = load_query_epitope_csv(query_epitope_csv)
    representative_atoms = parse_representative_atoms(query_pdb)
    query_residues = build_query_residue_set(epitope_rows, representative_atoms)
    query_matrix, query_connected_pairs, query_unknown_pairs = build_connectivity_matrix(
        query_residues, args.cutoff
    )
    write_query_residue_set_csv(
        query_matrix_dir / f"{query_name}.query_residue_set.csv",
        query_residues,
    )
    write_query_matrix_csv(
        query_matrix_dir / f"{query_name}.query_matrix.csv",
        query_matrix,
    )
    write_query_flat_csv(
        query_matrix_dir / f"{query_name}.query_flat.csv",
        query_name,
        query_matrix,
        len(query_residues),
        query_connected_pairs,
        query_unknown_pairs,
    )

    target_rows = load_target_flat_rows(target_flat_dir, args.target_match)
    scores = score_query_against_targets(query_name, query_matrix, target_rows)
    scores_sorted = sorted(scores, key=pcc_sort_key)
    best_per_target = select_best_patch_per_target(scores_sorted)

    write_pcc_scores_csv(
        pcc_scores_dir / f"{query_name}.pcc_scores.csv",
        scores,
    )
    write_pcc_scores_csv(
        pcc_ranked_dir / f"{query_name}.pcc_scores.ranked.csv",
        scores_sorted,
    )
    write_pcc_scores_csv(
        pcc_ranked_dir / f"{query_name}.best_patch_per_target.csv",
        best_per_target,
    )

    nan_count = sum(math.isnan(row.pcc) for row in scores)
    logging.info(
        "Done. query_residues=%d target_patches=%d nan_pcc=%d",
        len(query_residues),
        len(scores),
        nan_count,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
