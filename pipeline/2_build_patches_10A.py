#!/usr/bin/env python3
"""
Build 10 A surface patches from surface residue CSV outputs.

Default behavior:
- Input PDB directory: ./pdb_all
- Input surface residue CSV directory: ./1_surface_residues/surface_residue_csv
- Surface center coordinate output: ./2_patches_10A/surface_center_coords
- Patch member output: ./2_patches_10A/patch_members
- Patch summary output: ./2_patches_10A/patch_summary
- Logs: ./2_patches_10A/patch_logs

Patch definition:
- Candidate residues: only rows with is_surface_gt10A2 == True
- Representative atom: CB for non-glycine, CA for glycine
- Patch membership: representative-atom distance <= 10.0 A
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


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = Path.cwd().resolve()
STAGE_ROOT = RUN_ROOT / "2_patches_10A"
DEFAULT_PDB_DIR = RUN_ROOT / "pdb_all"
DEFAULT_SURFACE_CSV_DIR = RUN_ROOT / "1_surface_residues" / "surface_residue_csv"
DEFAULT_CENTER_DIR = STAGE_ROOT / "surface_center_coords"
DEFAULT_MEMBERS_DIR = STAGE_ROOT / "patch_members"
DEFAULT_SUMMARY_DIR = STAGE_ROOT / "patch_summary"
DEFAULT_LOG_DIR = STAGE_ROOT / "patch_logs"

DEFAULT_RADIUS = 10.0


@dataclass(frozen=True)
class SurfaceResidue:
    order_index: int
    source_pdb: str
    accession: str
    allergen: str
    chain_id: str
    residue_number: str
    residue_name: str


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
class CenterRow:
    order_index: int
    source_pdb: str
    accession: str
    allergen: str
    chain_id: str
    residue_number: str
    residue_name: str
    center_atom_name: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class PatchMember:
    center_order_index: int
    member_order_index: int
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
class PatchSummary:
    center_order_index: int
    patch_id: str
    source_pdb: str
    accession: str
    allergen: str
    center_chain_id: str
    center_residue_number: str
    center_residue_name: str
    center_atom_name: str
    patch_size: int


@dataclass(frozen=True)
class FileSummary:
    surface_csv: str
    source_pdb: str
    accession: str
    allergen: str
    surface_residue_count: int
    center_count: int
    missing_center_count: int
    patch_count: int
    status: str
    error: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build 10 A patches using only surface residues from FreeSASA CSV outputs."
        )
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=DEFAULT_PDB_DIR,
        help=f"Input PDB directory (default: {DEFAULT_PDB_DIR})",
    )
    parser.add_argument(
        "--surface-csv-dir",
        type=Path,
        default=DEFAULT_SURFACE_CSV_DIR,
        help=f"Input surface residue CSV directory (default: {DEFAULT_SURFACE_CSV_DIR})",
    )
    parser.add_argument(
        "--center-dir",
        type=Path,
        default=DEFAULT_CENTER_DIR,
        help=f"Surface center coordinate output directory (default: {DEFAULT_CENTER_DIR})",
    )
    parser.add_argument(
        "--members-dir",
        type=Path,
        default=DEFAULT_MEMBERS_DIR,
        help=f"Patch member output directory (default: {DEFAULT_MEMBERS_DIR})",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=DEFAULT_SUMMARY_DIR,
        help=f"Patch summary output directory (default: {DEFAULT_SUMMARY_DIR})",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Log output directory (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=DEFAULT_RADIUS,
        help=f"Patch radius in angstrom (default: {DEFAULT_RADIUS})",
    )
    parser.add_argument(
        "--match",
        default="*.csv",
        help="Glob pattern used to select surface residue CSV files (default: *.csv)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N matching CSV files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing per-file outputs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which files would be processed without building patches",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "patch_10A_run.log"
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


def sanitize_token(token: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", token).strip("_")
    return sanitized or "blank"


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "y"}


def iter_surface_csv_files(
    surface_csv_dir: Path, match: str, limit: int | None
) -> List[Path]:
    files = sorted(surface_csv_dir.glob(match), key=lambda path: path.name)
    if limit is not None:
        files = files[:limit]
    return files


def load_surface_rows(
    surface_csv_path: Path,
) -> Tuple[str, str, str, List[SurfaceResidue]]:
    with surface_csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        all_rows = list(reader)

    if not all_rows:
        raise ValueError(f"Surface residue CSV is empty: {surface_csv_path.name}")

    first = all_rows[0]
    source_pdb = first["source_pdb"]
    accession = first["accession"]
    allergen = first["allergen"]

    surface_rows: List[SurfaceResidue] = []
    for row_index, row in enumerate(all_rows):
        if not parse_bool(row["is_surface_gt10A2"]):
            continue
        surface_rows.append(
            SurfaceResidue(
                order_index=row_index,
                source_pdb=row["source_pdb"],
                accession=row["accession"],
                allergen=row["allergen"],
                chain_id=row["chain_id"],
                residue_number=row["residue_number"],
                residue_name=row["residue_name"],
            )
        )

    return source_pdb, accession, allergen, surface_rows


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


def build_center_rows(
    surface_rows: Sequence[SurfaceResidue],
    representative_atoms: Dict[Tuple[str, str], RepresentativeAtom],
) -> Tuple[List[CenterRow], List[SurfaceResidue]]:
    center_rows: List[CenterRow] = []
    missing_rows: List[SurfaceResidue] = []

    for row in surface_rows:
        key = (row.chain_id, row.residue_number)
        atom = representative_atoms.get(key)
        if atom is None:
            missing_rows.append(row)
            continue

        center_rows.append(
            CenterRow(
                order_index=row.order_index,
                source_pdb=row.source_pdb,
                accession=row.accession,
                allergen=row.allergen,
                chain_id=row.chain_id,
                residue_number=row.residue_number,
                residue_name=row.residue_name,
                center_atom_name=atom.atom_name,
                x=atom.x,
                y=atom.y,
                z=atom.z,
            )
        )

    return center_rows, missing_rows


def make_patch_id(center: CenterRow) -> str:
    return (
        f"{sanitize_token(center.chain_id)}_"
        f"{sanitize_token(center.residue_number)}_"
        f"{sanitize_token(center.residue_name)}"
    )


def euclidean_distance(a: CenterRow, b: CenterRow) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def build_patches(
    center_rows: Sequence[CenterRow], radius: float
) -> Tuple[List[PatchMember], List[PatchSummary]]:
    members: List[PatchMember] = []
    summaries: List[PatchSummary] = []

    for center in center_rows:
        patch_id = make_patch_id(center)
        patch_size = 0

        for member in center_rows:
            distance = euclidean_distance(center, member)
            if distance > radius:
                continue

            members.append(
                PatchMember(
                    center_order_index=center.order_index,
                    member_order_index=member.order_index,
                    patch_id=patch_id,
                    source_pdb=center.source_pdb,
                    accession=center.accession,
                    allergen=center.allergen,
                    center_chain_id=center.chain_id,
                    center_residue_number=center.residue_number,
                    center_residue_name=center.residue_name,
                    center_atom_name=center.center_atom_name,
                    center_x=center.x,
                    center_y=center.y,
                    center_z=center.z,
                    member_chain_id=member.chain_id,
                    member_residue_number=member.residue_number,
                    member_residue_name=member.residue_name,
                    member_atom_name=member.center_atom_name,
                    member_x=member.x,
                    member_y=member.y,
                    member_z=member.z,
                    distance_to_center=distance,
                    is_center=(
                        center.chain_id == member.chain_id
                        and center.residue_number == member.residue_number
                    ),
                )
            )
            patch_size += 1

        summaries.append(
            PatchSummary(
                center_order_index=center.order_index,
                patch_id=patch_id,
                source_pdb=center.source_pdb,
                accession=center.accession,
                allergen=center.allergen,
                center_chain_id=center.chain_id,
                center_residue_number=center.residue_number,
                center_residue_name=center.residue_name,
                center_atom_name=center.center_atom_name,
                patch_size=patch_size,
            )
        )

    members.sort(
        key=lambda row: (
            row.center_order_index,
            row.distance_to_center,
            row.member_order_index,
            row.member_chain_id,
            row.member_residue_number,
        )
    )
    summaries.sort(
        key=lambda row: (
            row.center_order_index,
            row.center_chain_id,
            row.center_residue_number,
        )
    )
    return members, summaries


def write_center_csv(output_path: Path, rows: Sequence[CenterRow]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "source_pdb",
                "accession",
                "allergen",
                "chain_id",
                "residue_number",
                "residue_name",
                "center_atom_name",
                "x",
                "y",
                "z",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.source_pdb,
                    row.accession,
                    row.allergen,
                    row.chain_id,
                    row.residue_number,
                    row.residue_name,
                    row.center_atom_name,
                    f"{row.x:.6f}",
                    f"{row.y:.6f}",
                    f"{row.z:.6f}",
                ]
            )


def write_patch_members_csv(output_path: Path, rows: Sequence[PatchMember]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_id",
                "source_pdb",
                "accession",
                "allergen",
                "center_chain_id",
                "center_residue_number",
                "center_residue_name",
                "center_atom_name",
                "center_x",
                "center_y",
                "center_z",
                "member_chain_id",
                "member_residue_number",
                "member_residue_name",
                "member_atom_name",
                "member_x",
                "member_y",
                "member_z",
                "distance_to_center",
                "is_center",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.patch_id,
                    row.source_pdb,
                    row.accession,
                    row.allergen,
                    row.center_chain_id,
                    row.center_residue_number,
                    row.center_residue_name,
                    row.center_atom_name,
                    f"{row.center_x:.6f}",
                    f"{row.center_y:.6f}",
                    f"{row.center_z:.6f}",
                    row.member_chain_id,
                    row.member_residue_number,
                    row.member_residue_name,
                    row.member_atom_name,
                    f"{row.member_x:.6f}",
                    f"{row.member_y:.6f}",
                    f"{row.member_z:.6f}",
                    f"{row.distance_to_center:.6f}",
                    str(row.is_center),
                ]
            )


def write_patch_summary_csv(output_path: Path, rows: Sequence[PatchSummary]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patch_id",
                "source_pdb",
                "accession",
                "allergen",
                "center_chain_id",
                "center_residue_number",
                "center_residue_name",
                "center_atom_name",
                "patch_size",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.patch_id,
                    row.source_pdb,
                    row.accession,
                    row.allergen,
                    row.center_chain_id,
                    row.center_residue_number,
                    row.center_residue_name,
                    row.center_atom_name,
                    row.patch_size,
                ]
            )


def write_file_summary_csv(output_path: Path, rows: Sequence[FileSummary]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "surface_csv",
                "source_pdb",
                "accession",
                "allergen",
                "surface_residue_count",
                "center_count",
                "missing_center_count",
                "patch_count",
                "status",
                "error",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.surface_csv,
                    row.source_pdb,
                    row.accession,
                    row.allergen,
                    row.surface_residue_count,
                    row.center_count,
                    row.missing_center_count,
                    row.patch_count,
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

    pdb_dir = args.pdb_dir.resolve()
    surface_csv_dir = args.surface_csv_dir.resolve()
    center_dir = args.center_dir.resolve()
    members_dir = args.members_dir.resolve()
    summary_dir = args.summary_dir.resolve()

    if not pdb_dir.exists():
        logging.error("PDB directory not found: %s", pdb_dir)
        return 1

    if not surface_csv_dir.exists():
        logging.error("Surface CSV directory not found: %s", surface_csv_dir)
        return 1

    surface_csv_files = iter_surface_csv_files(surface_csv_dir, args.match, args.limit)
    if not surface_csv_files:
        logging.error("No surface CSV files matched pattern '%s' in %s", args.match, surface_csv_dir)
        return 1

    logging.info("Script directory: %s", SCRIPT_DIR)
    logging.info("Run root: %s", RUN_ROOT)
    logging.info("Stage root: %s", STAGE_ROOT)
    logging.info("PDB directory: %s", pdb_dir)
    logging.info("Surface CSV directory: %s", surface_csv_dir)
    logging.info("Matched surface CSV files: %d", len(surface_csv_files))
    logging.info("Patch radius: %.3f A", args.radius)

    if args.dry_run:
        for surface_csv_path in surface_csv_files:
            logging.info("[DRY RUN] %s", surface_csv_path.name)
        return 0

    file_summaries: List[FileSummary] = []

    for index, surface_csv_path in enumerate(surface_csv_files, start=1):
        stem = sanitize_stem(surface_csv_path.stem)
        center_output_path = center_dir / f"{stem}.surface_centers.csv"
        members_output_path = members_dir / f"{stem}.patch_members.csv"
        summary_output_path = summary_dir / f"{stem}.patch_summary.csv"

        try:
            source_pdb, accession, allergen, surface_rows = load_surface_rows(surface_csv_path)
            surface_count = len(surface_rows)

            if (
                center_output_path.exists()
                and members_output_path.exists()
                and summary_output_path.exists()
                and not args.overwrite
            ):
                center_count = count_csv_rows(center_output_path)
                patch_count = count_csv_rows(summary_output_path)
                missing_count = max(surface_count - center_count, 0)
                logging.info(
                    "[%d/%d] Skipping existing outputs: %s",
                    index,
                    len(surface_csv_files),
                    stem,
                )
                file_summaries.append(
                    FileSummary(
                        surface_csv=surface_csv_path.name,
                        source_pdb=source_pdb,
                        accession=accession,
                        allergen=allergen,
                        surface_residue_count=surface_count,
                        center_count=center_count,
                        missing_center_count=missing_count,
                        patch_count=patch_count,
                        status="skipped_existing",
                        error="",
                    )
                )
                continue

            pdb_path = pdb_dir / source_pdb
            if not pdb_path.exists():
                raise FileNotFoundError(f"PDB file not found: {pdb_path}")

            representative_atoms = parse_representative_atoms(pdb_path)
            center_rows, missing_rows = build_center_rows(surface_rows, representative_atoms)
            patch_members, patch_summaries = build_patches(center_rows, args.radius)

            write_center_csv(center_output_path, center_rows)
            write_patch_members_csv(members_output_path, patch_members)
            write_patch_summary_csv(summary_output_path, patch_summaries)

            if missing_rows:
                preview = ", ".join(
                    f"{row.chain_id}:{row.residue_number}:{row.residue_name}"
                    for row in missing_rows[:5]
                )
                logging.warning(
                    "[%d/%d] Missing representative atom for %d surface residues in %s: %s",
                    index,
                    len(surface_csv_files),
                    len(missing_rows),
                    source_pdb,
                    preview,
                )

            file_summaries.append(
                FileSummary(
                    surface_csv=surface_csv_path.name,
                    source_pdb=source_pdb,
                    accession=accession,
                    allergen=allergen,
                    surface_residue_count=surface_count,
                    center_count=len(center_rows),
                    missing_center_count=len(missing_rows),
                    patch_count=len(patch_summaries),
                    status="ok",
                    error="",
                )
            )
            logging.info(
                "[%d/%d] %s -> surface=%d centers=%d patches=%d",
                index,
                len(surface_csv_files),
                surface_csv_path.name,
                surface_count,
                len(center_rows),
                len(patch_summaries),
            )
        except Exception as exc:  # pragma: no cover - file-level failure path
            logging.exception("Failed to build patches for %s", surface_csv_path.name)
            file_summaries.append(
                FileSummary(
                    surface_csv=surface_csv_path.name,
                    source_pdb="",
                    accession="",
                    allergen="",
                    surface_residue_count=0,
                    center_count=0,
                    missing_center_count=0,
                    patch_count=0,
                    status="error",
                    error=str(exc),
                )
            )

    write_file_summary_csv(summary_dir / "patch_generation_counts.csv", file_summaries)

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
