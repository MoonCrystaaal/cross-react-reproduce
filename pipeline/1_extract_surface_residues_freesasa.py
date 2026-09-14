#!/usr/bin/env python3
"""
Batch residue-level SASA extraction with FreeSASA.

Default behavior:
- Input PDB directory: ./pdb_all
- Per-file CSV output: ./1_surface_residues/surface_residue_csv
- Raw RSA output: ./1_surface_residues/surface_residue_rsa
- Master summary output: ./1_surface_residues/surface_residue_summary
- Logs: ./1_surface_residues/surface_residue_logs
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    import freesasa
except ImportError as exc:  # pragma: no cover - import guard for wrong env
    raise SystemExit(
        "freesasa module not found. Run this script in the correct environment."
    ) from exc


SCRIPT_DIR = Path(__file__).resolve().parent
RUN_ROOT = Path.cwd().resolve()
STAGE_ROOT = RUN_ROOT / "1_surface_residues"
DEFAULT_INPUT_DIR = RUN_ROOT / "pdb_all"
DEFAULT_OUTPUT_DIR = STAGE_ROOT / "surface_residue_csv"
DEFAULT_RSA_DIR = STAGE_ROOT / "surface_residue_rsa"
DEFAULT_SUMMARY_DIR = STAGE_ROOT / "surface_residue_summary"
DEFAULT_LOG_DIR = STAGE_ROOT / "surface_residue_logs"

DEFAULT_ALGORITHM = "LeeRichards"
DEFAULT_PROBE_RADIUS = 1.4
DEFAULT_SURFACE_CUTOFF = 10.0
DEFAULT_N_SLICES = 20
DEFAULT_N_POINTS = 100
DEFAULT_N_THREADS = 1


@dataclass(frozen=True)
class ResidueRow:
    source_pdb: str
    accession: str
    allergen: str
    chain_id: str
    residue_number: str
    residue_name: str
    residue_sasa_total: float
    is_surface_gt10A2: bool


@dataclass(frozen=True)
class FileSummary:
    source_pdb: str
    accession: str
    allergen: str
    residue_count: int
    surface_residue_count: int
    status: str
    error: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract residue-level SASA from PDB files with FreeSASA and mark "
            "surface residues using the paper-aligned cutoff Total SASA > 10 A^2."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Input PDB directory (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Per-file CSV output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--rsa-dir",
        type=Path,
        default=DEFAULT_RSA_DIR,
        help=f"Raw RSA output directory (default: {DEFAULT_RSA_DIR})",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=DEFAULT_SUMMARY_DIR,
        help=f"Summary CSV output directory (default: {DEFAULT_SUMMARY_DIR})",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Log output directory (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--algorithm",
        choices=["LeeRichards", "ShrakeRupley"],
        default=DEFAULT_ALGORITHM,
        help=f"FreeSASA algorithm (default: {DEFAULT_ALGORITHM})",
    )
    parser.add_argument(
        "--probe-radius",
        type=float,
        default=DEFAULT_PROBE_RADIUS,
        help=f"Solvent probe radius in angstrom (default: {DEFAULT_PROBE_RADIUS})",
    )
    parser.add_argument(
        "--surface-cutoff",
        type=float,
        default=DEFAULT_SURFACE_CUTOFF,
        help=f"Surface residue cutoff in A^2 (default: {DEFAULT_SURFACE_CUTOFF})",
    )
    parser.add_argument(
        "--n-slices",
        type=int,
        default=DEFAULT_N_SLICES,
        help=f"Lee-Richards resolution (default: {DEFAULT_N_SLICES})",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=DEFAULT_N_POINTS,
        help=f"Shrake-Rupley resolution (default: {DEFAULT_N_POINTS})",
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        default=DEFAULT_N_THREADS,
        help=f"FreeSASA thread count (default: {DEFAULT_N_THREADS})",
    )
    parser.add_argument(
        "--include-hetatm",
        action="store_true",
        help="Include HETATM records in FreeSASA structure parsing",
    )
    parser.add_argument(
        "--include-hydrogen",
        action="store_true",
        help="Include hydrogen atoms in FreeSASA structure parsing",
    )
    parser.add_argument(
        "--match",
        default="*.pdb",
        help="Glob pattern used to select input files (default: *.pdb)",
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
        help="Overwrite existing per-file CSV and RSA outputs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which files would be processed without running FreeSASA",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "surface_residue_run.log"
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


def split_filename(stem: str) -> Tuple[str, str]:
    if "_" in stem:
        accession, allergen = stem.split("_", 1)
    else:
        accession, allergen = stem, ""
    return accession, allergen


def iter_pdb_files(input_dir: Path, match: str, limit: int | None) -> List[Path]:
    files = sorted(input_dir.glob(match), key=lambda path: path.name)
    if limit is not None:
        files = files[:limit]
    return files


def build_parameters(args: argparse.Namespace) -> "freesasa.Parameters":
    algorithm = getattr(freesasa, args.algorithm)
    values = {
        "algorithm": algorithm,
        "probe-radius": args.probe_radius,
        "n-threads": args.n_threads,
    }
    if args.algorithm == "LeeRichards":
        values["n-slices"] = args.n_slices
    else:
        values["n-points"] = args.n_points
    return freesasa.Parameters(values)


def build_rsa_command(args: argparse.Namespace, pdb_path: Path) -> List[str]:
    command = ["freesasa"]
    if args.algorithm == "LeeRichards":
        command.extend(["--lee-richards", f"--resolution={args.n_slices}"])
    else:
        command.extend(["--shrake-rupley", f"--resolution={args.n_points}"])
    command.extend(
        [
            f"--probe-radius={args.probe_radius}",
            "--radii=protor",
            "--format=rsa",
        ]
    )
    if args.include_hetatm:
        command.append("--hetatm")
    if args.include_hydrogen:
        command.append("--hydrogen")
    command.append(str(pdb_path))
    return command


def read_residue_order(
    pdb_path: Path, include_hetatm: bool = False
) -> Dict[Tuple[str, str], int]:
    accepted_records = {"ATOM  "}
    if include_hetatm:
        accepted_records.add("HETATM")
    order: Dict[Tuple[str, str], int] = {}
    seen = set()
    index = 0
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if len(line) < 27:
                continue
            if line[:6] not in accepted_records:
                continue
            chain_id = line[21].strip()
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()
            if insertion_code:
                residue_number = f"{residue_number}{insertion_code}"
            key = (chain_id, residue_number)
            if key in seen:
                continue
            seen.add(key)
            order[key] = index
            index += 1
    return order


def extract_rows(
    pdb_path: Path,
    parameters: "freesasa.Parameters",
    surface_cutoff: float,
    include_hetatm: bool,
    include_hydrogen: bool,
) -> List[ResidueRow]:
    accession, allergen = split_filename(pdb_path.stem)
    structure_options = {
        "hetatm": include_hetatm,
        "hydrogen": include_hydrogen,
    }
    structure = freesasa.Structure(str(pdb_path), options=structure_options)
    result = freesasa.calc(structure, parameters)
    residue_areas = result.residueAreas()
    residue_order = read_residue_order(pdb_path, include_hetatm=include_hetatm)

    rows: List[ResidueRow] = []
    for chain_id, residues in residue_areas.items():
        for residue_number, area in residues.items():
            total = float(area.total)
            rows.append(
                ResidueRow(
                    source_pdb=pdb_path.name,
                    accession=accession,
                    allergen=allergen,
                    chain_id=chain_id,
                    residue_number=str(residue_number),
                    residue_name=area.residueType,
                    residue_sasa_total=total,
                    is_surface_gt10A2=total > surface_cutoff,
                )
            )

    def sort_key(row: ResidueRow) -> Tuple[int, str, str]:
        index = residue_order.get((row.chain_id, row.residue_number), 10**9)
        return index, row.chain_id, row.residue_number

    rows.sort(key=sort_key)
    return rows


def write_rsa_file(output_path: Path, args: argparse.Namespace, pdb_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    command = build_rsa_command(args, pdb_path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        result = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        stderr = (result.stderr or "").strip()
        detail = f": {stderr}" if stderr else ""
        raise RuntimeError(f"FreeSASA RSA export failed for {pdb_path.name}{detail}")


def write_residue_csv(output_path: Path, rows: Sequence[ResidueRow]) -> None:
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
                "residue_sasa_total",
                "is_surface_gt10A2",
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
                    f"{row.residue_sasa_total:.6f}",
                    str(row.is_surface_gt10A2),
                ]
            )


def read_existing_residue_csv(output_path: Path) -> List[ResidueRow]:
    rows: List[ResidueRow] = []
    with output_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                ResidueRow(
                    source_pdb=row["source_pdb"],
                    accession=row["accession"],
                    allergen=row["allergen"],
                    chain_id=row["chain_id"],
                    residue_number=row["residue_number"],
                    residue_name=row["residue_name"],
                    residue_sasa_total=float(row["residue_sasa_total"]),
                    is_surface_gt10A2=row["is_surface_gt10A2"].strip().lower() == "true",
                )
            )
    return rows


def write_master_csv(summary_dir: Path, rows: Sequence[ResidueRow]) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    output_path = summary_dir / "surface_residue_master.csv"
    write_residue_csv(output_path, rows)


def write_file_summary_csv(summary_dir: Path, summaries: Sequence[FileSummary]) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    output_path = summary_dir / "surface_residue_counts.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "source_pdb",
                "accession",
                "allergen",
                "residue_count",
                "surface_residue_count",
                "status",
                "error",
            ]
        )
        for summary in summaries:
            writer.writerow(
                [
                    summary.source_pdb,
                    summary.accession,
                    summary.allergen,
                    summary.residue_count,
                    summary.surface_residue_count,
                    summary.status,
                    summary.error,
                ]
            )


def main() -> int:
    args = parse_args()
    configure_logging(args.log_dir)

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    rsa_dir = args.rsa_dir.resolve()
    summary_dir = args.summary_dir.resolve()

    if not input_dir.exists():
        logging.error("Input directory not found: %s", input_dir)
        return 1

    pdb_files = iter_pdb_files(input_dir, args.match, args.limit)
    if not pdb_files:
        logging.error("No PDB files matched pattern '%s' in %s", args.match, input_dir)
        return 1

    logging.info("Script directory: %s", SCRIPT_DIR)
    logging.info("Run root: %s", RUN_ROOT)
    logging.info("Stage root: %s", STAGE_ROOT)
    logging.info("Input directory: %s", input_dir)
    logging.info("CSV output directory: %s", output_dir)
    logging.info("RSA output directory: %s", rsa_dir)
    logging.info("Matched PDB files: %d", len(pdb_files))
    logging.info(
        "Parameters: algorithm=%s probe_radius=%.3f surface_cutoff=%.3f",
        args.algorithm,
        args.probe_radius,
        args.surface_cutoff,
    )
    logging.info(
        "Structure options: include_hetatm=%s include_hydrogen=%s",
        args.include_hetatm,
        args.include_hydrogen,
    )

    if args.dry_run:
        for pdb_path in pdb_files:
            logging.info("[DRY RUN] %s", pdb_path.name)
        return 0

    parameters = build_parameters(args)
    all_rows: List[ResidueRow] = []
    summaries: List[FileSummary] = []

    for index, pdb_path in enumerate(pdb_files, start=1):
        accession, allergen = split_filename(pdb_path.stem)
        stem = sanitize_stem(pdb_path.stem)
        output_path = output_dir / f"{stem}.csv"
        rsa_output_path = rsa_dir / f"{stem}.rsa"

        if output_path.exists() and rsa_output_path.exists() and not args.overwrite:
            existing_rows = read_existing_residue_csv(output_path)
            all_rows.extend(existing_rows)
            surface_count = sum(row.is_surface_gt10A2 for row in existing_rows)
            logging.info(
                "[%d/%d] Skipping existing outputs: %s",
                index,
                len(pdb_files),
                stem,
            )
            summaries.append(
                FileSummary(
                    source_pdb=pdb_path.name,
                    accession=accession,
                    allergen=allergen,
                    residue_count=len(existing_rows),
                    surface_residue_count=surface_count,
                    status="skipped_existing",
                    error="",
                )
            )
            continue

        if (output_path.exists() or rsa_output_path.exists()) and not args.overwrite:
            logging.info(
                "[%d/%d] Regenerating incomplete outputs: %s",
                index,
                len(pdb_files),
                stem,
            )

        try:
            rows = extract_rows(
                pdb_path=pdb_path,
                parameters=parameters,
                surface_cutoff=args.surface_cutoff,
                include_hetatm=args.include_hetatm,
                include_hydrogen=args.include_hydrogen,
            )
            write_residue_csv(output_path, rows)
            write_rsa_file(rsa_output_path, args, pdb_path)
            all_rows.extend(rows)
            surface_count = sum(row.is_surface_gt10A2 for row in rows)
            summaries.append(
                FileSummary(
                    source_pdb=pdb_path.name,
                    accession=accession,
                    allergen=allergen,
                    residue_count=len(rows),
                    surface_residue_count=surface_count,
                    status="ok",
                    error="",
                )
            )
            logging.info(
                "[%d/%d] %s -> residues=%d surface=%d",
                index,
                len(pdb_files),
                pdb_path.name,
                len(rows),
                surface_count,
            )
        except Exception as exc:  # pragma: no cover - file-level failure path
            logging.exception("Failed to process %s", pdb_path.name)
            summaries.append(
                FileSummary(
                    source_pdb=pdb_path.name,
                    accession=accession,
                    allergen=allergen,
                    residue_count=0,
                    surface_residue_count=0,
                    status="error",
                    error=str(exc),
                )
            )

    write_file_summary_csv(summary_dir, summaries)
    if all_rows:
        write_master_csv(summary_dir, all_rows)

    ok_count = sum(summary.status == "ok" for summary in summaries)
    error_count = sum(summary.status == "error" for summary in summaries)
    skipped_count = sum(summary.status == "skipped_existing" for summary in summaries)
    logging.info(
        "Done. ok=%d error=%d skipped=%d total=%d",
        ok_count,
        error_count,
        skipped_count,
        len(summaries),
    )
    return 0 if error_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
