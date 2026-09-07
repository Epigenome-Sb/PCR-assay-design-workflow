#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""
02_in_silico_validation.py

Mode-aware in silico validation for VarVAMP assay designs.

Branches
--------
SINGLE: validate assigned primer pairs with MFEprimer.
QPCR: validate LEFT/RIGHT with MFEprimer, extract retained PCR amplicons,
expand any existing IUPAC ambiguity in the VarVAMP probe into concrete
A/C/G/T variants, and search those variants with blastn-short. No new
ambiguity is proposed.
TILED: validate assigned primer pairs with MFEprimer; no probe BLAST.

Input/output architecture
-------------------------
The validation workflow consumes two independent sources:
1. a completed Workflow 01 assay-design run, selected through its
   assay_design_manifest.json;
2. a larger validation FASTA, normally stored under data/validation/.

Each validation execution receives an isolated run ID. By default:

    work/<project>/validation/<run_id>/
        assay_inputs/       reconstructible primer/probe FASTA files
        database/           filtered/formatted FASTA + MFEprimer indexes

    results/<project>/validation/<run_id>/
        inputs/             provenance/handoff records
        mfeprimer/          raw + processed primer-pair validation results
        probe_blast/        qPCR probe BLAST results
        summary/            cross-assay summary tables
        validation_manifest.json
        validation.log

No validation_inputs/ directory is used.

Normal terminal mode is concise. Use --verbose to show commands and detailed
paths. Every external command is recorded in the run-specific validation.log.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from contextlib import redirect_stdout
from io import StringIO
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from math import ceil
from pathlib import Path

IUPAC_DNA = set("ACGTRYSWKMBDHVN")

VERBOSE = False
VALIDATION_LOG: Path | None = None

def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mode-aware in silico validation of VarVAMP SINGLE, QPCR or TILED "
            "assay designs with MFEprimer and, for QPCR, BLAST+ probe analysis."
        )
    )
    parser.add_argument(
        "--design-manifest",
        type=Path,
        default=None,
        help=(
            "Workflow 01 assay_design_manifest.json. When omitted, completed "
            "design runs under results/<project>/design/<run_id>/ are listed."
        ),
    )
    parser.add_argument(
        "--validation-db",
        type=Path,
        default=None,
        help=(
            "Validation FASTA. When omitted, FASTA files under data/validation/ "
            "are listed."
        ),
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Optional validation run identifier; default is a local timestamp.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help=(
            "Optional exact work directory override. Default: "
            "work/<project>/validation/<run_id>/"
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help=(
            "Optional exact results directory override. Default: "
            "results/<project>/validation/<run_id>/"
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show external commands and detailed paths in the terminal.",
    )
    return parser.parse_args()

def log_line(message: str) -> None:
    if VALIDATION_LOG is None:
        return
    VALIDATION_LOG.parent.mkdir(parents=True, exist_ok=True)
    with VALIDATION_LOG.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")

def compact_path(path: Path | str) -> str:
    path_obj = Path(path)
    try:
        return str(path_obj.resolve().relative_to(Path.cwd().resolve()))
    except Exception:
        return str(path_obj)


def portable_project_path(path: Path | str) -> str:
    """
    Serialize paths for manifests and provenance tables.

    Files located inside the current project are stored relative to the
    repository root. Truly external files remain absolute.
    """
    resolved = Path(path).expanduser().resolve()
    repository_root = Path.cwd().resolve()

    try:
        return resolved.relative_to(repository_root).as_posix()
    except ValueError:
        return str(resolved)


def optional_portable_project_path(value: object) -> str:
    """Serialize an optional path while preserving empty/not-applicable values."""
    if value is None or str(value) == "":
        return ""
    return portable_project_path(Path(str(value)))


def _legacy_project_relative_candidate(path: Path) -> Path | None:
    """
    Recover a project-relative suffix from an old absolute path.

    This is mainly for manifests created before portable paths were introduced,
    for example:
        /app/results/<project>/design/<run_id>/varvamp
        /home/user/repo/results/<project>/design/<run_id>/varvamp
    """
    parts = path.parts
    for marker in ("results", "work", "data"):
        if marker in parts:
            index = parts.index(marker)
            return Path(*parts[index:])
    return None


def resolve_manifest_path(value: str | Path, manifest_path: Path) -> Path:
    """
    Resolve a path stored in a Workflow 01 manifest.

    Resolution order:
    1. Existing absolute path.
    2. Legacy absolute project path remapped under the current repository root.
    3. Repository-relative path (schema v2 portable path policy).
    4. Legacy manifest-directory-relative path.

    This makes a design run portable between native execution and Docker while
    retaining compatibility with older manifests that stored /app/... paths.
    """
    repository_root = Path.cwd().resolve()
    path = Path(value).expanduser()

    if path.is_absolute():
        resolved = path.resolve()
        if resolved.exists():
            return resolved

        legacy_relative = _legacy_project_relative_candidate(resolved)
        if legacy_relative is not None:
            remapped = (repository_root / legacy_relative).resolve()
            if remapped.exists():
                return remapped

        return resolved

    repository_candidate = (repository_root / path).resolve()
    if repository_candidate.exists():
        return repository_candidate

    manifest_candidate = (manifest_path.parent / path).resolve()
    if manifest_candidate.exists():
        return manifest_candidate

    return repository_candidate


def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))

def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 48 - len(label))
    print(f"{label} {dots} {status}")


# ---------------------------------------------------------------------------
# Generic terminal helpers
# ---------------------------------------------------------------------------

def ask_yes_no(question: str) -> bool:
    while True:
        value = input(f"{question} [y/n]: ").strip().lower()
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please answer y or n.")


def ask_choice(title: str, choices: list[str]) -> int:
    print(f"\n{title}")
    print("=" * len(title))
    for idx, item in enumerate(choices, start=1):
        print(f"{idx}. {item}")

    while True:
        raw = input("Select an option by number: ").strip()
        try:
            selected = int(raw)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue
        if 1 <= selected <= len(choices):
            return selected
        print("Selection outside the available range.")


def ask_int(label: str, minimum: int, maximum: int) -> int:
    while True:
        raw = input(f"{label} [{minimum}-{maximum}]: ").strip()
        try:
            value = int(raw)
        except ValueError:
            print("Please enter an integer.")
            continue
        if minimum <= value <= maximum:
            return value
        print(f"Value must be between {minimum} and {maximum}.")


def ask_float(label: str, minimum: float, maximum: float) -> float:
    while True:
        raw = input(f"{label} [{minimum}-{maximum}]: ").strip()
        try:
            value = float(raw)
        except ValueError:
            print("Please enter a numeric value.")
            continue
        if minimum <= value <= maximum:
            return value
        print(f"Value must be between {minimum} and {maximum}.")


def run_command(command: list[str]) -> None:
    command_text = " ".join(command)
    log_line("COMMAND: " + command_text)
    if VERBOSE:
        print("\nCommand:", command_text)
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: " + command_text
        ) from error


def run_command_to_file(command: list[str], output_file: Path) -> None:
    command_text = " ".join(command) + " > " + str(output_file)
    log_line("COMMAND: " + command_text)
    if VERBOSE:
        print("\nCommand:", command_text)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_file.open("w", encoding="utf-8") as handle:
            subprocess.run(command, check=True, stdout=handle)
    except subprocess.CalledProcessError as error:
        output_file.unlink(missing_ok=True)
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: " + command_text
        ) from error



# ---------------------------------------------------------------------------
# Run/provenance helpers
# ---------------------------------------------------------------------------

def sanitize_identifier(value: str, label: str) -> str:
    """Convert a project/run label into a safe directory identifier."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError(f"{label} is empty or invalid.")
    return cleaned


def generate_run_id() -> str:
    """Generate a local timestamp suitable for one validation run."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path) -> str:
    """Return SHA-256 for provenance tracking."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def executable_version(executable: str) -> str:
    """Best-effort capture of a command-line tool version."""
    for command in (
        [executable, "--version"],
        [executable, "-version"],
        [executable, "version"],
    ):
        try:
            result = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except Exception:
            continue
        output = (result.stdout or result.stderr or "").strip()
        if output:
            return output.splitlines()[0][:300]
    return "unknown"


def load_design_manifest(manifest_path: Path) -> dict[str, object]:
    """Load and validate the Workflow 01 handoff manifest."""
    manifest_path = manifest_path.expanduser().resolve()
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Design manifest not found: {manifest_path}")

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"Invalid JSON in design manifest {manifest_path}: {error}"
        ) from error

    if manifest.get("workflow_stage") not in {None, "design"}:
        raise ValueError(
            f"Manifest is not an assay-design manifest: {manifest_path}"
        )

    project_name = manifest.get("project_name")
    if not isinstance(project_name, str) or not project_name.strip():
        raise ValueError("Design manifest does not contain a valid project_name.")

    varvamp = manifest.get("varvamp")
    if not isinstance(varvamp, dict):
        raise ValueError("Design manifest does not contain a VarVAMP section.")

    mode = str(varvamp.get("mode") or "").lower()
    if mode not in {"single", "qpcr", "tiled"}:
        raise ValueError(
            "Design manifest does not identify a completed SINGLE/QPCR/TILED run."
        )

    result_value = varvamp.get("result_dir")
    if not result_value:
        result_directories = manifest.get("result_directories")
        if isinstance(result_directories, dict):
            result_value = result_directories.get("varvamp")
    if not result_value:
        raise ValueError("Design manifest does not contain a VarVAMP result directory.")

    result_dir = resolve_manifest_path(str(result_value), manifest_path)
    if not result_dir.is_dir():
        raise FileNotFoundError(
            "VarVAMP result directory referenced by the design manifest does not "
            f"exist: {result_dir}"
        )

    manifest["_manifest_path"] = str(manifest_path)
    manifest["_resolved_varvamp_result_dir"] = str(result_dir)
    manifest["_resolved_mode"] = mode
    return manifest


def find_design_manifests() -> list[Path]:
    """Find run-level Workflow 01 manifests without selecting nested copies."""
    root = Path.cwd() / "results"
    if not root.is_dir():
        return []
    candidates = [
        path.resolve()
        for path in root.glob("*/design/*/assay_design_manifest.json")
        if path.is_file()
    ]
    return sorted(
        set(candidates),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )


def choose_design_run(
    requested_manifest: Path | None,
) -> tuple[str, Path, Path, dict[str, object], str, str]:
    """
    Select a Workflow 01 run through assay_design_manifest.json.

    Returns mode, VarVAMP result directory, manifest path, manifest data,
    project name and design run ID.
    """
    if requested_manifest is not None:
        paths = [requested_manifest.expanduser().resolve()]
    else:
        paths = find_design_manifests()
        if not paths:
            raise FileNotFoundError(
                "No assay_design_manifest.json was found under "
                "results/<project>/design/<run_id>/. Run Workflow 01 first or "
                "provide --design-manifest."
            )

    valid: list[tuple[Path, dict[str, object]]] = []
    errors: list[str] = []
    for path in paths:
        try:
            valid.append((path, load_design_manifest(path)))
        except (FileNotFoundError, ValueError) as error:
            errors.append(f"{path}: {error}")

    if not valid:
        raise RuntimeError(
            "No usable Workflow 01 design manifest was found.\n" + "\n".join(errors)
        )

    print("\n1. Assay-design run")
    print("===================")

    if requested_manifest is None:
        for idx, (path, manifest) in enumerate(valid, start=1):
            mode = str(manifest["_resolved_mode"]).upper()
            project = str(manifest["project_name"])
            design_run = str(manifest.get("run_id") or path.parent.name)
            design_input = Path(str(manifest.get("input_fasta", "unknown"))).name
            print(
                f"{idx}. {project} | {mode} | design run {design_run} | "
                f"input {design_input}"
            )

        while True:
            raw = input("Select the assay-design run by number: ").strip()
            try:
                selected = int(raw)
            except ValueError:
                print("Please enter one of the displayed numbers.")
                continue
            if 1 <= selected <= len(valid):
                manifest_path, manifest = valid[selected - 1]
                break
            print("Selection outside the available range.")
    else:
        manifest_path, manifest = valid[0]

    mode = str(manifest["_resolved_mode"])
    result_dir = Path(str(manifest["_resolved_varvamp_result_dir"]))
    project_name = sanitize_identifier(str(manifest["project_name"]), "Project name")
    design_run_id = sanitize_identifier(
        str(manifest.get("run_id") or manifest_path.parent.name),
        "Design run ID",
    )

    print(f"Project     : {project_name}")
    print(f"Design run  : {design_run_id}")
    print(f"Mode        : {mode.upper()}")
    print(f"VarVAMP dir : {compact_path(result_dir)}")
    print(f"Manifest    : {compact_path(manifest_path)}")

    if requested_manifest is None and not ask_yes_no("Use this design run?"):
        raise RuntimeError("Assay-design run selection cancelled.")

    return mode, result_dir, manifest_path, manifest, project_name, design_run_id


def resolve_validation_run_directories(
    cli: argparse.Namespace,
    project_name: str,
) -> tuple[str, Path, Path]:
    """Create an isolated work/results namespace for one validation run."""
    requested = (
        sanitize_identifier(cli.run_id, "Run ID")
        if cli.run_id is not None
        else generate_run_id()
    )

    if cli.workdir is not None or cli.results is not None:
        workdir = (
            cli.workdir.expanduser().resolve()
            if cli.workdir is not None
            else (Path("work") / project_name / "validation" / requested).resolve()
        )
        results_dir = (
            cli.results.expanduser().resolve()
            if cli.results is not None
            else (Path("results") / project_name / "validation" / requested).resolve()
        )
        if results_dir.exists() and any(results_dir.iterdir()):
            raise RuntimeError(
                f"Results directory is not empty; refusing to overwrite: {results_dir}"
            )
        return requested, workdir, results_dir

    candidate = requested
    suffix = 1
    while True:
        workdir = (Path("work") / project_name / "validation" / candidate).resolve()
        results_dir = (Path("results") / project_name / "validation" / candidate).resolve()
        if not workdir.exists() and not results_dir.exists():
            return candidate, workdir, results_dir
        suffix += 1
        candidate = f"{requested}_{suffix:02d}"


def write_assay_manifest(
    path: Path,
    assay_mode: str,
    schemes: dict[str, dict[str, str]],
) -> None:
    """Save the exact oligos that entered validation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "mode", "scheme", "left_sequence", "probe_sequence",
                "right_sequence", "left_length", "probe_length", "right_length",
            ]
        )
        for scheme in sorted(schemes):
            left = schemes[scheme]["LEFT"]
            right = schemes[scheme]["RIGHT"]
            probe = schemes[scheme].get("PROBE", "")
            writer.writerow(
                [
                    assay_mode, scheme, left, probe, right,
                    len(left), len(probe) if probe else "", len(right),
                ]
            )


def write_validation_database_summary(
    path: Path,
    *,
    source_database: Path,
    source_stats: dict[str, float | int],
    filter_description: str,
    retained_database: Path,
    retained_stats: dict[str, float | int],
    formatted_database: Path,
) -> None:
    """Write a compact provenance table for the validation database."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        rows = [
            ("source_fasta", portable_project_path(source_database)),
            ("source_sha256", sha256_file(source_database)),
            ("source_sequences", int(source_stats["sequences"])),
            ("source_length_range", f"{source_stats['min_length']}-{source_stats['max_length']}"),
            ("length_filter", filter_description),
            ("retained_fasta", portable_project_path(retained_database)),
            ("retained_sequences", int(retained_stats["sequences"])),
            ("retained_length_range", f"{retained_stats['min_length']}-{retained_stats['max_length']}"),
            ("formatted_mfeprimer_fasta", portable_project_path(formatted_database)),
        ]
        writer.writerows(rows)


def write_validation_manifest(
    path: Path,
    *,
    project_name: str,
    run_id: str,
    design_manifest_path: Path,
    design_manifest: dict[str, object],
    design_run_id: str,
    assay_mode: str,
    varvamp_result_dir: Path,
    workdir: Path,
    results_dir: Path,
    source_database: Path,
    raw_stats: dict[str, float | int],
    length_filter_description: str,
    retained_database: Path,
    retained_stats: dict[str, float | int],
    formatted_database: Path,
    index_info: dict[str, object],
    selected_schemes: list[str],
    mfeprimer_parameters: dict[str, float | int] | None,
    run_result: dict[str, object],
) -> None:
    """Write a machine-readable provenance record for Workflow 02."""
    index_files = [
        portable_project_path(Path(path))
        for path in index_info.get("files", [])
    ]
    manifest = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "workflow_stage": "validation",
        "path_base": "repository_root",
        "path_policy": (
            "repository-relative when inside the project; "
            "absolute only for external paths"
        ),
        "project_name": project_name,
        "run_id": run_id,
        "parent_design": {
            "manifest": portable_project_path(design_manifest_path),
            "run_id": design_run_id,
            "mode": assay_mode,
            "varvamp_result_dir": portable_project_path(varvamp_result_dir),
            "design_input_fasta": design_manifest.get("input_fasta"),
            "design_input_sha256": design_manifest.get("input_sha256"),
        },
        "input_role": "validation_database",
        "validation_database": {
            "source_fasta": portable_project_path(source_database),
            "source_sha256": sha256_file(source_database),
            "raw_statistics": raw_stats,
            "length_filter": length_filter_description,
            "retained_fasta": portable_project_path(retained_database),
            "retained_statistics": retained_stats,
            "mfeprimer_formatted_fasta": portable_project_path(formatted_database),
            "index": {
                "status": index_info.get("status"),
                "layout": index_info.get("detail"),
                "files": index_files,
            },
        },
        "workdir": portable_project_path(workdir),
        "results_dir": portable_project_path(results_dir),
        "result_directories": {
            "inputs": portable_project_path(results_dir / "inputs"),
            "mfeprimer": portable_project_path(results_dir / "mfeprimer"),
            "probe_blast": portable_project_path(results_dir / "probe_blast"),
            "summary": portable_project_path(results_dir / "summary"),
        },
        "selected_assays": selected_schemes,
        "mfeprimer_parameters": (
            {
                "mode": "defaults",
                "min_size": 0,
                "max_size": 2000,
                "tm_cutoff": 30,
            }
            if mfeprimer_parameters is None
            else {"mode": "custom", **mfeprimer_parameters}
        ),
        "coverage_denominator_sequences": int(retained_stats["sequences"]),
        "product_filter": {
            "target_only": bool(run_result["exclude_self_priming"]),
            "description": (
                "TARGET only (SELF_PRIMING excluded)"
                if run_result["exclude_self_priming"]
                else "TARGET + SELF_PRIMING retained"
            ),
        },
        "probe_validation": {
            "applicable": assay_mode == "qpcr",
            "executed": bool(run_result["probe_validation_enabled"]),
        },
        "summary_files": {
            "coverage_results": portable_project_path(Path(run_result["coverage_results"])),
            "mfeprimer_runs": portable_project_path(Path(run_result["mfeprimer_runs"])),
            "secondary_structure_summary": portable_project_path(Path(run_result["secondary_summary"])),
            "secondary_structure_details": portable_project_path(Path(run_result["secondary_details"])),
        },
        "tool_versions": {
            "seqkit": executable_version("seqkit"),
            "mfeprimer": executable_version("mfeprimer"),
            "blastn": executable_version("blastn") if assay_mode == "qpcr" else "not applicable",
        },
        "workflow_log": None if VALIDATION_LOG is None else portable_project_path(VALIDATION_LOG),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: list[tuple[str, str]] = []
    current: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current, seq_parts
        if current is None:
            return

        sequence = "".join(seq_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {current}")

        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide character(s) in {current}: "
                + ", ".join(invalid)
            )

        records.append((current, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                current = line[1:].split()[0]
                seq_parts = []
            else:
                if current is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )
                seq_parts.append(line)

    flush()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_fasta_with_full_headers(path: Path) -> list[tuple[str, str]]:
    """
    Read FASTA records while preserving the complete header text.

    This is used for database filtering/formatting checks so descriptive
    headers are not silently shortened.
    """
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal header, seq_parts
        if header is None:
            return

        sequence = "".join(seq_parts).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {header}")

        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise ValueError(
                f"Invalid nucleotide character(s) in {header}: "
                + ", ".join(invalid)
            )

        records.append((header, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )
                seq_parts.append(line)

    flush()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def write_fasta_with_full_headers(
    path: Path,
    records: list[tuple[str, str]],
) -> None:
    """Write FASTA records while preserving complete header descriptions."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def fasta_statistics(path: Path) -> dict[str, float | int]:
    records = read_fasta_with_full_headers(path)
    lengths = [len(seq) for _, seq in records]
    sorted_lengths = sorted(lengths)
    n = len(sorted_lengths)

    if n % 2:
        median_length = float(sorted_lengths[n // 2])
    else:
        median_length = (
            sorted_lengths[n // 2 - 1] + sorted_lengths[n // 2]
        ) / 2.0

    ambiguous = sum(
        sum(base not in {"A", "C", "G", "T"} for base in seq)
        for _, seq in records
    )

    return {
        "sequences": len(records),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "median_length": median_length,
        "ambiguous": ambiguous,
    }


def print_fasta_statistics(path: Path, title: str) -> None:
    stats = fasta_statistics(path)
    print(f"\n{title}")
    print("=" * len(title))
    print(f"FASTA:                 {path}")
    print(f"Sequences:             {stats['sequences']}")
    print(
        f"Sequence length range: {stats['min_length']}-"
        f"{stats['max_length']} nt"
    )
    print(f"Mean sequence length:  {stats['mean_length']:.2f} nt")
    print(f"Median sequence length:{stats['median_length']:.2f} nt")
    print(f"Ambiguous bases:       {stats['ambiguous']}")


def choose_sequence_length_filter(
    database: Path,
    database_workdir: Path,
) -> tuple[Path, str, dict[str, float | int]]:
    """
    Decide which sequence-length set enters MFEprimer preparation.

    The "near-complete" option is an operational length criterion relative
    to the longest sequence in the selected database; it does not infer
    biological completeness from annotation.
    """
    before_stats = fasta_statistics(database)
    max_length = int(before_stats["max_length"])

    print("\n3. Sequence-length filtering")
    print("============================")
    print(
        f"Current database: {before_stats['sequences']} sequences; "
        f"length range {before_stats['min_length']}-{before_stats['max_length']} nt."
    )

    mode = ask_choice(
        "How should the validation database be prepared?",
        [
            "Keep all sequences.",
            "Keep complete / nearly complete sequences using a relative-length threshold.",
            "Apply a custom minimum sequence length.",
        ],
    )

    threshold: int | None = None
    strategy: str

    if mode == 1:
        strategy = "keep all sequences"
        retained_database = database
    elif mode == 2:
        print(
            "\nNear-complete filtering is defined here by sequence length "
            "relative to the longest sequence in this database."
        )
        relative_mode = ask_choice(
            "Near-complete length threshold",
            [
                "Keep sequences >=90% of the maximum sequence length.",
                "Keep sequences >=95% of the maximum sequence length.",
                "Enter a custom percentage of the maximum sequence length.",
            ],
        )
        if relative_mode == 1:
            percentage = 90.0
        elif relative_mode == 2:
            percentage = 95.0
        else:
            percentage = ask_float(
                "Minimum relative length (% of maximum)",
                50.0,
                100.0,
            )
        threshold = ceil(max_length * percentage / 100.0)
        strategy = (
            f"near-complete length >= {percentage:g}% of maximum "
            f"({threshold} nt)"
        )
    else:
        threshold = ask_int("Minimum sequence length (nt)", 1, max_length)
        strategy = f"minimum sequence length >= {threshold} nt"

    if threshold is not None:
        source_records = read_fasta_with_full_headers(database)
        retained_records = [
            (header, sequence)
            for header, sequence in source_records
            if len(sequence) >= threshold
        ]
        if not retained_records:
            raise RuntimeError(
                "The selected length filter removed every sequence. "
                "Choose a less restrictive threshold."
            )
        database_workdir.mkdir(parents=True, exist_ok=True)
        safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
        retained_database = database_workdir / f"{safe_stem}_min{threshold}.fasta"
        write_fasta_with_full_headers(retained_database, retained_records)

    after_stats = fasta_statistics(retained_database)
    before_n = int(before_stats["sequences"])
    after_n = int(after_stats["sequences"])
    removed = before_n - after_n
    retained_pct = 100.0 * after_n / before_n if before_n else 0.0
    removed_pct = 100.0 * removed / before_n if before_n else 0.0

    print("\nSequence filtering results")
    print("--------------------------")
    print(f"Strategy           : {strategy}")
    print(f"Sequences before   : {before_n}")
    print(f"Sequences retained : {after_n} ({retained_pct:.2f} %)")
    print(f"Sequences removed  : {removed} ({removed_pct:.2f} %)")
    print(f"Retained length    : {after_stats['min_length']}-{after_stats['max_length']} nt")
    print(f"Mean length        : {after_stats['mean_length']:.2f} nt")
    print(f"Median length      : {after_stats['median_length']:.2f} nt")
    print(f"Ambiguous bases    : {after_stats['ambiguous']}")
    if threshold is not None:
        print(f"Filtered FASTA     : {compact_path(retained_database)}")
    else:
        print("Filtered FASTA     : not created (all sequences retained)")

    return retained_database.resolve(), strategy, after_stats

# ---------------------------------------------------------------------------
# VarVAMP result-mode detection
# ---------------------------------------------------------------------------

def find_varvamp_result_directories() -> list[tuple[str, Path]]:
    """Compatibility helper backed by Workflow 01 run manifests."""
    found: list[tuple[str, Path]] = []
    for manifest_path in find_design_manifests():
        try:
            manifest = load_design_manifest(manifest_path)
        except (FileNotFoundError, ValueError):
            continue
        found.append(
            (
                str(manifest["_resolved_mode"]),
                Path(str(manifest["_resolved_varvamp_result_dir"])),
            )
        )
    return found

def choose_varvamp_result() -> tuple[str, Path]:
    """Compatibility wrapper; new code should use choose_design_run()."""
    mode, result_dir, _, _, _, _ = choose_design_run(None)
    return mode, result_dir

def _pick_column(fieldnames: list[str], tokens: tuple[str, ...]) -> str | None:
    lowered = {name: name.lower() for name in fieldnames}
    for token in tokens:
        for name, lower in lowered.items():
            if lower == token:
                return name
    for token in tokens:
        for name, lower in lowered.items():
            if token in lower:
                return name
    return None


def read_varvamp_primer_table(primer_tsv: Path) -> dict[str, str]:
    if not primer_tsv.is_file():
        raise FileNotFoundError(f"VarVAMP primer table not found: {primer_tsv}")
    with primer_tsv.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise RuntimeError(f"No header was found in {primer_tsv.name}.")
        id_col = _pick_column(fieldnames, ("primer", "primer_name", "name", "id"))
        seq_col = _pick_column(
            fieldnames,
            ("sequence", "ambiguous_sequence", "primer_sequence", "seq"),
        )
        if seq_col is None:
            raise RuntimeError(
                "Could not identify the primer-sequence column in primer.tsv. "
                f"Columns found: {', '.join(fieldnames)}"
            )
        rows = list(reader)
    primers: dict[str, str] = {}
    for row_index, row in enumerate(rows, 1):
        seq = str(row.get(seq_col, "")).strip().upper()
        if not seq or not set(seq) <= IUPAC_DNA:
            continue
        primer_id = str(row.get(id_col, "")).strip() if id_col else ""
        if not primer_id:
            primer_id = f"primer_{row_index}"
        primers[primer_id] = seq
    if not primers:
        raise RuntimeError(f"No usable primer sequences were parsed from {primer_tsv}.")
    return primers


def read_varvamp_assignments(
    assignment_file: Path,
    primer_ids: set[str],
) -> list[tuple[str, str, str]]:
    if not assignment_file.is_file():
        raise FileNotFoundError(
            "VarVAMP primer-to-amplicon assignment file not found: "
            f"{assignment_file}"
        )
    pairs: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    with assignment_file.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            matched = [
                pid for pid in primer_ids
                if re.search(
                    rf"(?<![A-Za-z0-9_.-]){re.escape(pid)}(?![A-Za-z0-9_.-])",
                    line,
                )
            ]
            matched = sorted(set(matched), key=lambda pid: line.find(pid))
            if len(matched) < 2:
                continue
            left_id, right_id = matched[:2]
            if (left_id, right_id) in seen:
                continue
            seen.add((left_id, right_id))
            first = line.split("\t")[0].strip()
            if first and first not in {left_id, right_id} and first.lower() not in {"amplicon", "primer"}:
                scheme = re.sub(r"[^A-Za-z0-9._-]+", "_", first).strip("_")
            else:
                scheme = f"amplicon_{len(pairs)+1}"
            pairs.append((scheme or f"amplicon_{len(pairs)+1}", left_id, right_id))
    if not pairs:
        raise RuntimeError(
            "No primer pairs could be recovered from "
            "primer_to_amplicon_assignments.tabular."
        )
    return pairs


def prepare_single_or_tiled_files(
    result_dir: Path,
    assay_inputs_dir: Path,
    mode: str,
) -> dict[str, dict[str, str]]:
    primers = read_varvamp_primer_table(result_dir / "primer.tsv")
    pairs = read_varvamp_assignments(
        result_dir / "primer_to_amplicon_assignments.tabular",
        set(primers),
    )
    primer_dir = assay_inputs_dir / "primer_pairs"
    scheme_root = assay_inputs_dir / "schemes"
    primer_dir.mkdir(parents=True, exist_ok=True)
    scheme_root.mkdir(parents=True, exist_ok=True)
    schemes: dict[str, dict[str, str]] = {}
    for index, (raw_scheme, left_id, right_id) in enumerate(pairs, 1):
        scheme = re.sub(r"[^A-Za-z0-9._-]+", "_", raw_scheme).strip("_")
        if not scheme or scheme in schemes:
            scheme = f"amplicon_{index}"
        left_seq, right_seq = primers[left_id], primers[right_id]
        schemes[scheme] = {"LEFT": left_seq, "RIGHT": right_seq}
        records = [(f"{scheme}_LEFT", left_seq), (f"{scheme}_RIGHT", right_seq)]
        write_fasta(primer_dir / f"{scheme}_primers.fasta", records)
        write_fasta(scheme_root / scheme / "primers.fasta", records)
    section("VarVAMP assay preparation")
    print(f"Mode         : {mode.upper()}")
    print(f"Primer pairs : {len(schemes)}")
    progress("Preparing primer-pair validation files")
    return schemes


def find_qpcr_oligos(result_dir: Path) -> Path:
    candidate = result_dir / "oligos.fasta"
    if candidate.is_file():
        return candidate.resolve()
    matches = sorted(result_dir.glob("*oligo*.fasta"))
    if matches:
        return matches[0].resolve()
    raise FileNotFoundError(f"No oligos.fasta was found inside {result_dir}.")

# ---------------------------------------------------------------------------
# VarVAMP oligo preparation
# ---------------------------------------------------------------------------

def find_oligos_fasta() -> Path:
    cwd = Path.cwd()

    preferred_names = (
        "oligos.fasta",
        "oligos.fa",
        "qpcr_oligos.fasta",
        "varvamp_oligos.fasta",
    )

    for name in preferred_names:
        candidate = cwd / name
        if candidate.is_file():
            return candidate.resolve()

    candidates: list[Path] = []

    for root in (cwd / "results", cwd / "work"):
        if not root.is_dir():
            continue
        for pattern in ("*oligo*.fasta", "*oligo*.fa"):
            candidates.extend(root.rglob(pattern))

    candidates = sorted({p.resolve() for p in candidates if p.is_file()})

    if not candidates:
        raise FileNotFoundError(
            "No oligos FASTA file was found automatically.\n"
            "Use a Workflow 01 assay-design manifest that points to the VarVAMP result directory."
        )

    if len(candidates) == 1:
        return candidates[0]

    print("\nSeveral oligos FASTA files were found:")
    for idx, path in enumerate(candidates, start=1):
        try:
            display = path.relative_to(cwd)
        except ValueError:
            display = path
        print(f"{idx}. {display}")

    while True:
        raw = input("Select the oligos FASTA file by number: ").strip()
        try:
            idx = int(raw)
        except ValueError:
            print("Please enter one of the displayed numbers.")
            continue

        if 1 <= idx <= len(candidates):
            return candidates[idx - 1]

        print("Selection outside the available range.")


def parse_varvamp_header(header: str) -> tuple[str, str]:
    match = re.fullmatch(
        r"(.+)_(LEFT|PROBE|RIGHT)",
        header,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError(
            f"Unrecognized VarVAMP header: {header!r}. "
            "Expected _LEFT, _PROBE or _RIGHT."
        )

    return match.group(1), match.group(2).upper()


def group_schemes(
    records: list[tuple[str, str]],
) -> dict[str, dict[str, str]]:
    schemes: dict[str, dict[str, str]] = {}

    for header, sequence in records:
        scheme, kind = parse_varvamp_header(header)
        schemes.setdefault(scheme, {})

        if kind in schemes[scheme]:
            raise ValueError(f"Duplicate {kind} for scheme {scheme}")

        schemes[scheme][kind] = sequence

    return schemes


def prepare_varvamp_files(
    oligos_file: Path,
    assay_inputs_dir: Path,
) -> dict[str, dict[str, str]]:
    records = read_fasta(oligos_file)
    schemes = group_schemes(records)

    primer_dir = assay_inputs_dir / "primer_pairs"
    probe_dir = assay_inputs_dir / "probes"
    scheme_root = assay_inputs_dir / "schemes"

    primer_dir.mkdir(parents=True, exist_ok=True)
    probe_dir.mkdir(parents=True, exist_ok=True)
    scheme_root.mkdir(parents=True, exist_ok=True)

    all_primers: list[tuple[str, str]] = []
    all_probes: list[tuple[str, str]] = []
    manifest_rows: list[list[object]] = []

    section("VarVAMP assay preparation")
    print("Mode         : QPCR")
    print(f"Schemes      : {len(schemes)}")
    if VERBOSE:
        print(f"Oligos FASTA : {oligos_file}")

    for scheme in sorted(schemes):
        oligos = schemes[scheme]
        missing = [kind for kind in ("LEFT", "PROBE", "RIGHT") if kind not in oligos]
        if missing:
            raise ValueError(
                f"Scheme {scheme} is incomplete. Missing: {', '.join(missing)}"
            )

        pair_records = [
            (f"{scheme}_LEFT", oligos["LEFT"]),
            (f"{scheme}_RIGHT", oligos["RIGHT"]),
        ]
        probe_records = [(f"{scheme}_PROBE", oligos["PROBE"])]

        write_fasta(primer_dir / f"{scheme}_primers.fasta", pair_records)
        write_fasta(probe_dir / f"{scheme}_probe.fasta", probe_records)
        all_primers.extend(pair_records)
        all_probes.extend(probe_records)

        one_dir = scheme_root / scheme
        write_fasta(one_dir / "primers.fasta", pair_records)
        write_fasta(one_dir / "probe.fasta", probe_records)
        with (one_dir / "scheme.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["scheme", "oligo_type", "sequence", "length"])
            for kind in ("LEFT", "PROBE", "RIGHT"):
                seq = oligos[kind]
                writer.writerow([scheme, kind, seq, len(seq)])

        manifest_rows.append(
            [
                scheme,
                oligos["LEFT"], oligos["PROBE"], oligos["RIGHT"],
                len(oligos["LEFT"]), len(oligos["PROBE"]), len(oligos["RIGHT"]),
            ]
        )
        print(
            f"  {scheme}: LEFT {len(oligos['LEFT'])} nt | "
            f"PROBE {len(oligos['PROBE'])} nt | RIGHT {len(oligos['RIGHT'])} nt"
        )

    write_fasta(primer_dir / "all_primers.fasta", all_primers)
    write_fasta(probe_dir / "all_probes.fasta", all_probes)
    with (assay_inputs_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme", "left_sequence", "probe_sequence", "right_sequence",
                "left_length", "probe_length", "right_length",
            ]
        )
        writer.writerows(manifest_rows)

    return schemes

# ---------------------------------------------------------------------------
# Target database selection and MFEprimer indexing
# ---------------------------------------------------------------------------

def is_inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def find_candidate_target_fastas(validation_root: Path | None = None) -> list[Path]:
    """Return user validation FASTA files from data/validation only."""
    root = (Path.cwd() / "data" / "validation").resolve()
    if not root.is_dir():
        return []
    candidates: list[Path] = []
    for pattern in ("*.fasta", "*.fa", "*.fna"):
        candidates.extend(path.resolve() for path in root.rglob(pattern) if path.is_file())
    return sorted(set(candidates))

def choose_target_database(
    requested_database: Path | None = None,
) -> Path:
    """Select the independent validation FASTA, normally from data/validation/."""
    print("\n2. Validation database")
    print("======================")

    if requested_database is not None:
        candidate = requested_database.expanduser().resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Validation FASTA not found: {candidate}")
        return candidate

    candidates = find_candidate_target_fastas()
    if candidates:
        for idx, path in enumerate(candidates, start=1):
            print(f"{idx}. {compact_path(path)}")
        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            raw = input("Select the validation database by number: ").strip()
            try:
                selected = int(raw)
            except ValueError:
                print("Please enter one of the displayed numbers.")
                continue
            if 1 <= selected <= len(candidates):
                return candidates[selected - 1]
            if selected == len(candidates) + 1:
                break
            print("Selection outside the available range.")

    while True:
        raw = input("Name or path of the validation FASTA: ").strip().strip("'\"")
        if not raw:
            print("A FASTA path is required.")
            continue
        candidate = Path(raw).expanduser()
        paths = [candidate]
        if not candidate.is_absolute():
            paths.extend(
                [
                    Path.cwd() / candidate,
                    Path.cwd() / "data" / "validation" / candidate,
                ]
            )
        for path in paths:
            if path.is_file():
                return path.resolve()
        print("FASTA file not found. Preferred location: data/validation/.")

def mfeprimer_index_files(database: Path) -> list[Path]:
    """Return recognized MFEprimer index files that currently exist."""
    candidates = [
        Path(str(database) + ".primerqc.bin"),
        Path(str(database) + ".primerqc"),
        Path(str(database) + ".primerqc.fai"),
    ]
    return [path for path in candidates if path.is_file()]


def mfeprimer_index_status(database: Path) -> tuple[bool, str]:
    """
    Recognize common MFEprimer index layouts.

    This function only recognizes index files. The downstream report parser
    still intentionally requires the legacy MFEprimer 3.x full text report.
    """
    binary_index = Path(str(database) + ".primerqc.bin")
    legacy_index = Path(str(database) + ".primerqc")
    legacy_fai = Path(str(database) + ".primerqc.fai")

    if binary_index.is_file():
        return True, "binary index (.primerqc.bin)"
    if legacy_index.is_file() and legacy_fai.is_file():
        return True, "legacy index (.primerqc + .primerqc.fai)"
    return False, "not indexed"

def ensure_validation_tools_available(assay_mode: str) -> None:
    """Check tools required for the selected validation branch."""
    tools = ["seqkit", "mfeprimer"]
    if assay_mode == "qpcr":
        tools.append("blastn")
    missing = [tool for tool in tools if shutil.which(tool) is None]

    if missing:
        raise RuntimeError(
            "Required validation tool(s) not found in PATH: "
            + ", ".join(missing)
            + ". Install/activate them before running coverage validation."
        )


def format_target_database_for_mfeprimer(
    database: Path,
    database_workdir: Path,
) -> Path:
    """
    Write a one-sequence-per-line FASTA and verify that complete headers and
    nucleotide content are preserved. Stale adjacent indexes are removed when
    the formatted FASTA changes.
    """
    database_workdir.mkdir(parents=True, exist_ok=True)
    safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", database.stem)
    fixed_database = database_workdir / f"{safe_stem}_fixed.fasta"
    temporary = database_workdir / f".{safe_stem}_fixed.tmp.fasta"

    print("\n4. MFEprimer database formatting")
    print("================================")
    progress("Formatting database with SeqKit", "RUN")
    run_command_to_file(["seqkit", "seq", "-w", "0", str(database)], temporary)

    source_records = read_fasta_with_full_headers(database)
    formatted_records = read_fasta_with_full_headers(temporary)
    headers_preserved = [h for h, _ in source_records] == [h for h, _ in formatted_records]
    sequences_preserved = [s for _, s in source_records] == [s for _, s in formatted_records]
    if not headers_preserved or not sequences_preserved:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            "SeqKit formatting changed FASTA headers or sequence content unexpectedly. "
            "The formatted database will not be used."
        )

    reused_existing = fixed_database.is_file() and fixed_database.read_bytes() == temporary.read_bytes()
    if reused_existing:
        temporary.unlink(missing_ok=True)
    else:
        # Prevent stale indexes from being associated with changed FASTA content.
        for index_file in mfeprimer_index_files(fixed_database):
            index_file.unlink(missing_ok=True)
        temporary.replace(fixed_database)

    source_stats = fasta_statistics(database)
    fixed_stats = fasta_statistics(fixed_database)
    progress("Formatting database with SeqKit", "OK")

    print("\nFormatting validation")
    print("---------------------")
    print(f"Input FASTA        : {compact_path(database)}")
    print(f"Formatted FASTA    : {compact_path(fixed_database)}")
    print(f"Sequences          : {source_stats['sequences']} -> {fixed_stats['sequences']}")
    print("Headers preserved  : YES")
    print("Sequences preserved: YES")
    print("IUPAC preserved    : YES")
    print("Modification       : FASTA line wrapping only")
    print(
        "Status             : "
        + (
            "EXISTING FORMATTED FASTA REUSED"
            if reused_existing
            else "FORMATTED FASTA CREATED/UPDATED"
        )
    )
    return fixed_database.resolve()

def index_target_database(database: Path) -> dict[str, object]:
    """
    Ensure that the formatted validation database has an MFEprimer index and
    explicitly report whether it was created or reused.
    """
    print("\n5. MFEprimer database indexing")
    print("==============================")

    indexed, detail = mfeprimer_index_status(database)
    existing_files = mfeprimer_index_files(database)

    if indexed:
        print("Existing MFEprimer index detected.")
        print("Status   : EXISTING INDEX REUSED")
        print(f"Database : {compact_path(database)}")
        print(f"Layout   : {detail}")
        print("\nIndex files detected")
        print("--------------------")
        for path in existing_files:
            print(f"✓ {compact_path(path)}")
        print("\nNo re-indexing was necessary.")
        return {
            "status": "REUSED",
            "detail": detail,
            "files": existing_files,
        }

    print("No complete MFEprimer index was detected for this formatted FASTA.")
    if not ask_yes_no("Index this database with MFEprimer now?"):
        raise RuntimeError(
            "Coverage validation requires an indexed MFEprimer database."
        )

    index_mode = ask_choice(
        "MFEprimer indexing parameters",
        [
            "Use MFEprimer default indexing parameters.",
            "Set a custom k-mer seed length (-k).",
        ],
    )

    command = ["mfeprimer", "index", "-i", str(database)]
    if index_mode == 2:
        k_value = ask_int("MFEprimer index k-mer length (-k)", 9, 15)
        command.extend(["-k", str(k_value)])

    progress("Indexing database with MFEprimer", "RUN")
    run_command(command)
    progress("Indexing database with MFEprimer", "OK")

    indexed, detail = mfeprimer_index_status(database)
    generated_files = mfeprimer_index_files(database)

    if not indexed:
        raise RuntimeError(
            "MFEprimer finished indexing, but the expected index layout "
            "could not be recognized. Check the installed MFEprimer version."
        )

    print("\nIndexing result")
    print("---------------")
    print("Status   : INDEX CREATED SUCCESSFULLY")
    print(f"Database : {compact_path(database)}")
    print(f"Layout   : {detail}")

    print("\nGenerated index files")
    print("---------------------")
    for path in generated_files:
        print(f"✓ {compact_path(path)}")

    print("\nMFEprimer database is ready for validation.")

    return {
        "status": "CREATED",
        "detail": detail,
        "files": generated_files,
    }


# ---------------------------------------------------------------------------
# Pair selection and MFEprimer execution
# ---------------------------------------------------------------------------

def choose_schemes(schemes: dict[str, dict[str, str]]) -> list[str]:
    """Select one, several, or all VarVAMP schemes."""
    scheme_names = sorted(schemes)

    print("\nAssays to validate")
    print("------------------")
    for idx, scheme in enumerate(scheme_names, start=1):
        probe_text = (
            f" | PROBE {len(schemes[scheme]['PROBE'])} nt"
            if "PROBE" in schemes[scheme]
            else ""
        )
        print(
            f"{idx}. {scheme}: LEFT {len(schemes[scheme]['LEFT'])} nt | "
            f"RIGHT {len(schemes[scheme]['RIGHT'])} nt{probe_text}"
        )
    print(f"{len(scheme_names) + 1}. ALL")

    while True:
        raw = input("Select one/multiple numbers (e.g. 1,3) or ALL: ").strip()
        try:
            numbers = [int(item.strip()) for item in raw.split(",") if item.strip()]
        except ValueError:
            print("Use numbers separated by commas, e.g. 1,2.")
            continue
        if not numbers:
            print("At least one selection is required.")
            continue

        all_option = len(scheme_names) + 1
        if all_option in numbers:
            if len(numbers) > 1:
                print("Choose ALL alone, or specific assay numbers.")
                continue
            return scheme_names

        if any(number < 1 or number > len(scheme_names) for number in numbers):
            print("One or more selected numbers are outside the available range.")
            continue

        selected: list[str] = []
        for number in numbers:
            scheme = scheme_names[number - 1]
            if scheme not in selected:
                selected.append(scheme)
        return selected


def choose_mfeprimer_parameters() -> dict[str, float | int] | None:
    print("\n6. Validation parameters")
    print("========================")
    mode = ask_choice(
        "MFEprimer search filters",
        [
            "Use defaults (0-2000 bp; Tm >=30 °C).",
            "Set custom amplicon size and Tm thresholds.",
        ],
    )
    if mode == 1:
        return None

    min_size = ask_int("Minimum predicted amplicon size (-s)", 0, 100000)
    max_size = ask_int("Maximum predicted amplicon size (-S)", min_size, 100000)
    tm_cutoff = ask_float("Minimum Tm cutoff (-t)", 0.0, 100.0)
    return {
        "min_size": min_size,
        "max_size": max_size,
        "tm_cutoff": tm_cutoff,
    }



# ---------------------------------------------------------------------------
# MFEprimer 3.x report parsing and target-amplicon extraction
# ---------------------------------------------------------------------------

def choose_coverage_amplicon_filter() -> tuple[int | None, int | None]:
    """
    Choose the size range used to count an amplicon as a target-coverage hit.

    This filter is applied AFTER MFEprimer has generated its report, so
    non-target/self-priming products remain visible in the QC output.
    """
    mode = ask_choice(
        "Target-amplicon extraction filter",
        [
            "qPCR range: keep target amplicons from 70 to 250 bp.",
            "Enter a custom target-amplicon size range.",
            "Do not apply an additional size filter.",
        ],
    )

    if mode == 1:
        return 70, 250

    if mode == 2:
        min_size = ask_int(
            "Minimum target amplicon size",
            0,
            100000,
        )
        max_size = ask_int(
            "Maximum target amplicon size",
            min_size,
            100000,
        )
        return min_size, max_size

    return None, None


def normalize_primer_identifier(value: str) -> str:
    """
    Normalize an identifier for robust comparisons.

    Examples:
        varVAMP_0      -> varvamp0
        VARVAMP-0      -> varvamp0
        varVAMP_0_LEFT -> varvamp0left
    """
    return re.sub(
        r"[^A-Za-z0-9]+",
        "",
        value,
    ).lower()


def parse_primer_identity(
    primer_name: str,
) -> tuple[str | None, str | None]:
    """
    Extract the VarVAMP scheme and LEFT/RIGHT role from a MFEprimer
    primer identifier.

    This intentionally tolerates harmless decorations added around the
    FASTA identifier, for example:
        varVAMP_0_LEFT
        varVAMP_0_LEFT(+)
        varVAMP_0-LEFT
        varVAMP_0_LEFT_1

    LEFT/RIGHT must still appear as a real token, not as part of another
    ordinary word.
    """
    raw = primer_name.strip().lstrip(">")

    role_match = re.search(
        r"(^|[^A-Za-z0-9])(LEFT|RIGHT)(?=$|[^A-Za-z0-9])",
        raw,
        flags=re.IGNORECASE,
    )

    if not role_match:
        return None, None

    role = role_match.group(2).upper()

    # Everything before the LEFT/RIGHT token is treated as the scheme
    # component. Remove punctuation/separators surrounding the token.
    scheme = raw[:role_match.start(2)]
    scheme = re.sub(
        r"[^A-Za-z0-9]+$",
        "",
        scheme,
    ).strip()

    if not scheme:
        return None, role

    return scheme, role


def same_scheme(
    observed: str | None,
    expected: str,
) -> bool:
    """Compare scheme identifiers while ignoring separators and case."""
    if observed is None:
        return False

    return (
        normalize_primer_identifier(observed)
        == normalize_primer_identifier(expected)
    )


def classify_mfeprimer_product(
    primer_1: str,
    primer_2: str,
    expected_scheme: str,
) -> tuple[str, str, str]:
    """
    Classify a predicted MFEprimer product.

    TARGET:
        one LEFT and one RIGHT primer from the selected VarVAMP scheme.

    SELF_PRIMING:
        LEFT+LEFT or RIGHT+RIGHT from the selected scheme.

    CROSS_SCHEME:
        LEFT/RIGHT roles are recognized, but the primers do not belong to
        the same scheme.

    OTHER:
        at least one raw primer identifier cannot be interpreted safely.
    """
    scheme_1, role_1 = parse_primer_identity(primer_1)
    scheme_2, role_2 = parse_primer_identity(primer_2)

    if role_1 is None or role_2 is None:
        return (
            "OTHER",
            role_1 or "UNKNOWN",
            role_2 or "UNKNOWN",
        )

    # If roles were recognized but either scheme prefix was not recoverable,
    # do not guess.
    if scheme_1 is None or scheme_2 is None:
        return "OTHER", role_1, role_2

    if normalize_primer_identifier(
        scheme_1
    ) != normalize_primer_identifier(
        scheme_2
    ):
        return "CROSS_SCHEME", role_1, role_2

    if not same_scheme(
        scheme_1,
        expected_scheme,
    ):
        return "OTHER", role_1, role_2

    roles = {role_1, role_2}

    if roles == {"LEFT", "RIGHT"}:
        return "TARGET", role_1, role_2

    if (
        role_1 == role_2
        and role_1 in {"LEFT", "RIGHT"}
    ):
        return "SELF_PRIMING", role_1, role_2

    return "OTHER", role_1, role_2



def parse_mfeprimer3_report(report_file: Path) -> dict[str, object]:
    """
    Parse a legacy MFEprimer 3.x specificity text report.

    The parser deliberately avoids fixed-width column positions.

    It uses two independent structures:
    1. the AmpID/HitID/Size summary table;
    2. FASTA-like blocks beginning with:
           >Amp_N PRIMER1 + PRIMER2 ==> HitID

    The sequence block is considered authoritative for extraction, while
    the summary table is used for consistency checks.
    """
    if not report_file.is_file():
        raise FileNotFoundError(
            f"MFEprimer report not found: {report_file}"
        )

    lines = report_file.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    expected_count: int | None = None

    for line in lines:
        match = re.search(
            r"Descriptions\s+of\s+\[\s*(\d+)\s*\]\s+potential\s+amplicons",
            line,
            flags=re.IGNORECASE,
        )
        if match:
            expected_count = int(match.group(1))
            break

    # Parse the summary table without relying on fixed-width spacing.
    summary_by_ampid: dict[int, dict[str, object]] = {}
    in_table = False

    for line in lines:
        stripped = line.strip()

        if re.match(
            r"^AmpID\s+HitID\s+Size\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            in_table = True
            continue

        if in_table and re.match(
            r"^Amplicon\s+details\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            break

        if not in_table:
            continue

        row = re.match(
            r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+",
            line,
        )
        if not row:
            continue

        amp_id = int(row.group(1))
        hit_id = row.group(2)
        size = int(row.group(3))

        if amp_id in summary_by_ampid:
            raise RuntimeError(
                f"Duplicate AmpID {amp_id} in MFEprimer summary table."
            )

        summary_by_ampid[amp_id] = {
            "hit_id": hit_id,
            "size": size,
        }

    # Parse FASTA-like amplicon blocks.
    header_re = re.compile(
        r"^>Amp_(\d+)\s+(\S+)\s+\+\s+(\S+)\s+==>\s+(\S+)"
    )
    sequence_re = re.compile(
        r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$"
    )

    parsed_blocks: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    sequence_parts: list[str] = []

    def flush_current() -> None:
        nonlocal current, sequence_parts

        if current is None:
            return

        sequence = "".join(sequence_parts).upper()

        if not sequence:
            raise RuntimeError(
                f"No amplicon sequence found for Amp_{current['amp_id']}."
            )

        current["sequence"] = sequence
        current["sequence_length"] = len(sequence)
        parsed_blocks.append(current)

        current = None
        sequence_parts = []

    for raw_line in lines:
        stripped = raw_line.strip()

        header_match = header_re.match(stripped)

        if header_match:
            flush_current()

            current = {
                "amp_id": int(header_match.group(1)),
                "primer_1": header_match.group(2),
                "primer_2": header_match.group(3),
                "hit_id": header_match.group(4),
            }
            sequence_parts = []
            continue

        if current is None:
            continue

        if stripped and sequence_re.fullmatch(stripped):
            sequence_parts.append(stripped)
            continue

        # A blank line or a non-sequence line terminates the FASTA block.
        if sequence_parts:
            flush_current()

    flush_current()

    # If no amplicons were reported, allow a clean empty result.
    if expected_count == 0 and not parsed_blocks:
        return {
            "expected_count": 0,
            "summary_by_ampid": summary_by_ampid,
            "amplicons": [],
        }

    if not parsed_blocks:
        raise RuntimeError(
            "No '>Amp_N ... ==> HitID' amplicon sequence blocks were found "
            "in the MFEprimer report. This parser expects the legacy "
            "MFEprimer 3.x text report."
        )

    amp_ids = [int(item["amp_id"]) for item in parsed_blocks]

    if len(amp_ids) != len(set(amp_ids)):
        raise RuntimeError(
            "Duplicate AmpID values were found in the MFEprimer sequence blocks."
        )

    if (
        expected_count is not None
        and len(parsed_blocks) != expected_count
    ):
        raise RuntimeError(
            "MFEprimer report consistency check failed: "
            f"report announces {expected_count} potential amplicons, "
            f"but {len(parsed_blocks)} sequence blocks were parsed."
        )

    # Cross-check HitID and Size between the table and sequence blocks.
    for item in parsed_blocks:
        amp_id = int(item["amp_id"])
        table_item = summary_by_ampid.get(amp_id)

        if table_item is None:
            raise RuntimeError(
                f"Amp_{amp_id} is present in amplicon details but missing "
                "from the AmpID/HitID/Size table."
            )

        if str(table_item["hit_id"]) != str(item["hit_id"]):
            raise RuntimeError(
                f"HitID mismatch for Amp_{amp_id}: "
                f"table={table_item['hit_id']} vs "
                f"details={item['hit_id']}."
            )

        reported_size = int(table_item["size"])
        sequence_length = int(item["sequence_length"])

        if reported_size != sequence_length:
            raise RuntimeError(
                f"Size mismatch for Amp_{amp_id}: "
                f"table reports {reported_size} bp but extracted sequence "
                f"contains {sequence_length} nt."
            )

        item["size"] = reported_size

    return {
        "expected_count": expected_count,
        "summary_by_ampid": summary_by_ampid,
        "amplicons": parsed_blocks,
    }


def analyse_and_extract_mfeprimer_amplicons(
    report_file: Path,
    scheme: str,
    scheme_dir: Path,
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
    exclude_self_priming: bool = True,
) -> dict[str, object]:
    """
    Classify MFEprimer products, extract valid target amplicons and calculate
    unique primer-pair target coverage.
    """
    parsed = parse_mfeprimer3_report(report_file)
    amplicons = list(parsed["amplicons"])

    all_rows: list[dict[str, object]] = []
    valid_rows: list[dict[str, object]] = []
    non_target_rows: list[dict[str, object]] = []

    for item in amplicons:
        classification, role_1, role_2 = classify_mfeprimer_product(
            str(item["primer_1"]),
            str(item["primer_2"]),
            scheme,
        )

        size = int(item["size"])

        if min_target_size is None or max_target_size is None:
            size_pass = True
        else:
            size_pass = min_target_size <= size <= max_target_size

        allowed_classes = (
            {"TARGET"}
            if exclude_self_priming
            else {"TARGET", "SELF_PRIMING"}
        )
        coverage_eligible = (
            classification in allowed_classes
            and size_pass
        )

        row = {
            "AmpID": int(item["amp_id"]),
            "HitID": str(item["hit_id"]),
            "Primer1": str(item["primer_1"]),
            "Primer2": str(item["primer_2"]),
            "Role1": role_1,
            "Role2": role_2,
            "ProductClass": classification,
            "Size": size,
            "SizePass": "YES" if size_pass else "NO",
            "CoverageEligible": "YES" if coverage_eligible else "NO",
            "Sequence": str(item["sequence"]).upper(),
        }

        all_rows.append(row)

        if coverage_eligible:
            valid_rows.append(row)

        if classification != "TARGET":
            non_target_rows.append(row)

    # Stable first-occurrence order, without double-counting a HitID.
    unique_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    valid_amplicons = len(valid_rows)
    unique_positive_targets = len(unique_positive_hitids)
    duplicate_target_products = (
        valid_amplicons - unique_positive_targets
    )

    if total_target_sequences <= 0:
        raise ValueError(
            "Total target sequence count must be greater than zero."
        )

    coverage_pct = (
        100.0
        * unique_positive_targets
        / total_target_sequences
    )

    # Output files.
    all_tsv = scheme_dir / "amplicons_all.tsv"
    valid_tsv = scheme_dir / "amplicons_valid.tsv"
    valid_fasta = scheme_dir / "amplicons_valid.fasta"
    hitids_file = scheme_dir / "positive_hitids.txt"
    non_target_tsv = scheme_dir / "non_target_products.tsv"
    summary_file = scheme_dir / "coverage_summary.txt"

    tsv_columns = [
        "AmpID",
        "HitID",
        "Primer1",
        "Primer2",
        "Role1",
        "Role2",
        "ProductClass",
        "Size",
        "SizePass",
        "CoverageEligible",
    ]

    def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=tsv_columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_rows(all_tsv, all_rows)
    write_rows(valid_tsv, valid_rows)
    write_rows(non_target_tsv, non_target_rows)

    with valid_fasta.open("w", encoding="utf-8") as handle:
        for row in valid_rows:
            pair = f"{row['Role1']}+{row['Role2']}"
            handle.write(
                f">Amp_{row['AmpID']}|HitID={row['HitID']}|"
                f"Scheme={scheme}|Pair={pair}|Size={row['Size']}\n"
            )
            sequence = str(row["Sequence"])
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    hitids_file.write_text(
        "".join(f"{hit_id}\n" for hit_id in unique_positive_hitids),
        encoding="utf-8",
    )

    class_counts: dict[str, int] = {}
    for row in all_rows:
        key = str(row["ProductClass"])
        class_counts[key] = class_counts.get(key, 0) + 1

    unrecognized_examples: list[str] = []
    seen_unrecognized: set[tuple[str, str, str]] = set()

    for row in all_rows:
        product_class = str(row["ProductClass"])

        if product_class not in {"OTHER", "CROSS_SCHEME"}:
            continue

        key = (
            str(row["Primer1"]),
            str(row["Primer2"]),
            product_class,
        )

        if key in seen_unrecognized:
            continue

        seen_unrecognized.add(key)
        unrecognized_examples.append(
            f"{product_class}: "
            f"{row['Primer1']} + {row['Primer2']}"
        )

        if len(unrecognized_examples) >= 10:
            break

    if min_target_size is None or max_target_size is None:
        size_filter_text = "none"
    else:
        size_filter_text = f"{min_target_size}-{max_target_size} bp"

    summary_lines = [
        f"MFEprimer target-coverage summary — {scheme}",
        "=" * (38 + len(scheme)),
        "",
        f"MFEprimer report: {compact_path(report_file)}",
        f"Total validation sequences: {total_target_sequences}",
        f"Potential amplicons parsed: {len(all_rows)}",
        f"TARGET products (LEFT+RIGHT or RIGHT+LEFT): "
        f"{class_counts.get('TARGET', 0)}",
        f"SELF_PRIMING products (LEFT+LEFT or RIGHT+RIGHT): "
        f"{class_counts.get('SELF_PRIMING', 0)}",
        f"CROSS_SCHEME products: {class_counts.get('CROSS_SCHEME', 0)}",
        f"OTHER products: {class_counts.get('OTHER', 0)}",
        f"Target amplicon size filter: {size_filter_text}",
        f"Retained amplicons after filtering: {valid_amplicons}",
        f"Unique positive HitIDs: {unique_positive_targets}",
        f"Additional valid amplicons on already-positive HitIDs: "
        f"{duplicate_target_products}",
        (
            "Primer-pair target coverage: "
            f"{unique_positive_targets} / {total_target_sequences} "
            f"= {coverage_pct:.2f} %"
        ),
        "",
        "Coverage-set definition:",
        (
            "- TARGET only (LEFT+RIGHT / RIGHT+LEFT);"
            if exclude_self_priming
            else "- TARGET plus SELF_PRIMING products;"
        ),
        "- optional target-amplicon size filter passed.",
    ]

    if unrecognized_examples:
        summary_lines.extend(
            [
                "",
                "Examples of products not recognized as TARGET/SELF_PRIMING:",
                *[
                    f"- {example}"
                    for example in unrecognized_examples
                ],
            ]
        )

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"MFEprimer results — {scheme}")
    print(f"Validation sequences       {total_target_sequences}")
    print(f"Potential products         {len(all_rows)}")
    print(f"Retained coverage products {valid_amplicons}")
    print(f"Self-priming products      {class_counts.get('SELF_PRIMING', 0)}")
    print(f"Unique PCR-positive HitIDs {unique_positive_targets}")
    print(f"PCR coverage               {coverage_pct:.2f} %")
    if class_counts.get("CROSS_SCHEME", 0):
        print(f"Cross-scheme products      {class_counts.get('CROSS_SCHEME', 0)}")
    if class_counts.get("OTHER", 0):
        print(f"Unclassified products      {class_counts.get('OTHER', 0)}")
    if unrecognized_examples and VERBOSE:
        print("\nUnrecognized primer-name examples:")
        for example in unrecognized_examples[:5]:
            print(f"- {example}")
    if VERBOSE:
        print(f"Valid amplicons TSV        {valid_tsv}")
        print(f"Valid amplicons FASTA      {valid_fasta}")
        print(f"Positive HitIDs            {hitids_file}")
        print(f"Non-target products        {non_target_tsv}")
        print(f"Coverage summary           {summary_file}")

    return {
        "scheme": scheme,
        "potential_amplicons": len(all_rows),
        "target_products": class_counts.get("TARGET", 0),
        "self_priming_products": class_counts.get("SELF_PRIMING", 0),
        "valid_target_amplicons": valid_amplicons,
        "unique_positive_hitids": unique_positive_targets,
        "coverage_pct": coverage_pct,
        "valid_tsv": valid_tsv,
        "valid_fasta": valid_fasta,
        "positive_hitids": hitids_file,
        "summary": summary_file,
        "valid_rows": valid_rows,
        "self_priming_excluded": exclude_self_priming,
    }



# ---------------------------------------------------------------------------
# BLAST+ analysis of the original VarVAMP probe
# ---------------------------------------------------------------------------

PROBE_IUPAC_TO_BASES = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}


def expand_iupac_probe(
    sequence: str,
    max_variants: int = 256,
) -> list[str]:
    """
    Expand an already-degenerate VarVAMP probe into concrete A/C/G/T queries.

    This is NOT probe optimization. It only allows BLAST to evaluate an
    existing IUPAC probe because short-query BLAST works best with concrete
    nucleotide queries.
    """
    sequence = sequence.upper()

    invalid = sorted(set(sequence) - set(PROBE_IUPAC_TO_BASES))
    if invalid:
        raise ValueError(
            "Unsupported probe character(s): " + ", ".join(invalid)
        )

    degeneracy = 1
    for symbol in sequence:
        degeneracy *= len(PROBE_IUPAC_TO_BASES[symbol])

    if degeneracy > max_variants:
        raise RuntimeError(
            f"Probe degeneracy is {degeneracy}, which would require more "
            f"than {max_variants} concrete BLAST queries. "
            "This BLAST-only version stops instead of silently truncating "
            "the probe variants."
        )

    variants = [""]

    for symbol in sequence:
        variants = [
            prefix + base
            for prefix in variants
            for base in PROBE_IUPAC_TO_BASES[symbol]
        ]

    return variants


def write_probe_blast_queries(
    scheme: str,
    probe: str,
    output_file: Path,
) -> list[str]:
    """Write one or more concrete probe variants for blastn-short."""
    variants = expand_iupac_probe(probe)

    with output_file.open("w", encoding="utf-8") as handle:
        for idx, variant in enumerate(variants, start=1):
            handle.write(
                f">{scheme}_PROBE_variant_{idx}\n"
                f"{variant}\n"
            )

    return variants


def parse_amplicon_subject_id(subject_id: str) -> int | None:
    """
    Extract AmpID from headers written by amplicons_valid.fasta.

    Example:
        Amp_2|HitID=OR613708.1|Scheme=varVAMP_0|...
        -> 2
    """
    match = re.match(r"^Amp_(\d+)(?:\||$)", subject_id)

    if not match:
        return None

    return int(match.group(1))


def run_probe_blast(
    scheme: str,
    probe: str,
    valid_amplicons_fasta: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
) -> tuple[Path, Path, int]:
    """
    Run blastn-short directly against amplicons_valid.fasta using -subject.

    No second BLAST database is required. The raw tabular BLAST output is
    preserved for later manual/scientific analysis.
    """
    query_file = scheme_dir / "probe_blast_queries.fasta"
    variants = write_probe_blast_queries(
        scheme,
        probe,
        query_file,
    )

    raw_file = scheme_dir / "probe_blast_raw.tsv"

    # Keep enough target sequences even when a future dataset contains
    # substantially more than BLAST's default number of reported targets.
    max_targets = max(1000, len(valid_rows) * 2)

    outfmt = (
        "6 "
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qseq sseq"
    )

    command = [
        "blastn",
        "-task",
        "blastn-short",
        "-query",
        str(query_file),
        "-subject",
        str(valid_amplicons_fasta),
        "-strand",
        "both",
        "-dust",
        "no",
        "-evalue",
        "1000",
        "-max_target_seqs",
        str(max_targets),
        "-max_hsps",
        "20",
        "-outfmt",
        outfmt,
        "-out",
        str(raw_file),
    ]

    ambiguous_positions = [
        f"{idx}:{base}"
        for idx, base in enumerate(probe.upper(), start=1)
        if base not in {"A", "C", "G", "T"}
    ]

    section(f"Probe BLAST — {scheme}")
    print(f"Original probe    : {probe}")
    print(f"Length            : {len(probe)} nt")
    print(
        "IUPAC ambiguity   : "
        + (", ".join(ambiguous_positions) if ambiguous_positions else "none")
    )
    print(f"Concrete variants : {len(variants)}")
    print(f"Target amplicons  : {len(valid_rows)}")
    if VERBOSE:
        print("BLAST task        : blastn-short")
        print("Strands           : both")
        print("DUST              : disabled")
        print("E-value           : 1000")

    run_command(command)

    return raw_file, query_file, len(variants)


def read_probe_blast_raw(
    raw_file: Path,
) -> list[dict[str, object]]:
    """Parse the custom BLAST outfmt 6 produced by run_probe_blast()."""
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    hits: list[dict[str, object]] = []

    if not raw_file.is_file():
        raise FileNotFoundError(
            f"BLAST output file not found: {raw_file}"
        )

    with raw_file.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != len(columns):
                raise RuntimeError(
                    f"Unexpected BLAST field count on line {line_number}: "
                    f"expected {len(columns)}, found {len(fields)}."
                )

            row = dict(zip(columns, fields))

            amp_id = parse_amplicon_subject_id(
                str(row["sseqid"])
            )

            row.update(
                {
                    "AmpID": amp_id,
                    "pident": float(row["pident"]),
                    "length": int(row["length"]),
                    "mismatch": int(row["mismatch"]),
                    "gapopen": int(row["gapopen"]),
                    "qstart": int(row["qstart"]),
                    "qend": int(row["qend"]),
                    "sstart": int(row["sstart"]),
                    "send": int(row["send"]),
                    "evalue": float(row["evalue"]),
                    "bitscore": float(row["bitscore"]),
                }
            )

            hits.append(row)

    return hits


def classify_full_length_probe_hit(
    hit: dict[str, object],
    probe_length: int,
) -> tuple[bool, str]:
    """
    Define the transparent BLAST classes used by this workflow.

    Full-length ungapped means:
    - alignment starts at query base 1;
    - alignment ends at the complete probe length;
    - aligned length equals probe length;
    - zero gap openings.

    No biological PASS/FAIL label is imposed.
    """
    full_length = (
        int(hit["qstart"]) == 1
        and int(hit["qend"]) == probe_length
        and int(hit["length"]) == probe_length
        and int(hit["gapopen"]) == 0
    )

    if not full_length:
        return False, "PARTIAL_OR_GAPPED"

    mismatches = int(hit["mismatch"])

    if mismatches == 0:
        return True, "0_MISMATCH"

    if mismatches == 1:
        return True, "1_MISMATCH"

    if mismatches == 2:
        return True, "2_MISMATCHES"

    return True, "GE3_MISMATCHES"



def build_blast_match_line(
    query_aligned: str,
    subject_aligned: str,
) -> str:
    """
    Create a BLAST-like visual match line.

    |  exact nucleotide identity
    .  aligned non-identical bases
       gap
    """
    symbols: list[str] = []

    for query_base, subject_base in zip(
        query_aligned.upper(),
        subject_aligned.upper(),
    ):
        if query_base == "-" or subject_base == "-":
            symbols.append(" ")
        elif query_base == subject_base:
            symbols.append("|")
        else:
            symbols.append(".")

    return "".join(symbols)


def extract_full_length_mismatches(
    row: dict[str, object],
) -> list[dict[str, object]]:
    """
    Extract mismatch positions from a full-length ungapped BLAST alignment.

    ProbePosition is 1-based relative to the concrete BLAST query variant.
    SubjectCoordinate respects the BLAST subject orientation.
    """
    if str(row["FullLengthUngapped"]) != "YES":
        return []

    query = str(row["QueryAligned"]).upper()
    subject = str(row["SubjectAligned"]).upper()

    if len(query) != len(subject):
        return []

    subject_start = int(row["SubjectStart"])
    subject_end = int(row["SubjectEnd"])
    subject_step = 1 if subject_end >= subject_start else -1
    subject_coord = subject_start

    mismatch_rows: list[dict[str, object]] = []

    for probe_position, (query_base, subject_base) in enumerate(
        zip(query, subject),
        start=1,
    ):
        if query_base == "-" or subject_base == "-":
            if subject_base != "-":
                subject_coord += subject_step
            continue

        current_subject_coord = subject_coord

        if query_base != subject_base:
            mismatch_rows.append(
                {
                    "AmpID": int(row["AmpID"]),
                    "HitID": str(row["HitID"]),
                    "BlastCategory": str(row["BlastCategory"]),
                    "QueryVariant": str(row["QueryVariant"]),
                    "ProbePosition": probe_position,
                    "ProbeBase": query_base,
                    "TargetBase": subject_base,
                    "SubjectCoordinate": current_subject_coord,
                    "SubjectStrand": (
                        "+"
                        if subject_step == 1
                        else "-"
                    ),
                }
            )

        subject_coord += subject_step

    return mismatch_rows


def format_probe_blast_alignment_block(
    row: dict[str, object],
) -> str:
    """Create one human-readable BLAST alignment block."""
    category = str(row["BlastCategory"])
    hit_id = str(row["HitID"])
    amp_id = int(row["AmpID"])
    amp_size = int(row["AmpliconSize"])

    lines = [
        "=" * 78,
        f"HitID: {hit_id} | AmpID: {amp_id} | Amplicon size: {amp_size} bp",
        f"Category: {category}",
    ]

    if category == "NO_HIT":
        lines.extend(
            [
                "No BLAST HSP was reported for this amplicon.",
                "",
            ]
        )
        return "\n".join(lines)

    qseq = str(row["QueryAligned"])
    sseq = str(row["SubjectAligned"])
    match_line = build_blast_match_line(qseq, sseq)

    sstart = int(row["SubjectStart"])
    send = int(row["SubjectEnd"])
    strand = "+" if send >= sstart else "-"

    lines.extend(
        [
            (
                f"Identity: {row['PercentIdentity']}% | "
                f"Alignment: {row['AlignmentLength']} nt | "
                f"Mismatches: {row['Mismatches']} | "
                f"Gap openings: {row['GapOpenings']}"
            ),
            (
                f"Query variant: {row['QueryVariant']} | "
                f"Query coordinates: {row['QueryStart']}..{row['QueryEnd']}"
            ),
            (
                f"Subject coordinates: {sstart}..{send} | "
                f"Subject strand: {strand}"
            ),
            f"E-value: {row['Evalue']} | Bit score: {row['BitScore']}",
            "",
            f"Probe  {str(row['QueryStart']).rjust(4)}  {qseq}  {row['QueryEnd']}",
            f"             {match_line}",
            f"Target {str(sstart).rjust(4)}  {sseq}  {send}",
        ]
    )

    mismatch_rows = extract_full_length_mismatches(row)

    if str(row["FullLengthUngapped"]) == "YES":
        if mismatch_rows:
            mismatch_text = ", ".join(
                (
                    f"pos{item['ProbePosition']} "
                    f"{item['ProbeBase']}>{item['TargetBase']} "
                    f"(target:{item['SubjectCoordinate']})"
                )
                for item in mismatch_rows
            )
            lines.extend(
                [
                    "",
                    f"Mismatch positions: {mismatch_text}",
                ]
            )
        else:
            lines.extend(
                [
                    "",
                    "Mismatch positions: none",
                ]
            )
    else:
        lines.extend(
            [
                "",
                (
                    "Note: this is not a full-length ungapped alignment; "
                    "mismatch-position interpretation is not used."
                ),
            ]
        )

    lines.append("")
    return "\n".join(lines)


def write_probe_blast_alignment_reports(
    amplicon_summary: list[dict[str, object]],
    hitid_summary: list[dict[str, object]],
    scheme_dir: Path,
) -> tuple[Path, Path, Path]:
    """
    Write human-readable alignments and a machine-readable mismatch table.
    """
    category_order = {
        "0_MISMATCH": 0,
        "1_MISMATCH": 1,
        "2_MISMATCHES": 2,
        "GE3_MISMATCHES": 3,
        "PARTIAL_ONLY": 4,
        "NO_HIT": 5,
    }

    def sort_rows(
        rows: list[dict[str, object]],
    ) -> list[dict[str, object]]:
        return sorted(
            rows,
            key=lambda row: (
                category_order.get(
                    str(row["BlastCategory"]),
                    99,
                ),
                str(row["HitID"]),
                int(row["AmpID"]),
            ),
        )

    amplicon_report = (
        scheme_dir / "probe_blast_alignments_by_amplicon.txt"
    )
    hitid_report = (
        scheme_dir / "probe_blast_alignments_by_hitid.txt"
    )
    mismatch_file = (
        scheme_dir / "probe_blast_mismatch_positions.tsv"
    )

    with amplicon_report.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "BLAST probe alignments — one block per valid MFEprimer amplicon\n"
        )
        handle.write(
            "Symbols: | exact identity, . mismatch, blank = gap\n\n"
        )

        for row in sort_rows(amplicon_summary):
            handle.write(
                format_probe_blast_alignment_block(row)
            )

    with hitid_report.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "BLAST probe alignments — best amplicon per unique HitID\n"
        )
        handle.write(
            "Symbols: | exact identity, . mismatch, blank = gap\n\n"
        )

        for row in sort_rows(hitid_summary):
            handle.write(
                format_probe_blast_alignment_block(row)
            )

    mismatch_rows: list[dict[str, object]] = []

    for row in hitid_summary:
        mismatch_rows.extend(
            extract_full_length_mismatches(row)
        )

    mismatch_columns = [
        "AmpID",
        "HitID",
        "BlastCategory",
        "QueryVariant",
        "ProbePosition",
        "ProbeBase",
        "TargetBase",
        "SubjectCoordinate",
        "SubjectStrand",
    ]

    with mismatch_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=mismatch_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(mismatch_rows)

    return (
        amplicon_report,
        hitid_report,
        mismatch_file,
    )




def summarize_probe_variants(
    query_file: Path,
    raw_hits: list[dict[str, object]],
    valid_rows: list[dict[str, object]],
    probe_length: int,
    scheme_dir: Path,
    category_rank: dict[str, int],
) -> dict[str, object]:
    """
    Analyse every concrete IUPAC-expanded probe variant independently.

    Statistics are calculated at unique HitID level:
      - each concrete probe variant is evaluated against every PCR-positive HitID;
      - if one HitID has multiple valid MFEprimer amplicons, the best result for
        that concrete variant is retained once;
      - variants are NOT merged in this table.

    A separate best-of-all-variants summary is still produced elsewhere and is
    used for the final qPCR coverage.
    """
    query_records = read_fasta(query_file)

    variants = [
        {
            "query_id": header,
            "sequence": sequence,
        }
        for header, sequence in query_records
    ]

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    # Stable first-occurrence order of unique biological PCR-positive targets.
    pcr_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    # Hits have already been classified as full-length / partial by the caller.
    hits_by_variant_amp: dict[
        tuple[str, int],
        list[dict[str, object]],
    ] = {}

    for hit in raw_hits:
        amp_id = hit.get("AmpID")
        if amp_id is None:
            continue

        amp_id = int(amp_id)
        if amp_id not in valid_by_ampid:
            continue

        query_id = str(hit["qseqid"])
        hits_by_variant_amp.setdefault(
            (query_id, amp_id),
            [],
        ).append(hit)

    def choose_best_amp_result(
        query_id: str,
        amp_id: int,
    ) -> dict[str, object]:
        amp_hits = hits_by_variant_amp.get(
            (query_id, amp_id),
            [],
        )

        full_hits = [
            hit
            for hit in amp_hits
            if hit.get("FullLengthUngapped") == "YES"
        ]

        if full_hits:
            full_hits.sort(
                key=lambda hit: (
                    int(hit["mismatch"]),
                    -float(hit["pident"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = full_hits[0]
            return {
                "category": str(best["BlastCategory"]),
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        if amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            return {
                "category": "PARTIAL_ONLY",
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        return {
            "category": "NO_HIT",
            "mismatches": 999999,
            "alignment_length": 0,
            "bitscore": 0.0,
            "amp_id": amp_id,
        }

    def result_rank(result: dict[str, object]) -> tuple:
        return (
            category_rank[str(result["category"])],
            int(result["mismatches"]),
            -int(result["alignment_length"]),
            -float(result["bitscore"]),
            int(result["amp_id"]),
        )

    # Map each HitID to all its valid amplicons.
    ampids_by_hitid: dict[str, list[int]] = {}
    for row in valid_rows:
        ampids_by_hitid.setdefault(
            str(row["HitID"]),
            [],
        ).append(int(row["AmpID"]))

    matrix: dict[str, dict[str, str]] = {
        hit_id: {}
        for hit_id in pcr_positive_hitids
    }

    stats_rows: list[dict[str, object]] = []

    for variant in variants:
        query_id = str(variant["query_id"])
        sequence = str(variant["sequence"])

        # A. Statistics by unique PCR-positive HitID.
        #    These remain the basis of biological coverage so one biological
        #    target is not counted several times when MFEprimer reports more
        #    than one retained amplicon for the same HitID.
        counts = {
            category: 0
            for category in category_rank
        }

        for hit_id in pcr_positive_hitids:
            amp_results = [
                choose_best_amp_result(
                    query_id,
                    amp_id,
                )
                for amp_id in ampids_by_hitid[hit_id]
            ]

            best = min(
                amp_results,
                key=result_rank,
            )
            category = str(best["category"])

            counts[category] += 1
            matrix[hit_id][query_id] = category

        full_length = (
            counts["0_MISMATCH"]
            + counts["1_MISMATCH"]
            + counts["2_MISMATCHES"]
            + counts["GE3_MISMATCHES"]
        )

        denominator = len(pcr_positive_hitids)

        # B. Statistics by retained TARGET amplicon.
        #    This is displayed immediately after probe BLAST so the user can
        #    see, for every concrete probe version, how many of the retained
        #    amplicons are exact, 1-mismatch, 2-mismatch, etc.
        amplicon_counts = {
            category: 0
            for category in category_rank
        }

        for row in valid_rows:
            amp_id = int(row["AmpID"])
            best_amp = choose_best_amp_result(query_id, amp_id)
            amplicon_counts[str(best_amp["category"])] += 1

        amplicon_full_length = (
            amplicon_counts["0_MISMATCH"]
            + amplicon_counts["1_MISMATCH"]
            + amplicon_counts["2_MISMATCHES"]
            + amplicon_counts["GE3_MISMATCHES"]
        )

        stats_rows.append(
            {
                "QueryVariant": query_id,
                "ConcreteSequence": sequence,

                # Retained TARGET amplicon-level counts
                "TargetAmplicons": len(valid_rows),
                "AmpliconExact0Mismatch": amplicon_counts["0_MISMATCH"],
                "AmpliconOneMismatch": amplicon_counts["1_MISMATCH"],
                "AmpliconTwoMismatches": amplicon_counts["2_MISMATCHES"],
                "AmpliconGE3Mismatches": amplicon_counts["GE3_MISMATCHES"],
                "AmpliconPartialOnly": amplicon_counts["PARTIAL_ONLY"],
                "AmpliconNoHit": amplicon_counts["NO_HIT"],
                "AmpliconFullLengthSite": amplicon_full_length,

                # Unique-HitID counts used for coverage
                "PCRPositiveHitIDs": denominator,
                "Exact0Mismatch": counts["0_MISMATCH"],
                "OneMismatch": counts["1_MISMATCH"],
                "TwoMismatches": counts["2_MISMATCHES"],
                "GE3Mismatches": counts["GE3_MISMATCHES"],
                "PartialOnly": counts["PARTIAL_ONLY"],
                "NoHit": counts["NO_HIT"],
                "FullLengthSite": full_length,
                "ExactPercentPCRPositive": (
                    100.0 * counts["0_MISMATCH"] / denominator
                    if denominator
                    else 0.0
                ),
                "LE1MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "LE2MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                        + counts["2_MISMATCHES"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "FullLengthPercentPCRPositive": (
                    100.0 * full_length / denominator
                    if denominator
                    else 0.0
                ),
            }
        )

    stats_file = scheme_dir / "probe_variant_statistics.tsv"
    matrix_file = scheme_dir / "probe_variant_hitid_matrix.tsv"

    stats_columns = [
        "QueryVariant",
        "ConcreteSequence",

        "TargetAmplicons",
        "AmpliconExact0Mismatch",
        "AmpliconOneMismatch",
        "AmpliconTwoMismatches",
        "AmpliconGE3Mismatches",
        "AmpliconPartialOnly",
        "AmpliconNoHit",
        "AmpliconFullLengthSite",

        "PCRPositiveHitIDs",
        "Exact0Mismatch",
        "OneMismatch",
        "TwoMismatches",
        "GE3Mismatches",
        "PartialOnly",
        "NoHit",
        "FullLengthSite",
        "ExactPercentPCRPositive",
        "LE1MismatchPercentPCRPositive",
        "LE2MismatchPercentPCRPositive",
        "FullLengthPercentPCRPositive",
    ]

    with stats_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=stats_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(stats_rows)

    variant_ids = [
        str(variant["query_id"])
        for variant in variants
    ]

    with matrix_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )
        writer.writerow(
            ["HitID", *variant_ids, "BestAcrossAllVariants"]
        )

        for hit_id in pcr_positive_hitids:
            categories = [
                matrix[hit_id][variant_id]
                for variant_id in variant_ids
            ]

            best_category = min(
                categories,
                key=lambda category: category_rank[category],
            )

            writer.writerow(
                [
                    hit_id,
                    *categories,
                    best_category,
                ]
            )

    # Display this table when expansion generated more than one concrete probe.
    if VERBOSE and len(stats_rows) > 1:
        print("\nConcrete probe variant statistics — unique PCR-positive HitIDs")
        print(
            f"{'Variant':<25} {'Sequence':<{probe_length + 2}} "
            f"{'Exact':>6} {'1MM':>6} {'2MM':>6} "
            f"{'>=3':>6} {'Partial':>8} {'No hit':>7} {'Full-site':>10}"
        )
        print(
            "-" * (
                25 + probe_length + 2
                + 6 + 6 + 6 + 6 + 8 + 7 + 6 + 8
            )
        )

        for row in stats_rows:
            print(
                f"{str(row['QueryVariant']):<25} "
                f"{str(row['ConcreteSequence']):<{probe_length + 2}} "
                f"{int(row['Exact0Mismatch']):6d} "
                f"{int(row['OneMismatch']):6d} "
                f"{int(row['TwoMismatches']):6d} "
                f"{int(row['GE3Mismatches']):6d} "
                f"{int(row['PartialOnly']):8d} "
                f"{int(row['NoHit']):7d} "
                f"{int(row['FullLengthSite']):6d}"
            )

        print(
            "\nNote: each row is one concrete version of the original "
            "ambiguous probe; targets are unique HitIDs."
        )

    return {
        "stats_rows": stats_rows,
        "stats_file": stats_file,
        "matrix_file": matrix_file,
    }




def short_probe_variant_label(query_id: str) -> str:
    """Convert scheme_PROBE_variant_N to a compact terminal label vN."""
    match = re.search(r"_variant_(\d+)$", query_id)
    if match:
        return f"v{match.group(1)}"
    return query_id


def print_probe_variant_amplicon_table(
    scheme: str,
    probe: str,
    stats_file: Path,
    retained_amplicons: int,
) -> None:
    """
    Display concrete probe versions and their match counts among the retained
    TARGET amplicons that entered the probe-search step.
    """
    if not stats_file.is_file():
        print(f"\nProbe variant statistics unavailable for {scheme}.")
        return

    with stats_file.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    print(f"\nProbe variants — {scheme}")
    print("=" * (18 + len(scheme)))
    print(f"Original IUPAC probe      : {probe}")
    print(f"Retained TARGET amplicons : {retained_amplicons}")
    print(f"Concrete probe variants   : {len(rows)}")

    print("\nConcrete sequences")
    print("------------------")
    print(f"{'Var':<5} Sequence")
    print("-" * max(38, len(probe) + 8))
    for row in rows:
        label = short_probe_variant_label(str(row["QueryVariant"]))
        print(f"{label:<5} {row['ConcreteSequence']}")

    print("\nMatches among retained TARGET amplicons")
    print("----------------------------------------")
    print(
        f"{'Var':<5} {'Exact':>6} {'1MM':>5} {'<=1':>5} "
        f"{'2MM':>5} {'<=2':>5} {'>=3':>5} {'Part':>5} "
        f"{'NoHit':>6} {'Full-site':>9}"
    )
    print("-" * 64)

    for row in rows:
        label = short_probe_variant_label(str(row["QueryVariant"]))
        exact = int(row["AmpliconExact0Mismatch"])
        mm1 = int(row["AmpliconOneMismatch"])
        mm2 = int(row["AmpliconTwoMismatches"])
        mm3 = int(row["AmpliconGE3Mismatches"])
        partial = int(row["AmpliconPartialOnly"])
        no_hit = int(row["AmpliconNoHit"])
        full = int(row["AmpliconFullLengthSite"])

        print(
            f"{label:<5} {exact:>6} {mm1:>5} {exact + mm1:>5} "
            f"{mm2:>5} {exact + mm1 + mm2:>5} {mm3:>5} "
            f"{partial:>5} {no_hit:>6} {full:>9}"
        )

    print(
        f"\nDenominator for every row above: {retained_amplicons} retained "
        "TARGET amplicons."
    )
    print(
        "Exact/1MM/2MM are mutually exclusive categories; <=1 and <=2 are "
        "cumulative counts."
    )
    print(
        "Coverage across the complete validation database is calculated "
        "separately from UNIQUE HitIDs to avoid double-counting biological targets."
    )

def summarize_probe_blast(
    scheme: str,
    probe: str,
    raw_file: Path,
    query_file: Path,
    concrete_query_count: int,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
) -> dict[str, object]:
    """
    Produce BLAST-only summaries.

    There is intentionally:
    - no IUPAC optimization;
    - no ambiguity proposal;
    - no PASS/FAIL probe decision.
    """
    probe_length = len(probe)
    raw_hits = read_probe_blast_raw(raw_file)

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    hits_by_ampid: dict[int, list[dict[str, object]]] = {}

    for hit in raw_hits:
        amp_id = hit["AmpID"]

        if amp_id is None:
            continue

        amp_id = int(amp_id)

        if amp_id not in valid_by_ampid:
            continue

        full_length, category = classify_full_length_probe_hit(
            hit,
            probe_length,
        )

        hit["FullLengthUngapped"] = (
            "YES" if full_length else "NO"
        )
        hit["BlastCategory"] = category

        hits_by_ampid.setdefault(amp_id, []).append(hit)

    category_rank = {
        "0_MISMATCH": 0,
        "1_MISMATCH": 1,
        "2_MISMATCHES": 2,
        "GE3_MISMATCHES": 3,
        "PARTIAL_ONLY": 4,
        "NO_HIT": 5,
    }

    variant_analysis = summarize_probe_variants(
        query_file,
        raw_hits,
        valid_rows,
        probe_length,
        scheme_dir,
        category_rank,
    )

    amplicon_summary: list[dict[str, object]] = []

    for amp_id in sorted(valid_by_ampid):
        source = valid_by_ampid[amp_id]
        amp_hits = hits_by_ampid.get(amp_id, [])

        full_hits = [
            hit
            for hit in amp_hits
            if hit.get("FullLengthUngapped") == "YES"
        ]

        if full_hits:
            full_hits.sort(
                key=lambda hit: (
                    int(hit["mismatch"]),
                    -float(hit["pident"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = full_hits[0]
            category = str(best["BlastCategory"])
        elif amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            category = "PARTIAL_ONLY"
        else:
            best = None
            category = "NO_HIT"

        if best is None:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": "NO",
                "QueryVariant": "",
                "PercentIdentity": "",
                "AlignmentLength": "",
                "Mismatches": "",
                "GapOpenings": "",
                "QueryStart": "",
                "QueryEnd": "",
                "SubjectStart": "",
                "SubjectEnd": "",
                "Evalue": "",
                "BitScore": "",
                "QueryAligned": "",
                "SubjectAligned": "",
                "RawHSPCount": 0,
            }
        else:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": str(
                    best.get("FullLengthUngapped", "NO")
                ),
                "QueryVariant": str(best["qseqid"]),
                "PercentIdentity": float(best["pident"]),
                "AlignmentLength": int(best["length"]),
                "Mismatches": int(best["mismatch"]),
                "GapOpenings": int(best["gapopen"]),
                "QueryStart": int(best["qstart"]),
                "QueryEnd": int(best["qend"]),
                "SubjectStart": int(best["sstart"]),
                "SubjectEnd": int(best["send"]),
                "Evalue": float(best["evalue"]),
                "BitScore": float(best["bitscore"]),
                "QueryAligned": str(best["qseq"]),
                "SubjectAligned": str(best["sseq"]),
                "RawHSPCount": len(amp_hits),
            }

        amplicon_summary.append(row)

    # One best amplicon result per unique target HitID.
    hitid_best: dict[str, dict[str, object]] = {}

    def result_rank(row: dict[str, object]) -> tuple:
        category = str(row["BlastCategory"])

        mismatch_value = (
            int(row["Mismatches"])
            if row["Mismatches"] != ""
            else 999999
        )

        alignment_length = (
            int(row["AlignmentLength"])
            if row["AlignmentLength"] != ""
            else 0
        )

        bitscore = (
            float(row["BitScore"])
            if row["BitScore"] != ""
            else 0.0
        )

        return (
            category_rank[category],
            mismatch_value,
            -alignment_length,
            -bitscore,
            int(row["AmpID"]),
        )

    for row in amplicon_summary:
        hit_id = str(row["HitID"])

        if (
            hit_id not in hitid_best
            or result_rank(row) < result_rank(hitid_best[hit_id])
        ):
            hitid_best[hit_id] = row

    hitid_summary = list(hitid_best.values())
    hitid_summary.sort(
        key=lambda row: str(row["HitID"])
    )

    (
        alignment_amplicon_file,
        alignment_hitid_file,
        mismatch_positions_file,
    ) = write_probe_blast_alignment_reports(
        amplicon_summary,
        hitid_summary,
        scheme_dir,
    )

    counts = {
        category: 0
        for category in category_rank
    }

    for row in hitid_summary:
        counts[str(row["BlastCategory"])] += 1

    amplicon_counts = {
        category: 0
        for category in category_rank
    }

    for row in amplicon_summary:
        amplicon_counts[str(row["BlastCategory"])] += 1

    full_length_amplicons = (
        amplicon_counts["0_MISMATCH"]
        + amplicon_counts["1_MISMATCH"]
        + amplicon_counts["2_MISMATCHES"]
        + amplicon_counts["GE3_MISMATCHES"]
    )

    exact_probe_amplicons = amplicon_counts["0_MISMATCH"]

    full_length_hitids = (
        counts["0_MISMATCH"]
        + counts["1_MISMATCH"]
        + counts["2_MISMATCHES"]
        + counts["GE3_MISMATCHES"]
    )

    exact_probe_hitids = counts["0_MISMATCH"]
    pcr_positive_hitids = len(hitid_summary)
    total_valid_amplicons = len(amplicon_summary)

    # -------------------- write processed outputs --------------------
    columns = [
        "AmpID",
        "HitID",
        "AmpliconSize",
        "BlastCategory",
        "FullLengthUngapped",
        "QueryVariant",
        "PercentIdentity",
        "AlignmentLength",
        "Mismatches",
        "GapOpenings",
        "QueryStart",
        "QueryEnd",
        "SubjectStart",
        "SubjectEnd",
        "Evalue",
        "BitScore",
        "QueryAligned",
        "SubjectAligned",
        "RawHSPCount",
    ]

    amplicon_file = scheme_dir / "probe_blast_amplicon_summary.tsv"
    hitid_file = scheme_dir / "probe_blast_hitid_summary.tsv"
    summary_file = scheme_dir / "probe_blast_summary.txt"

    def write_table(
        path: Path,
        rows: list[dict[str, object]],
    ) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_table(amplicon_file, amplicon_summary)
    write_table(hitid_file, hitid_summary)

    summary_lines = [
        f"BLAST probe summary — {scheme}",
        "=" * (22 + len(scheme)),
        "",
        f"Original VarVAMP probe: {probe}",
        f"Probe length: {probe_length} nt",
        f"Concrete BLAST query variants: {concrete_query_count}",
        f"Valid MFEprimer target amplicons: {len(valid_rows)}",
        f"Unique PCR-positive HitIDs: {pcr_positive_hitids}",
        f"Raw BLAST HSPs: {len(raw_hits)}",
        "",
        "AMPlicon-level probe statistics:",
        f"Valid target amplicons analysed: {total_valid_amplicons}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * exact_probe_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * full_length_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        f"0 mismatch: {amplicon_counts['0_MISMATCH']}",
        f"1 mismatch: {amplicon_counts['1_MISMATCH']}",
        f"2 mismatches: {amplicon_counts['2_MISMATCHES']}",
        f">=3 mismatches: {amplicon_counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {amplicon_counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {amplicon_counts['NO_HIT']}",
        "",
        "HITID-level probe statistics (unique PCR-positive targets):",
        f"Unique PCR-positive HitIDs analysed: {pcr_positive_hitids}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * exact_probe_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        f"0 mismatch: {counts['0_MISMATCH']}",
        f"1 mismatch: {counts['1_MISMATCH']}",
        f"2 mismatches: {counts['2_MISMATCHES']}",
        f">=3 mismatches: {counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {counts['NO_HIT']}",
        "",
        (
            "Full-length ungapped BLAST matches among PCR-positive HitIDs: "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= "
            f"{(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length ungapped probe matches relative to the complete "
            f"validation database: {full_length_hitids} / "
            f"{total_target_sequences} = "
            f"{(100.0 * full_length_hitids / total_target_sequences if total_target_sequences else 0.0):.2f} %"
        ),
        "",
        "Full-length ungapped definition:",
        f"- query starts at 1 and ends at {probe_length};",
        f"- alignment length = {probe_length} nt;",
        "- zero gap openings.",
        "",
        "This version reports BLAST results only.",
        "No probe optimization or IUPAC substitution is proposed.",
        f"Raw BLAST output: {compact_path(raw_file)}",
        f"BLAST query FASTA: {compact_path(query_file)}",
        f"Visual alignments by amplicon: {alignment_amplicon_file}",
        f"Visual alignments by HitID: {alignment_hitid_file}",
        f"Mismatch positions TSV: {mismatch_positions_file}",
        f"Per-variant statistics TSV: {variant_analysis['stats_file']}",
        f"HitID x probe-variant matrix TSV: {variant_analysis['matrix_file']}",
    ]

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"Probe results — {scheme}")

    def _pct(value: int, total: int) -> float:
        return 100.0 * value / total if total else 0.0

    print("Biological coverage — unique PCR-positive targets")
    print(f"  PCR-positive HitIDs   {pcr_positive_hitids}")
    print(
        f"  Exact probe           {counts['0_MISMATCH']:>4}  "
        f"{_pct(counts['0_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  1 mismatch            {counts['1_MISMATCH']:>4}  "
        f"{_pct(counts['1_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  2 mismatches          {counts['2_MISMATCHES']:>4}  "
        f"{_pct(counts['2_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  >=3 mismatches        {counts['GE3_MISMATCHES']:>4}  "
        f"{_pct(counts['GE3_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Partial/gapped        {counts['PARTIAL_ONLY']:>4}  "
        f"{_pct(counts['PARTIAL_ONLY'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  No hit                {counts['NO_HIT']:>4}  "
        f"{_pct(counts['NO_HIT'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Full-length site      {full_length_hitids:>4}  "
        f"{_pct(full_length_hitids, pcr_positive_hitids):6.2f} %"
    )

    print("\nTechnical amplicon statistics")
    print(f"  Valid amplicons       {total_valid_amplicons}")
    print(f"  Exact probe           {amplicon_counts['0_MISMATCH']}")
    print(f"  1 mismatch            {amplicon_counts['1_MISMATCH']}")
    print(f"  2 mismatches          {amplicon_counts['2_MISMATCHES']}")
    print(f"  >=3 mismatches        {amplicon_counts['GE3_MISMATCHES']}")
    print(f"  Partial/gapped        {amplicon_counts['PARTIAL_ONLY']}")
    print(f"  No hit                {amplicon_counts['NO_HIT']}")

    if VERBOSE:
        print(f"\nRaw BLAST               {raw_file}")
        print(f"Amplicon summary         {amplicon_file}")
        print(f"HitID summary            {hitid_file}")
        print(f"Alignment by amplicon    {alignment_amplicon_file}")
        print(f"Alignment by HitID       {alignment_hitid_file}")
        print(f"Mismatch positions       {mismatch_positions_file}")
        print(f"Variant statistics       {variant_analysis['stats_file']}")
        print(f"HitID x variant matrix   {variant_analysis['matrix_file']}")
        print(f"BLAST summary            {summary_file}")

    return {
        "raw_hsps": len(raw_hits),
        "amplicon_total": total_valid_amplicons,
        "amplicon_exact_0": exact_probe_amplicons,
        "amplicon_full_length": full_length_amplicons,
        "amplicon_1_mismatch": amplicon_counts["1_MISMATCH"],
        "amplicon_2_mismatch": amplicon_counts["2_MISMATCHES"],
        "amplicon_ge3_mismatch": amplicon_counts["GE3_MISMATCHES"],
        "exact_0": counts["0_MISMATCH"],
        "mismatch_1": counts["1_MISMATCH"],
        "mismatch_2": counts["2_MISMATCHES"],
        "mismatch_ge3": counts["GE3_MISMATCHES"],
        "partial_only": counts["PARTIAL_ONLY"],
        "no_hit": counts["NO_HIT"],
        "full_length_hitids": full_length_hitids,
        "pcr_positive_hitids": pcr_positive_hitids,
        "amplicon_summary": amplicon_file,
        "hitid_summary": hitid_file,
        "blast_summary": summary_file,
        "alignment_amplicon_file": alignment_amplicon_file,
        "alignment_hitid_file": alignment_hitid_file,
        "mismatch_positions_file": mismatch_positions_file,
        "variant_statistics_file": variant_analysis["stats_file"],
        "variant_hitid_matrix_file": variant_analysis["matrix_file"],
    }


def analyse_probe_with_blast(
    scheme: str,
    probe: str,
    valid_amplicons_fasta: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
) -> dict[str, object]:
    """Run and summarize BLAST analysis for the original VarVAMP probe."""
    raw_file, query_file, variant_count = run_probe_blast(
        scheme,
        probe,
        valid_amplicons_fasta,
        valid_rows,
        scheme_dir,
    )

    return summarize_probe_blast(
        scheme,
        probe,
        raw_file,
        query_file,
        variant_count,
        valid_rows,
        scheme_dir,
        total_target_sequences,
    )





def parse_mfeprimer_secondary_structures(report_file: Path) -> dict[str, object]:
    """
    Parse Hairpin and Dimer sections from a legacy MFEprimer 3.x text report.

    Reported counts are authoritative. Detailed structure blocks are parsed
    when their textual format is recognized. If MFEprimer reports structures
    but their detail format is not recognized, the counts are still preserved
    and the raw report remains the authoritative source.
    """
    if not report_file.is_file():
        raise FileNotFoundError(f"MFEprimer report not found: {report_file}")

    lines = report_file.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count_patterns = {
        "hairpin": re.compile(r"^\s*Hairpin\s+List\s*\(\s*(\d+)\s*\)\s*$", re.IGNORECASE),
        "dimer": re.compile(r"^\s*Dimer\s+List\s*\(\s*(\d+)\s*\)\s*$", re.IGNORECASE),
    }

    counts: dict[str, int | None] = {"hairpin": None, "dimer": None}
    section_start: dict[str, int | None] = {"hairpin": None, "dimer": None}

    for index, line in enumerate(lines):
        for kind, pattern in count_patterns.items():
            match = pattern.match(line)
            if match:
                counts[kind] = int(match.group(1))
                section_start[kind] = index + 1

    # MFEprimer 3.x reports normally contain both sections. Keep "not found"
    # distinct from a true zero so the workflow does not silently invent data.
    hairpin_count = counts["hairpin"]
    dimer_count = counts["dimer"]

    def section_lines(kind: str) -> list[str]:
        start = section_start[kind]
        if start is None:
            return []

        end = len(lines)
        if kind == "hairpin" and section_start["dimer"] is not None:
            end = int(section_start["dimer"]) - 1
        elif kind == "dimer":
            for idx in range(start, len(lines)):
                if re.search(
                    r"Descriptions\s+of\s+\[\s*\d+\s*\]\s+potential\s+amplicons",
                    lines[idx],
                    flags=re.IGNORECASE,
                ):
                    end = idx
                    break
        return lines[start:end]

    def parse_numeric(block: list[str], patterns: list[re.Pattern[str]]) -> float | None:
        joined = "\n".join(block)
        for pattern in patterns:
            match = pattern.search(joined)
            if match:
                try:
                    return float(match.group(1))
                except ValueError:
                    return None
        return None

    score_patterns = [
        re.compile(r"\bScore\s*[:=]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE),
    ]
    tm_patterns = [
        re.compile(r"\bTm\s*[:=]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE),
        re.compile(
            r"\bmelting\s+temperature\s*[:=]\s*(-?\d+(?:\.\d+)?)",
            re.IGNORECASE,
        ),
    ]
    dg_patterns = [
        re.compile(
            r"(?:Delta\s*G|Dg|ΔG)\s*[:=]\s*(-?\d+(?:\.\d+)?)",
            re.IGNORECASE,
        ),
    ]

    def parse_items(kind: str) -> list[dict[str, object]]:
        block_lines = section_lines(kind)
        if not block_lines:
            return []

        label = "Hairpin" if kind == "hairpin" else "Dimer"
        item_re = re.compile(
            rf"^\s*{label}\s+(\d+)\s*:\s*(.*?)\s*$",
            re.IGNORECASE,
        )

        starts: list[tuple[int, re.Match[str]]] = []
        for idx, raw in enumerate(block_lines):
            match = item_re.match(raw)
            if match:
                starts.append((idx, match))

        items: list[dict[str, object]] = []
        for pos, (start_idx, match) in enumerate(starts):
            end_idx = starts[pos + 1][0] if pos + 1 < len(starts) else len(block_lines)
            item_block = block_lines[start_idx:end_idx]
            subject = match.group(2).strip()

            oligo_1 = subject
            oligo_2 = ""
            if kind == "dimer":
                split = re.split(r"\s+(?:x|×)\s+", subject, maxsplit=1, flags=re.IGNORECASE)
                if len(split) == 2:
                    oligo_1, oligo_2 = split[0].strip(), split[1].strip()

            items.append(
                {
                    "index": int(match.group(1)),
                    "subject": subject,
                    "oligo_1": oligo_1,
                    "oligo_2": oligo_2,
                    "score": parse_numeric(item_block, score_patterns),
                    "tm_c": parse_numeric(item_block, tm_patterns),
                    "delta_g": parse_numeric(item_block, dg_patterns),
                    "raw_block": "\n".join(item_block).strip(),
                }
            )
        return items

    hairpins = parse_items("hairpin")
    dimers = parse_items("dimer")

    return {
        "hairpin_section_found": hairpin_count is not None,
        "dimer_section_found": dimer_count is not None,
        "hairpin_count": hairpin_count,
        "dimer_count": dimer_count,
        "hairpins": hairpins,
        "dimers": dimers,
    }


def print_secondary_structure_qc(
    secondary_by_scheme: dict[str, dict[str, object]],
    selected_schemes: list[str],
    assay_mode: str,
) -> None:
    """Display MFEprimer hairpin/dimer data only when those report sections exist."""
    available = [
        scheme
        for scheme in selected_schemes
        if secondary_by_scheme[scheme]["hairpin_count"] is not None
        or secondary_by_scheme[scheme]["dimer_count"] is not None
    ]

    print("\nSecondary-structure QC — MFEprimer")
    print("----------------------------------")

    if not available:
        print(
            "Hairpin/dimer sections are not present in these MFEprimer specificity "
            "reports, so no secondary-structure conclusion is drawn here."
        )
        if assay_mode == "qpcr":
            print(
                "Only LEFT/RIGHT primers are present in the MFEprimer QC input. "
                "Probe hairpins/dimers are therefore not assessed here; the probe is "
                "evaluated separately by BLAST in the probe-validation step."
            )
        return

    print(f"{'Scheme':<16} {'Hairpins':>10} {'Dimers':>10}")
    print("-" * 38)

    for scheme in selected_schemes:
        result = secondary_by_scheme[scheme]
        hairpins = (
            str(result["hairpin_count"])
            if result["hairpin_count"] is not None
            else "N/A"
        )
        dimers = (
            str(result["dimer_count"])
            if result["dimer_count"] is not None
            else "N/A"
        )
        print(f"{scheme:<16} {hairpins:>10} {dimers:>10}")

    print(
        "\nCounts mean structures reported by MFEprimer under the selected "
        "analysis parameters; zero does not prove that a structure can never form."
    )

    positive_found = False
    for scheme in selected_schemes:
        result = secondary_by_scheme[scheme]
        reported_positive = (
            (result["hairpin_count"] or 0) > 0
            or (result["dimer_count"] or 0) > 0
        )
        if not reported_positive:
            continue

        positive_found = True
        print(f"\nDetected structures — {scheme}")
        print("-" * (23 + len(scheme)))

        details = list(result["hairpins"]) + list(result["dimers"])
        if not details:
            print(
                "MFEprimer reported one or more structures, but the detailed "
                "block could not be parsed. See the raw MFEprimer report."
            )
            continue

        print(
            f"{'Type':<8} {'Oligo(s)':<31} {'Score':>7} "
            f"{'Tm°C':>8} {'ΔG':>9}"
        )
        print("-" * 68)

        for kind, items in (
            ("Hairpin", result["hairpins"]),
            ("Dimer", result["dimers"]),
        ):
            for item in items:
                if kind == "Dimer" and item["oligo_2"]:
                    oligos = f"{item['oligo_1']} × {item['oligo_2']}"
                else:
                    oligos = str(item["oligo_1"])

                score = (
                    f"{float(item['score']):.2f}"
                    if item["score"] is not None
                    else "N/A"
                )
                tm_c = (
                    f"{float(item['tm_c']):.2f}"
                    if item["tm_c"] is not None
                    else "N/A"
                )
                delta_g = (
                    f"{float(item['delta_g']):.2f}"
                    if item["delta_g"] is not None
                    else "N/A"
                )
                print(
                    f"{kind:<8} {oligos[:31]:<31} "
                    f"{score:>7} {tm_c:>8} {delta_g:>9}"
                )

    if not positive_found:
        print("\nNo hairpins or dimers were reported in the available QC sections.")


def save_secondary_structure_qc(
    coverage_root: Path,
    secondary_by_scheme: dict[str, dict[str, object]],
    selected_schemes: list[str],
) -> tuple[Path, Path]:
    """Save MFEprimer hairpin/dimer summary and any parsed detailed structures."""
    summary_file = coverage_root / "secondary_structure_summary.tsv"
    details_file = coverage_root / "secondary_structure_details.tsv"

    with summary_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "hairpin_section_found",
                "hairpins_reported",
                "hairpin_details_parsed",
                "dimer_section_found",
                "dimers_reported",
                "dimer_details_parsed",
            ]
        )
        for scheme in selected_schemes:
            result = secondary_by_scheme[scheme]
            writer.writerow(
                [
                    scheme,
                    result["hairpin_section_found"],
                    "" if result["hairpin_count"] is None else result["hairpin_count"],
                    len(result["hairpins"]),
                    result["dimer_section_found"],
                    "" if result["dimer_count"] is None else result["dimer_count"],
                    len(result["dimers"]),
                ]
            )

    with details_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "structure_type",
                "structure_index",
                "oligo_1",
                "oligo_2",
                "score",
                "tm_c",
                "delta_g",
                "raw_block",
            ]
        )
        for scheme in selected_schemes:
            result = secondary_by_scheme[scheme]
            for kind, items in (
                ("HAIRPIN", result["hairpins"]),
                ("DIMER", result["dimers"]),
            ):
                for item in items:
                    writer.writerow(
                        [
                            scheme,
                            kind,
                            item["index"],
                            item["oligo_1"],
                            item["oligo_2"],
                            "" if item["score"] is None else item["score"],
                            "" if item["tm_c"] is None else item["tm_c"],
                            "" if item["delta_g"] is None else item["delta_g"],
                            item["raw_block"],
                        ]
                    )

    return summary_file, details_file

def parse_mfeprimer_primer_properties(report_file: Path) -> list[dict[str, object]]:
    """Parse the compact primer-property rows from a legacy MFEprimer report."""
    row_re = re.compile(
        r"^\s*(\S+)\s+([ACGTRYSWKMBDHVNacgtryswkmbdhvn]+)\s+"
        r"(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+(-?[0-9.]+)\s+"
        r"(\d+)\s+(\d+)\s*$"
    )
    rows: list[dict[str, object]] = []
    for raw in report_file.read_text(encoding="utf-8", errors="replace").splitlines():
        match = row_re.match(raw)
        if not match:
            continue
        rows.append({
            "id": match.group(1),
            "sequence": match.group(2).upper(),
            "length": int(match.group(3)),
            "gc": float(match.group(4)),
            "tm": float(match.group(5)),
            "dg": float(match.group(6)),
            "bind_plus": int(match.group(7)),
            "bind_minus": int(match.group(8)),
        })
    return rows


def compact_primer_label(primer_id: str, scheme: str) -> str:
    """Remove a repeated scheme prefix while preserving LEFT/RIGHT variant labels."""
    label = primer_id
    prefix_re = re.compile(rf"^{re.escape(scheme)}[_-]?", flags=re.IGNORECASE)
    label = prefix_re.sub("", label)
    return label or primer_id


def print_compact_primer_table(report_file: Path, scheme: str) -> None:
    """Display a terminal-width friendly MFEprimer primer table."""
    rows = parse_mfeprimer_primer_properties(report_file)
    if not rows:
        if VERBOSE:
            print(f"No primer-property rows parsed from {report_file.name}.")
        return

    print(f"\nConcrete primer variants — {scheme}")
    print("-" * (27 + len(scheme)))
    print(
        f"{'Primer':<12} {'Sequence':<30} {'Len':>3} "
        f"{'GC%':>6} {'Tm°C':>7} {'ΔG':>8} {'Bind+':>7} {'Bind-':>7}"
    )
    print("-" * 86)
    for row in rows:
        label = compact_primer_label(str(row['id']), scheme)
        print(
            f"{label:<12} {str(row['sequence']):<30} {int(row['length']):>3} "
            f"{float(row['gc']):>6.2f} {float(row['tm']):>7.2f} "
            f"{float(row['dg']):>8.2f} {int(row['bind_plus']):>7} "
            f"{int(row['bind_minus']):>7}"
        )


def inspect_mfeprimer_products(report_file: Path, scheme: str) -> dict[str, object]:
    """
    Inspect ALL MFEprimer potential products before any workflow filtering.

    The returned structure keeps:
    - total potential amplicons;
    - product classes;
    - unique HitIDs represented in each class;
    - exact concrete primer-pair combinations.
    """
    parsed = parse_mfeprimer3_report(report_file)

    classes = ("TARGET", "SELF_PRIMING", "CROSS_SCHEME", "OTHER")
    class_counts = {key: 0 for key in classes}
    class_hitids: dict[str, set[str]] = {key: set() for key in classes}

    pair_counts = {
        key: 0
        for key in ("LEFT×RIGHT", "RIGHT×LEFT", "LEFT×LEFT", "RIGHT×RIGHT")
    }

    all_hitids: set[str] = set()
    concrete: dict[tuple[str, str, str], dict[str, object]] = {}

    for item in parsed["amplicons"]:
        primer_1 = str(item["primer_1"])
        primer_2 = str(item["primer_2"])
        hit_id = str(item["hit_id"])

        product_class, role_1, role_2 = classify_mfeprimer_product(
            primer_1,
            primer_2,
            scheme,
        )

        class_counts[product_class] = class_counts.get(product_class, 0) + 1
        class_hitids.setdefault(product_class, set()).add(hit_id)
        all_hitids.add(hit_id)

        role_pair = f"{role_1}×{role_2}"
        if role_pair in pair_counts:
            pair_counts[role_pair] += 1

        label_1 = compact_primer_label(primer_1, scheme)
        label_2 = compact_primer_label(primer_2, scheme)
        key = (label_1, label_2, product_class)

        if key not in concrete:
            concrete[key] = {
                "primer_1": label_1,
                "primer_2": label_2,
                "class": product_class,
                "amplicons": 0,
                "hitids": set(),
            }

        concrete[key]["amplicons"] = int(concrete[key]["amplicons"]) + 1
        concrete[key]["hitids"].add(hit_id)

    concrete_rows: list[dict[str, object]] = []
    for item in concrete.values():
        concrete_rows.append(
            {
                "primer_1": item["primer_1"],
                "primer_2": item["primer_2"],
                "class": item["class"],
                "amplicons": int(item["amplicons"]),
                "unique_hitids": len(item["hitids"]),
            }
        )

    concrete_rows.sort(
        key=lambda row: (
            0 if row["class"] == "TARGET" else
            1 if row["class"] == "SELF_PRIMING" else 2,
            str(row["primer_1"]),
            str(row["primer_2"]),
        )
    )

    return {
        "potential": len(parsed["amplicons"]),
        "all_unique_hitids": len(all_hitids),
        "class_counts": class_counts,
        "class_unique_hitids": {
            key: len(class_hitids.get(key, set()))
            for key in classes
        },
        "pair_counts": pair_counts,
        "concrete_pairs": concrete_rows,
    }



def print_general_potential_amplicon_statistics(
    previews: dict[str, dict[str, object]],
    selected_schemes: list[str],
    total_target_sequences: int,
) -> None:
    """
    Display the complete MFEprimer product picture BEFORE any workflow filter.
    """
    print("\nMFEprimer potential-product overview — BEFORE filtering")
    print("-------------------------------------------------------")
    print(
        "These are all potential products reported by MFEprimer. "
        "No TARGET/SELF_PRIMING exclusion has been applied yet."
    )

    for scheme in selected_schemes:
        preview = previews[scheme]
        counts = preview["class_counts"]
        hitids = preview["class_unique_hitids"]

        any_hitids = int(preview["all_unique_hitids"])
        any_pct = (
            100.0 * any_hitids / total_target_sequences
            if total_target_sequences else 0.0
        )

        print(f"\n{scheme}")
        print("=" * len(scheme))
        print(f"Validation database                 : {total_target_sequences} sequences")
        print(f"Potential amplicons reported        : {int(preview['potential'])}")
        print(
            f"Sequences with >=1 potential product: "
            f"{any_hitids}/{total_target_sequences} ({any_pct:.2f} %)"
        )

        print("\nProduct classes")
        print("---------------")
        print(
            f"{'Class':<18} {'Amplicons':>10} "
            f"{'Unique HitIDs':>14} {'% database':>12}"
        )
        print("-" * 58)

        for product_class, label in (
            ("TARGET", "TARGET"),
            ("SELF_PRIMING", "SELF_PRIMING"),
            ("CROSS_SCHEME", "CROSS_SCHEME"),
            ("OTHER", "OTHER"),
        ):
            n_amp = int(counts.get(product_class, 0))
            n_hits = int(hitids.get(product_class, 0))
            pct = (
                100.0 * n_hits / total_target_sequences
                if total_target_sequences else 0.0
            )
            print(
                f"{label:<18} {n_amp:>10} "
                f"{n_hits:>14} {pct:>11.2f}%"
            )

        print(
            "\nTARGET = LEFT×RIGHT or RIGHT×LEFT. "
            "SELF_PRIMING = LEFT×LEFT or RIGHT×RIGHT."
        )
        print(
            "Note: class-specific unique HitID counts are not additive; the same "
            "sequence may generate products in more than one class."
        )

def print_potential_amplicons_by_pair(
    previews: dict[str, dict[str, object]],
    selected_schemes: list[str],
) -> None:
    """Show potential MFEprimer amplicons at the exact primer-pair level."""
    print("\nPotential amplicons by primer pair")
    print("----------------------------------")
    print(
        "Counts below are raw MFEprimer potential products before the workflow "
        "classification/filtering step."
    )

    for scheme in selected_schemes:
        preview = previews[scheme]
        rows = list(preview.get("concrete_pairs", []))

        print(f"\n{scheme}")
        print("-" * len(scheme))
        print(
            f"{'Primer pair':<31} {'Class':<14} "
            f"{'Amplicons':>10} {'HitIDs':>8}"
        )
        print("-" * 67)

        if not rows:
            print("(no potential amplicons)")
            continue

        for row in rows:
            pair_label = f"{row['primer_1']} × {row['primer_2']}"
            print(
                f"{pair_label:<31} {str(row['class']):<14} "
                f"{int(row['amplicons']):>10} "
                f"{int(row['unique_hitids']):>8}"
            )

        print(
            f"{'TOTAL':<31} {'':<14} "
            f"{int(preview['potential']):>10}"
        )


def empty_probe_analysis() -> dict[str, object]:
    """Return a schema-compatible empty probe result."""
    return {
        "amplicon_total": 0,
        "amplicon_exact_0": 0,
        "amplicon_full_length": 0,
        "amplicon_1_mismatch": 0,
        "amplicon_2_mismatch": 0,
        "amplicon_ge3_mismatch": 0,
        "exact_0": 0,
        "mismatch_1": 0,
        "mismatch_2": 0,
        "mismatch_ge3": 0,
        "partial_only": 0,
        "no_hit": 0,
        "full_length_hitids": 0,
        "pcr_positive_hitids": 0,
        "hitid_summary": "",
        "blast_summary": "",
        "alignment_hitid_file": "",
        "mismatch_positions_file": "",
        "variant_statistics_file": "",
        "variant_hitid_matrix_file": "",
    }

def run_mfeprimer_pairs(
    schemes: dict[str, dict[str, str]],
    selected_schemes: list[str],
    assay_inputs_dir: Path,
    mfeprimer_results_dir: Path,
    probe_results_dir: Path,
    summary_results_dir: Path,
    database: Path,
    parameters: dict[str, float | int] | None,
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
    assay_mode: str,
) -> dict[str, object]:
    """
    Run validation in organized phases.

    All selected primer pairs are processed first. Product classes are then
    compared together, one global SELF_PRIMING decision is made only when
    necessary, and qPCR probes are handled in a separate batch phase.
    """
    mfeprimer_results_dir.mkdir(parents=True, exist_ok=True)
    probe_results_dir.mkdir(parents=True, exist_ok=True)
    summary_results_dir.mkdir(parents=True, exist_ok=True)

    run_manifest: list[list[str]] = []
    coverage_manifest: list[list[object]] = []
    reports: dict[str, Path] = {}
    scheme_dirs: dict[str, Path] = {}
    previews: dict[str, dict[str, object]] = {}
    secondary_by_scheme: dict[str, dict[str, object]] = {}

    print("\n8. Primer-pair validation")
    print("=========================")
    print(f"Database : {compact_path(database)} (N={total_target_sequences})")
    print(f"Assays   : {len(selected_schemes)}")

    # Phase A — run MFEprimer for every selected assay before asking downstream questions.
    for scheme in selected_schemes:
        pair_file = assay_inputs_dir / "primer_pairs" / f"{scheme}_primers.fasta"
        scheme_dir = mfeprimer_results_dir / scheme
        scheme_dir.mkdir(parents=True, exist_ok=True)
        scheme_dirs[scheme] = scheme_dir
        output_prefix = scheme_dir / scheme

        # Run FULL MFEprimer quality control.
        # Do not use the `spec` subcommand here: `spec` performs specificity
        # only, whereas the top-level `mfeprimer` command also reports
        # hairpins and dimers in MFEprimer 3.x.
        command = [
            "mfeprimer",
            "-i", str(pair_file),
            "-d", str(database),
            "-o", str(output_prefix),
        ]
        if parameters is not None:
            command.extend([
                "-s", str(parameters["min_size"]),
                "-S", str(parameters["max_size"]),
                "-t", str(parameters["tm_cutoff"]),
            ])

        progress(f"MFEprimer full QC {scheme}", "RUN")
        run_command(command)

        expected_tsv = Path(str(output_prefix) + ".spec.tsv")
        legacy_output = output_prefix
        if legacy_output.is_file():
            result_path = legacy_output
        elif expected_tsv.is_file():
            result_path = expected_tsv
        else:
            raise RuntimeError(
                "MFEprimer completed but its output file could not be found "
                f"for scheme {scheme}."
            )

        if result_path == expected_tsv:
            raise RuntimeError(
                "A .spec.tsv result was produced, but amplicon-sequence "
                "extraction in this workflow currently requires the legacy "
                "MFEprimer 3.x text report."
            )

        reports[scheme] = result_path
        previews[scheme] = inspect_mfeprimer_products(result_path, scheme)
        secondary_by_scheme[scheme] = parse_mfeprimer_secondary_structures(result_path)

        secondary_qc = secondary_by_scheme[scheme]
        if (
            not bool(secondary_qc["hairpin_section_found"])
            or not bool(secondary_qc["dimer_section_found"])
        ):
            raise RuntimeError(
                "MFEprimer full QC completed, but the expected Hairpin List "
                "and/or Dimer List section is missing from the legacy text report "
                f"for {scheme}. The workflow stops BEFORE product filtering so an "
                "incomplete QC report is not silently accepted. Raw report: "
                f"{result_path}"
            )

        run_manifest.append(
            [
                scheme,
                portable_project_path(pair_file),
                portable_project_path(database),
                portable_project_path(result_path),
            ]
        )
        progress(f"MFEprimer full QC {scheme}", "OK")

    # Show concrete primer variants first, immediately followed by the
    # number of potential products reported by MFEprimer for that scheme.
    for scheme in selected_schemes:
        print_compact_primer_table(reports[scheme], scheme)
        print(
            f"Potential amplicons reported by MFEprimer: "
            f"{int(previews[scheme]['potential'])}"
        )

    print_secondary_structure_qc(
        secondary_by_scheme,
        selected_schemes,
        assay_mode,
    )
    secondary_summary_file, secondary_details_file = save_secondary_structure_qc(
        summary_results_dir,
        secondary_by_scheme,
        selected_schemes,
    )

    # First show the complete product landscape BEFORE any filtering.
    print_general_potential_amplicon_statistics(
        previews,
        selected_schemes,
        total_target_sequences,
    )
    print_potential_amplicons_by_pair(previews, selected_schemes)

    print("\nFilter decision")
    print("---------------")
    print(
        "The downstream PCR-coverage/probe step can retain TARGET products only, "
        "or keep TARGET + SELF_PRIMING products."
    )
    print("TARGET       : LEFT×RIGHT / RIGHT×LEFT")
    print("SELF_PRIMING : LEFT×LEFT / RIGHT×RIGHT")
    print(
        "CROSS_SCHEME / OTHER are never treated as valid products for the "
        "selected assay and remain excluded."
    )

    exclude_self_priming = ask_yes_no(
        "Apply the TARGET-only filter before coverage and probe analysis?"
    )

    # Apply the user's explicit decision only after the unfiltered MFEprimer
    # results have been displayed.
    analyses: dict[str, dict[str, object]] = {}
    for scheme in selected_schemes:
        if VERBOSE:
            analysis = analyse_and_extract_mfeprimer_amplicons(
                reports[scheme],
                scheme,
                scheme_dirs[scheme],
                total_target_sequences,
                min_target_size,
                max_target_size,
                exclude_self_priming=exclude_self_priming,
            )
        else:
            with redirect_stdout(StringIO()):
                analysis = analyse_and_extract_mfeprimer_amplicons(
                    reports[scheme],
                    scheme,
                    scheme_dirs[scheme],
                    total_target_sequences,
                    min_target_size,
                    max_target_size,
                    exclude_self_priming=exclude_self_priming,
                )
        analyses[scheme] = analysis

    print("\nResults AFTER filtering")
    print("-----------------------")
    print(
        "Filter applied: "
        + (
            "TARGET only (SELF_PRIMING excluded)"
            if exclude_self_priming
            else "TARGET + SELF_PRIMING retained"
        )
    )

    print(
        f"\n{'Scheme':<16} {'Potential':>9} {'TARGET':>8} "
        f"{'SELF':>7} {'Retained':>9} {'Removed':>8} "
        f"{'PCR+':>7} {'Coverage':>10}"
    )
    print("-" * 84)

    retained_assays = 0
    for scheme in selected_schemes:
        preview = previews[scheme]
        analysis = analyses[scheme]

        potential = int(preview["potential"])
        target = int(preview["class_counts"].get("TARGET", 0))
        self_products = int(preview["class_counts"].get("SELF_PRIMING", 0))
        retained = int(analysis["valid_target_amplicons"])
        removed = potential - retained
        pcr_positive = int(analysis["unique_positive_hitids"])

        if pcr_positive > 0:
            retained_assays += 1

        print(
            f"{scheme:<16} {potential:>9} {target:>8} "
            f"{self_products:>7} {retained:>9} {removed:>8} "
            f"{pcr_positive:>7} "
            f"{float(analysis['coverage_pct']):>9.2f}%"
        )

    if len(selected_schemes) > 1:
        print(
            f"\nAssays with at least one retained PCR-positive target: "
            f"{retained_assays}/{len(selected_schemes)}"
        )
    else:
        scheme = selected_schemes[0]
        analysis = analyses[scheme]
        print(
            f"\nRetained PCR-positive targets: "
            f"{int(analysis['unique_positive_hitids'])}/"
            f"{total_target_sequences} "
            f"({float(analysis['coverage_pct']):.2f} %)"
        )

    if VERBOSE:
        print(
            "\nCoverage set: "
            + (
                "TARGET only (SELF_PRIMING excluded)."
                if exclude_self_priming
                else "TARGET + SELF_PRIMING."
            )
        )

    # Phase C — one probe decision, then probe batch.
    probe_results: dict[str, dict[str, object]] = {
        scheme: empty_probe_analysis() for scheme in selected_schemes
    }
    probe_validation_enabled = False

    if assay_mode == "qpcr":
        eligible = [
            scheme for scheme in selected_schemes
            if int(analyses[scheme]["valid_target_amplicons"]) > 0
        ]
        print("\n9. Probe validation")
        print("===================")
        if not eligible:
            print("No retained PCR amplicons; probe validation skipped.")
        else:
            probe_validation_enabled = ask_yes_no(
                "Test the existing VarVAMP probe(s) on the retained PCR products?"
            )
            if probe_validation_enabled:
                for scheme in eligible:
                    progress(f"Probe BLAST {scheme}", "RUN")
                    probe_scheme_dir = probe_results_dir / scheme
                    probe_scheme_dir.mkdir(parents=True, exist_ok=True)
                    if VERBOSE:
                        blast_analysis = analyse_probe_with_blast(
                            scheme,
                            schemes[scheme]["PROBE"],
                            analyses[scheme]["valid_fasta"],
                            analyses[scheme]["valid_rows"],
                            probe_scheme_dir,
                            total_target_sequences,
                        )
                    else:
                        with redirect_stdout(StringIO()):
                            blast_analysis = analyse_probe_with_blast(
                                scheme,
                                schemes[scheme]["PROBE"],
                                analyses[scheme]["valid_fasta"],
                                analyses[scheme]["valid_rows"],
                                probe_scheme_dir,
                                total_target_sequences,
                            )
                    probe_results[scheme] = blast_analysis
                    progress(f"Probe BLAST {scheme}", "OK")
            else:
                print("Probe validation skipped by user.")

        if probe_validation_enabled:
            # First show exactly what entered probe validation and how EACH
            # concrete IUPAC-expanded probe version behaved against the retained
            # TARGET amplicons.
            for scheme in eligible:
                result = probe_results[scheme]
                print_probe_variant_amplicon_table(
                    scheme,
                    schemes[scheme]["PROBE"],
                    Path(str(result["variant_statistics_file"])),
                    int(analyses[scheme]["valid_target_amplicons"]),
                )

            # Then show biological coverage after choosing the best result across
            # all concrete versions for each UNIQUE PCR-positive HitID.
            print("\nBest probe result across all concrete variants")
            print("----------------------------------------------")
            print(
                f"{'Scheme':<16} {'PCR+':>6} {'Exact':>7} "
                f"{'<=1MM':>7} {'<=2MM':>7} {'Full-site':>10}"
            )
            print("-" * 58)

            probe_positive_assays = 0
            for scheme in selected_schemes:
                result = probe_results[scheme]
                pcr_positive = int(result["pcr_positive_hitids"])
                exact = int(result["exact_0"])
                le1 = exact + int(result["mismatch_1"])
                le2 = le1 + int(result["mismatch_2"])
                full = int(result["full_length_hitids"])

                if full > 0:
                    probe_positive_assays += 1

                print(
                    f"{scheme:<16} {pcr_positive:>6} "
                    f"{exact:>7} {le1:>7} {le2:>7} {full:>7}"
                )

            print(
                "\nThis second table is deduplicated by UNIQUE HitID and is the "
                "one used for assay coverage."
            )

            if len(selected_schemes) > 1:
                print(
                    f"Assays with at least one PCR-positive HitID carrying a "
                    f"full-length probe site: "
                    f"{probe_positive_assays}/{len(selected_schemes)}"
                )

    # Build the same machine-readable coverage manifest as before.
    for scheme in selected_schemes:
        analysis = analyses[scheme]
        blast_analysis = probe_results[scheme]
        coverage_manifest.append([
            assay_mode,
            scheme,
            analysis["potential_amplicons"],
            analysis["target_products"],
            analysis["self_priming_products"],
            analysis["valid_target_amplicons"],
            analysis["unique_positive_hitids"],
            f"{analysis['coverage_pct']:.4f}",
            blast_analysis["amplicon_total"],
            blast_analysis["amplicon_exact_0"],
            blast_analysis["amplicon_full_length"],
            blast_analysis["amplicon_1_mismatch"],
            blast_analysis["amplicon_2_mismatch"],
            blast_analysis["amplicon_ge3_mismatch"],
            blast_analysis["exact_0"],
            blast_analysis["mismatch_1"],
            blast_analysis["mismatch_2"],
            blast_analysis["mismatch_ge3"],
            blast_analysis["partial_only"],
            blast_analysis["no_hit"],
            blast_analysis["full_length_hitids"],
            portable_project_path(Path(analysis["valid_fasta"])),
            portable_project_path(Path(analysis["positive_hitids"])),
            portable_project_path(Path(analysis["summary"])),
            optional_portable_project_path(blast_analysis["hitid_summary"]),
            optional_portable_project_path(blast_analysis["blast_summary"]),
            optional_portable_project_path(blast_analysis["alignment_hitid_file"]),
            optional_portable_project_path(blast_analysis["mismatch_positions_file"]),
            (
                f"{(100.0 * blast_analysis['exact_0'] / total_target_sequences if total_target_sequences else 0.0):.4f}"
                if assay_mode == "qpcr" and probe_validation_enabled else ""
            ),
            (
                f"{(100.0 * (blast_analysis['exact_0'] + blast_analysis['mismatch_1']) / total_target_sequences if total_target_sequences else 0.0):.4f}"
                if assay_mode == "qpcr" and probe_validation_enabled else ""
            ),
            (
                f"{(100.0 * (blast_analysis['exact_0'] + blast_analysis['mismatch_1'] + blast_analysis['mismatch_2']) / total_target_sequences if total_target_sequences else 0.0):.4f}"
                if assay_mode == "qpcr" and probe_validation_enabled else ""
            ),
            (
                f"{(100.0 * blast_analysis['full_length_hitids'] / total_target_sequences if total_target_sequences else 0.0):.4f}"
                if assay_mode == "qpcr" and probe_validation_enabled else ""
            ),
            optional_portable_project_path(blast_analysis["variant_statistics_file"]),
            optional_portable_project_path(blast_analysis["variant_hitid_matrix_file"]),
        ])

    final_step = 10 if assay_mode == "qpcr" else 9
    print(f"\n{final_step}. Final assay comparison")
    print("=" * len(f"{final_step}. Final assay comparison"))
    if assay_mode == "qpcr" and probe_validation_enabled:
        print(f"Denominator: complete validation database (N={total_target_sequences})")
        print(
            f"{'Scheme':<16} {'PCR':>8} {'Exact':>9} "
            f"{'<=1 MM':>9} {'<=2 MM':>9} {'Full-site':>10}"
        )
        print("-" * 64)
        for row in coverage_manifest:
            print(
                f"{str(row[1]):<16} {float(row[7]):>7.2f}% "
                f"{float(row[28]):>8.2f}% {float(row[29]):>8.2f}% "
                f"{float(row[30]):>8.2f}% {float(row[31]):>8.2f}%"
            )
        print(
            "\nPCR = unique primer-pair-positive HitIDs. Probe columns require "
            "PCR positivity plus the indicated probe criterion."
        )
    else:
        print(f"{'Scheme':<20} {'PCR coverage':>14} {'PCR+ HitIDs':>12}")
        print("-" * 50)
        for row in coverage_manifest:
            print(
                f"{str(row[1]):<20} {float(row[7]):>13.2f}% {int(row[6]):>12d}"
            )
        if assay_mode == "qpcr":
            print("Probe columns were not calculated because probe validation was skipped.")

    # Save manifests.
    mfeprimer_runs_file = summary_results_dir / "mfeprimer_runs.tsv"
    with mfeprimer_runs_file.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["scheme", "primer_pair_fasta", "target_database", "mfeprimer_result"])
        writer.writerows(run_manifest)

    coverage_manifest_file = summary_results_dir / "coverage_results.tsv"
    with coverage_manifest_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "mode", "scheme", "potential_amplicons", "target_products",
            "self_priming_products", "valid_target_amplicons",
            "unique_positive_hitids", "primer_pair_coverage_percent",
            "probe_blast_total_valid_amplicons",
            "probe_blast_exact_0_mismatch_amplicons",
            "probe_blast_full_length_amplicons",
            "probe_blast_1_mismatch_amplicons",
            "probe_blast_2_mismatch_amplicons",
            "probe_blast_ge3_mismatch_amplicons",
            "probe_blast_0_mismatch_hitids",
            "probe_blast_1_mismatch_hitids",
            "probe_blast_2_mismatch_hitids",
            "probe_blast_ge3_mismatch_hitids",
            "probe_blast_partial_only_hitids",
            "probe_blast_no_hit_hitids",
            "probe_blast_full_length_ungapped_hitids",
            "valid_amplicons_fasta", "positive_hitids_file", "coverage_summary",
            "probe_blast_hitid_summary", "probe_blast_summary",
            "probe_blast_visual_alignment_by_hitid",
            "probe_blast_mismatch_positions", "qpcr_exact_global_percent",
            "qpcr_le1_mismatch_global_percent", "qpcr_le2_mismatch_global_percent",
            "qpcr_full_probe_site_global_percent", "probe_variant_statistics",
            "probe_variant_hitid_matrix",
        ])
        writer.writerows(coverage_manifest)

    print("\nValidation completed")
    print("====================")
    print(f"Mode    : {assay_mode.upper()}")
    print(f"Results : {compact_path(summary_results_dir.parent)}")
    print(f"Summary : {compact_path(coverage_manifest_file)}")
    print(f"Hairpin/dimer summary : {compact_path(secondary_summary_file)}")
    if secondary_details_file.is_file() and secondary_details_file.stat().st_size > 0:
        print(f"Hairpin/dimer details : {compact_path(secondary_details_file)}")
    if VALIDATION_LOG is not None:
        print(f"Log     : {compact_path(VALIDATION_LOG)}")
    if not VERBOSE:
        print("Details : saved to files; use --verbose for expanded terminal output.")

    return {
        "exclude_self_priming": exclude_self_priming,
        "probe_validation_enabled": probe_validation_enabled,
        "coverage_results": coverage_manifest_file,
        "mfeprimer_runs": mfeprimer_runs_file,
        "secondary_summary": secondary_summary_file,
        "secondary_details": secondary_details_file,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    global VERBOSE, VALIDATION_LOG

    cli = parse_cli()
    VERBOSE = bool(cli.verbose)

    try:
        (
            assay_mode,
            result_dir,
            design_manifest_path,
            design_manifest,
            project_name,
            design_run_id,
        ) = choose_design_run(cli.design_manifest)

        run_id, workdir, results_dir = resolve_validation_run_directories(
            cli,
            project_name,
        )

        assay_inputs_dir = workdir / "assay_inputs"
        database_workdir = workdir / "database"
        inputs_results_dir = results_dir / "inputs"
        mfeprimer_results_dir = results_dir / "mfeprimer"
        probe_results_dir = results_dir / "probe_blast"
        summary_results_dir = results_dir / "summary"

        for directory in (
            assay_inputs_dir,
            database_workdir,
            inputs_results_dir,
            mfeprimer_results_dir,
            probe_results_dir,
            summary_results_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        VALIDATION_LOG = results_dir / "validation.log"
        VALIDATION_LOG.write_text(
            "VarVAMP in silico validation log\n"
            f"Project: {project_name}\n"
            f"Validation run ID: {run_id}\n"
            f"Parent design run ID: {design_run_id}\n"
            f"Design manifest: {compact_path(design_manifest_path)}\n",
            encoding="utf-8",
        )

        print("╔══════════════════════════════════════════════╗")
        print("║       VarVAMP IN SILICO VALIDATION          ║")
        print("╚══════════════════════════════════════════════╝")
        print(f"Project     : {project_name}")
        print(f"Design run  : {design_run_id}")
        print(f"Valid. run  : {run_id}")
        print(f"Mode        : {assay_mode.upper()}")
        print(f"Work        : {compact_path(workdir)}")
        print(f"Results     : {compact_path(results_dir)}")
        print(
            "Flow: design manifest -> validation database -> length filter -> "
            "formatting -> indexing -> MFEprimer full QC -> coverage"
        )

        # Preserve the exact parent handoff inside this validation run.
        shutil.copy2(
            design_manifest_path,
            inputs_results_dir / "assay_design_manifest.json",
        )

        if assay_mode == "qpcr":
            oligos_file = find_qpcr_oligos(result_dir)
            schemes = prepare_varvamp_files(oligos_file, assay_inputs_dir)
        elif assay_mode in {"single", "tiled"}:
            schemes = prepare_single_or_tiled_files(
                result_dir,
                assay_inputs_dir,
                assay_mode,
            )
        else:
            raise RuntimeError(f"Unsupported VarVAMP mode: {assay_mode}")

        write_assay_manifest(
            inputs_results_dir / "assay_manifest.tsv",
            assay_mode,
            schemes,
        )
        selected_schemes = choose_schemes(schemes)

        ensure_validation_tools_available(assay_mode)
        original_database = choose_target_database(cli.validation_db)
        raw_stats = fasta_statistics(original_database)

        print("\nRaw database summary")
        print("--------------------")
        print(f"FASTA        : {compact_path(original_database)}")
        print(f"Sequences    : {raw_stats['sequences']}")
        print(f"Length       : {raw_stats['min_length']}-{raw_stats['max_length']} nt")
        print(f"Mean length  : {raw_stats['mean_length']:.2f} nt")
        print(f"Median length: {raw_stats['median_length']:.2f} nt")
        print(f"Ambiguous    : {raw_stats['ambiguous']}")

        if not ask_yes_no("Use this FASTA as the source validation database?"):
            print("Validation cancelled.")
            return 0

        filtered_database, length_filter_description, filtered_stats = (
            choose_sequence_length_filter(
                original_database,
                database_workdir,
            )
        )

        database = format_target_database_for_mfeprimer(
            filtered_database,
            database_workdir,
        )
        index_info = index_target_database(database)

        write_validation_database_summary(
            inputs_results_dir / "validation_database.tsv",
            source_database=original_database,
            source_stats=raw_stats,
            filter_description=length_filter_description,
            retained_database=filtered_database,
            retained_stats=filtered_stats,
            formatted_database=database,
        )

        parameters = choose_mfeprimer_parameters()
        min_target_size = None
        max_target_size = None
        total_target_sequences = int(filtered_stats["sequences"])

        print("\n7. Validation plan")
        print("==================")
        print(f"Project           : {project_name}")
        print(f"Parent design run : {design_run_id}")
        print(f"Validation run    : {run_id}")
        print(f"Mode              : {assay_mode.upper()}")
        print(f"Source database   : {compact_path(original_database)}")
        print(f"Length filter     : {length_filter_description}")
        print(
            f"Validation set    : {filtered_database.name} "
            f"(N={total_target_sequences})"
        )
        print(f"MFEprimer FASTA   : {compact_path(database)}")
        print(
            f"MFEprimer index   : {index_info['status']} "
            f"({index_info['detail']})"
        )
        print(f"Assays            : {', '.join(selected_schemes)}")

        if parameters is None:
            print("MFEprimer search  : defaults (0-2000 bp; Tm >=30 °C)")
        else:
            print(
                "MFEprimer search  : "
                f"{parameters['min_size']}-{parameters['max_size']} bp; "
                f"Tm >={parameters['tm_cutoff']} °C"
            )

        print("Coverage          : unique HitIDs / final retained validation database")
        if assay_mode == "qpcr":
            print("Probe             : optional, after PCR-product filtering")
            print(
                "Secondary QC      : MFEprimer hairpin/dimer sections concern "
                "LEFT/RIGHT primer inputs; probe structure is not tested there"
            )
        else:
            print("Probe             : not applicable")

        if not ask_yes_no("Run validation now?"):
            print("Validation cancelled.")
            return 0

        run_result = run_mfeprimer_pairs(
            schemes,
            selected_schemes,
            assay_inputs_dir,
            mfeprimer_results_dir,
            probe_results_dir,
            summary_results_dir,
            database,
            parameters,
            total_target_sequences,
            min_target_size,
            max_target_size,
            assay_mode,
        )

        validation_manifest = results_dir / "validation_manifest.json"
        write_validation_manifest(
            validation_manifest,
            project_name=project_name,
            run_id=run_id,
            design_manifest_path=design_manifest_path,
            design_manifest=design_manifest,
            design_run_id=design_run_id,
            assay_mode=assay_mode,
            varvamp_result_dir=result_dir,
            workdir=workdir,
            results_dir=results_dir,
            source_database=original_database,
            raw_stats=raw_stats,
            length_filter_description=length_filter_description,
            retained_database=filtered_database,
            retained_stats=filtered_stats,
            formatted_database=database,
            index_info=index_info,
            selected_schemes=selected_schemes,
            mfeprimer_parameters=parameters,
            run_result=run_result,
        )

        print(f"Validation manifest : {compact_path(validation_manifest)}")
        print(f"Permanent results   : {compact_path(results_dir)}")
        print(f"Working files       : {compact_path(workdir)}")
        return 0

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        KeyboardInterrupt,
    ) as error:
        log_line("ERROR: " + str(error))
        print(f"\nERROR: {error}", file=sys.stderr)
        return 1



if __name__ == "__main__":
    raise SystemExit(main())