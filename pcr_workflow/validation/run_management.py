# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Workflow 02 execution, logging, paths, and reproducibility outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

VERBOSE = False
VALIDATION_LOG: Path | None = None


def available_cpu_threads() -> int:
    """Return logical CPU threads available to this workflow process."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def automatic_thread_count(cpu_threads: int) -> int:
    """Keep two logical CPU threads free when possible."""
    if cpu_threads <= 2:
        return 1
    return max(1, cpu_threads - 2)


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


def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))


def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 48 - len(label))
    print(f"{label} {dots} {status}")


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
        for directory in (workdir, results_dir):
            if directory.exists() and any(directory.iterdir()):
                raise RuntimeError(
                    f"Run directory is not empty; refusing to overwrite: {directory}"
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


def _has_files(directory: Path) -> bool:
    """Report only result directories containing actual output files."""
    return directory.is_dir() and any(path.is_file() for path in directory.rglob("*"))


def _permanent_retained_fasta(
    results_dir: Path,
    source_database: Path,
    retained_database: Path,
) -> Path | None:
    """Find the retained FASTA copy without claiming one was produced."""
    copy = results_dir / "database" / "filtered" / retained_database.name
    if copy.is_file():
        return copy
    if retained_database.resolve() == source_database.resolve():
        return source_database
    return None


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
        # Preserve the existing TSV fields and append provenance only when a
        # permanent FASTA is really available (filtered copy or source FASTA).
        permanent = _permanent_retained_fasta(
            path.parent.parent, source_database, retained_database
        )
        if permanent is not None:
            rows.append(("retained_fasta_permanent", portable_project_path(permanent)))
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
    threads: int,
    run_result: dict[str, object],
) -> None:
    """Write a machine-readable provenance record for Workflow 02."""
    index_files = [
        portable_project_path(Path(path))
        for path in index_info.get("files", [])
    ]
    permanent_retained = _permanent_retained_fasta(
        results_dir, source_database, retained_database
    )

    # Existing keys remain available to older consumers; inactive stages are
    # represented as null rather than pointing to non-existent directories.
    result_paths = {
        "inputs": results_dir / "inputs",
        "database": results_dir / "database",
        "mfeprimer": results_dir / "mfeprimer",
        "probe_blast": results_dir / "probe_blast",
        "coverage": results_dir / "coverage",
        "summary": results_dir / "summary",  # Legacy layout until reporting.py migrates.
    }
    result_directories = {
        key: portable_project_path(directory) if _has_files(directory) else None
        for key, directory in result_paths.items()
    }

    stage_paths = {
        "inputs": results_dir / "inputs",
        "inputs.assay_sequences": results_dir / "inputs" / "assay_sequences",
        "database": results_dir / "database",
        "database.filtered": results_dir / "database" / "filtered",
        "mfeprimer": results_dir / "mfeprimer",
        "mfeprimer.quality_control": results_dir / "mfeprimer" / "quality_control",
        "probe_blast": results_dir / "probe_blast",
        "coverage": results_dir / "coverage",
        "summary": results_dir / "summary",
    }
    for scheme in selected_schemes:
        stage_paths[f"mfeprimer.{scheme}"] = results_dir / "mfeprimer" / scheme
        stage_paths[f"probe_blast.{scheme}"] = results_dir / "probe_blast" / scheme
        stage_paths[f"coverage.{scheme}"] = results_dir / "coverage" / scheme
    stage_directories = {
        key: portable_project_path(directory)
        for key, directory in stage_paths.items()
        if _has_files(directory)
    }

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
        "result_directories": result_directories,
        "stage_directories": stage_directories,
        "selected_assays": selected_schemes,
        "threads": {
            "supported_tools": int(threads),
            "blastn": int(threads) if assay_mode == "qpcr" else "not applicable",
            "seqkit_formatting": min(int(threads), 4),
            "mfeprimer_full_qc": "no user CPU flag in MFEprimer 3.1 full-QC command",
        },
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
    if permanent_retained is not None:
        manifest["validation_database"]["retained_fasta_permanent"] = (
            portable_project_path(permanent_retained)
        )
    permanent_assay_inputs = results_dir / "inputs" / "assay_sequences"
    if _has_files(permanent_assay_inputs):
        manifest["assay_sequences_dir"] = portable_project_path(permanent_assay_inputs)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
