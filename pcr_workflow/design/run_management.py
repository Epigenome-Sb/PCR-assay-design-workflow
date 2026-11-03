# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Execution, logging, directories and provenance for Workflow 01."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


VERBOSE = False
WORKFLOW_LOG: Path | None = None


def available_cpu_threads() -> int:
    """Return the logical CPU threads available to this workflow process."""

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
    """Append one line to the workflow log when logging is active."""
    if WORKFLOW_LOG is None:
        return
    WORKFLOW_LOG.parent.mkdir(parents=True, exist_ok=True)
    with WORKFLOW_LOG.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")


def compact_path(path: Path | str) -> str:
    """Prefer a path relative to the current working directory for display."""
    path_obj = Path(path)
    try:
        return str(path_obj.resolve().relative_to(Path.cwd().resolve()))
    except Exception:
        return str(path_obj)


def portable_manifest_path(path: Path | str) -> str:
    """
    Serialize a path for machine-readable manifests.

    Paths located inside the repository/current working directory are stored
    relative to that directory so the same manifest can be reused after moving
    the project, including between native execution and Docker (/app).

    Truly external paths remain absolute because there is no safe repository-
    relative representation for them.
    """
    resolved = Path(path).expanduser().resolve()
    repository_root = Path.cwd().resolve()

    try:
        return resolved.relative_to(repository_root).as_posix()
    except ValueError:
        return str(resolved)


def section(title: str) -> None:
    print(f"\n{title}")
    print("─" * len(title))


def progress(label: str, status: str = "OK") -> None:
    dots = "." * max(2, 52 - len(label))
    print(f"{label} {dots} {status}")


def executable_version(executable: str) -> str:
    """
    Best-effort capture of an external tool version.

    Version lookup failures never stop the scientific workflow.
    """
    candidates = (
        [executable, "--version"],
        [executable, "-version"],
        [executable, "version"],
    )

    for command in candidates:
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


def write_project_manifest(
    path: Path,
    *,
    project_name: str,
    run_id: str,
    input_file: Path,
    workdir: Path,
    results_dir: Path,
    downstream_alignment: Path,
    initial_count: int,
    after_redundancy: int,
    args: argparse.Namespace,
    varvamp_results: Path | None,
    varvamp_attempts: list[dict[str, object]] | None = None,
) -> None:
    """
    Write a machine-readable handoff between assay design and validation.

    02_in_silico_validation.py can use this file to know which VarVAMP mode
    and result directory belong to this project.
    """
    selected_tools = {
        "mafft": executable_version(args.mafft_executable),
    }

    if args.topology == "circular":
        selected_tools["mars"] = executable_version(args.mars_executable)

    if args.redundancy == "seqkit":
        selected_tools["seqkit"] = executable_version(args.seqkit_executable)
    elif args.redundancy == "cdhit":
        selected_tools["cd-hit-est"] = executable_version(args.cdhit_executable)

    if args.trimming != "none":
        selected_tools["trimal"] = executable_version(args.trimal_executable)

    if not args.skip_varvamp:
        selected_tools["varvamp"] = executable_version(args.varvamp_executable)

    # Only publish stage directories that contain output files.
    def contains_files(directory: Path) -> bool:
        return directory.is_dir() and any(
            path.is_file() for path in directory.rglob("*")
        )

    # Resolve directories from artifacts actually produced, rather than assuming
    # that a stage was executed or that the current layout is already migrated.
    alignment_source = "trimal" if args.trimming != "none" else "mafft"
    conservation_candidates = (
        results_dir / "alignment" / alignment_source / "conservation",
        results_dir / "conservation",  # Legacy Workflow 01 runs.
    )
    conservation_directory = next(
        (p for p in conservation_candidates if contains_files(p)), None
    )

    configuration_candidates = (
        results_dir / "configuration" / "varvamp",
        results_dir / "config",  # Legacy Workflow 01 runs.
    )
    configuration_directory = next(
        (p for p in configuration_candidates if contains_files(p)), None
    )

    # These entries identify the scientific stage and the actual source tool.
    # Only directories that exist in this run are recorded.
    stage_paths = {
        "input": results_dir / "input",
        "preprocessing.orientation_mafft": (
            results_dir / "preprocessing" / "orientation_mafft"
        ),
        "preprocessing.topology_mars": (
            results_dir / "preprocessing" / "topology_mars"
        ),
        "preprocessing.redundancy.seqkit": (
            results_dir / "preprocessing" / "redundancy" / "seqkit"
        ),
        "preprocessing.redundancy.cdhit_est": (
            results_dir / "preprocessing" / "redundancy" / "cdhit_est"
        ),
        "alignment.mafft": results_dir / "alignment" / "mafft",
        "alignment.mafft.quality_control": (
            results_dir / "alignment" / "mafft" / "quality_control"
        ),
        "alignment.mafft.conservation": (
            results_dir / "alignment" / "mafft" / "conservation"
        ),
        "alignment.trimal": results_dir / "alignment" / "trimal",
        "alignment.trimal.quality_control": (
            results_dir / "alignment" / "trimal" / "quality_control"
        ),
        "alignment.trimal.conservation": (
            results_dir / "alignment" / "trimal" / "conservation"
        ),
        "configuration.varvamp": results_dir / "configuration" / "varvamp",
    }
    stage_directories = {
        label: portable_manifest_path(directory)
        for label, directory in stage_paths.items()
        if contains_files(directory)
    }

    manifest = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "workflow_stage": "design",
        "path_base": "repository_root",
        "path_policy": (
            "repository-relative when inside the project; "
            "absolute only for external paths"
        ),
        "project_name": project_name,
        "run_id": run_id,
        "input_role": "design_database",
        "input_fasta": portable_manifest_path(input_file),
        "input_sha256": sha256_file(input_file),
        "workdir": portable_manifest_path(workdir),
        "results_dir": portable_manifest_path(results_dir),
        "result_directories": {
            "preprocessing": portable_manifest_path(results_dir / "preprocessing"),
            "alignment": portable_manifest_path(results_dir / "alignment"),
            "conservation": (
                portable_manifest_path(conservation_directory)
                if conservation_directory is not None else None
            ),
            "config": (
                portable_manifest_path(configuration_directory)
                if configuration_directory is not None else None
            ),
            "varvamp": (
                None
                if varvamp_results is None
                else portable_manifest_path(varvamp_results)
            ),
        },
        "stage_directories": stage_directories,
        "selected_alignment_tool": alignment_source,
        "conservation_input_alignment": portable_manifest_path(downstream_alignment),
        "downstream_alignment": portable_manifest_path(downstream_alignment),
        "sequence_counts": {
            "initial_input": int(initial_count),
            "entering_final_mafft": int(after_redundancy),
        },
        "preprocessing": {
            "orientation": args.orientation,
            "topology": args.topology,
            "redundancy": args.redundancy,
            "cdhit_identity": args.identity,
            "cdhit_mode": args.cdhit_mode,
            "cdhit_strand": args.cdhit_strand,
            "mafft_strategy": args.mafft_strategy,
            "trimming": args.trimming,
            "trimal_gap_threshold": args.trimal_gap_threshold,
            "min_occupancy": args.min_occupancy,
            "min_major_frequency": args.min_major_frequency,
            "threads": int(args.threads),
        },
        "varvamp": {
            "executed": bool(varvamp_attempts)
            if varvamp_attempts is not None
            else not args.skip_varvamp,
            "mode": None if args.skip_varvamp else args.varvamp_mode,
            "result_dir": (
                None
                if varvamp_results is None
                else portable_manifest_path(varvamp_results)
            ),
            "consensus_threshold": args.varvamp_threshold,
            "primer_ambiguity": args.primer_ambiguity,
            "probe_ambiguity": args.probe_ambiguity,
            "single_opt_length": args.single_opt_length,
            "single_max_length": args.single_max_length,
            "single_report_n": args.single_report_n,
            "tiled_opt_length": args.tiled_opt_length,
            "tiled_max_length": args.tiled_max_length,
            "tiled_overlap": args.tiled_overlap,
            "qpcr_test_n": args.qpcr_test_n,
            "qpcr_deltag": args.qpcr_deltag,
            "custom_config": (
                None
                if args.varvamp_config is None
                else portable_manifest_path(args.varvamp_config)
            ),
            "attempts": varvamp_attempts or [],
        },
        "tool_versions": selected_tools,
        "workflow_log": (
            None
            if WORKFLOW_LOG is None
            else portable_manifest_path(WORKFLOW_LOG)
        ),
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )


def write_workflow_summary(
    path: Path,
    *,
    project_name: str,
    run_id: str,
    input_file: Path,
    downstream_alignment: Path,
    workdir: Path,
    results_dir: Path,
    initial_count: int,
    after_redundancy: int,
    args: argparse.Namespace,
    varvamp_results: Path | None,
) -> None:
    lines = [
        "VarVAMP ASSAY-DESIGN WORKFLOW SUMMARY",
        "=" * 55,
        f"Project: {project_name}",
        f"Run ID: {run_id}",
        "Workflow stage: design",
        f"Design FASTA: {compact_path(input_file)}",
        f"Initial sequences: {initial_count}",
        f"Sequences entering final MAFFT: {after_redundancy}",
        f"Orientation: {args.orientation}",
        f"Topology: {args.topology}",
        f"Redundancy: {args.redundancy}",
        f"Final MAFFT strategy: {args.mafft_strategy}",
        f"Threads: {args.threads}",
        f"Trimming: {args.trimming}",
        f"Downstream alignment: {compact_path(downstream_alignment)}",
        "",
        "VarVAMP:",
        f"Executed: {'yes' if not args.skip_varvamp else 'no'}",
        f"Mode: {args.varvamp_mode if not args.skip_varvamp else 'N/A'}",
        (
            f"Result directory: {compact_path(varvamp_results)}"
            if varvamp_results is not None
            else "Result directory: N/A"
        ),
        "",
        f"Work directory: {compact_path(workdir)}",
        f"Results directory: {compact_path(results_dir)}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sanitize_project_name(value: str) -> str:
    """Convert a project name into a safe file/directory identifier."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("The project name is empty or invalid.")
    return cleaned


def sanitize_run_id(value: str) -> str:
    """Validate/sanitize a run identifier used as a directory name."""

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("The run ID is empty or invalid.")
    return cleaned


def generate_run_id() -> str:
    """Generate a human-readable local timestamp for one workflow execution."""

    return datetime.now().strftime("%Y%m%d_%H%M%S")


def resolve_run_directories(
    args: argparse.Namespace,
    project_name: str,
) -> tuple[str, Path, Path]:
    """
    Resolve isolated work/results directories for one assay-design run.

    The default layout is:
      work/<project>/design/<run_id>/
      results/<project>/design/<run_id>/

    When the automatically generated timestamp already exists, a numeric
    suffix is added so an earlier run is never overwritten. Explicit
    --workdir/--results values remain exact per-run overrides.
    """

    requested_run_id = (
        sanitize_run_id(args.run_id)
        if args.run_id is not None
        else generate_run_id()
    )

    if args.workdir is not None or args.results is not None:
        workdir = (
            args.workdir.expanduser().resolve()
            if args.workdir is not None
            else (Path("work") / project_name / "design" / requested_run_id).resolve()
        )
        results_dir = (
            args.results.expanduser().resolve()
            if args.results is not None
            else (Path("results") / project_name / "design" / requested_run_id).resolve()
        )
        return requested_run_id, workdir, results_dir

    candidate = requested_run_id
    suffix = 1
    while True:
        workdir = (Path("work") / project_name / "design" / candidate).resolve()
        results_dir = (Path("results") / project_name / "design" / candidate).resolve()
        if not workdir.exists() and not results_dir.exists():
            return candidate, workdir, results_dir
        suffix += 1
        candidate = f"{requested_run_id}_{suffix:02d}"


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file for provenance tracking."""

    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_executable(executable: str) -> str | None:
    """Return a usable executable path or None."""

    path = Path(executable).expanduser()
    if path.parent != Path(".") or path.is_absolute():
        return str(path.resolve()) if path.is_file() else None
    return shutil.which(executable)


def check_required_tools(tools: Iterable[tuple[str, str]]) -> None:
    """Check that selected external tools are available."""

    missing: list[str] = []
    for label, executable in tools:
        if resolve_executable(executable) is None:
            missing.append(f"{label} ({executable})")

    if missing:
        raise RuntimeError(
            "Required program(s) not found: " + ", ".join(missing) + "."
        )


def run_command(
    command: list[str],
    *,
    stdout_file: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    """Run an external command, log it, and stop on failure."""
    command_text = " ".join(command)
    if stdout_file is not None:
        command_text += f" > {stdout_file}"

    log_line("COMMAND: " + command_text)

    if VERBOSE:
        print("\nCommand:", command_text)

    try:
        if stdout_file is None:
            subprocess.run(command, check=True, env=env)
        else:
            stdout_file.parent.mkdir(parents=True, exist_ok=True)
            with stdout_file.open("w", encoding="utf-8") as handle:
                subprocess.run(
                    command,
                    check=True,
                    stdout=handle,
                    env=env,
                )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Command failed with exit code {error.returncode}: "
            + command_text
        ) from error
