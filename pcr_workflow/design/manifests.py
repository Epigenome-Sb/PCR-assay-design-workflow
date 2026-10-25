# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import json

from ..core import runtime
from ..core.paths import compact_path, portable_manifest_path, sha256_file
from ..core.runtime import executable_version

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
            "conservation": portable_manifest_path(results_dir / "conservation"),
            "config": portable_manifest_path(results_dir / "config"),
            "varvamp": (
                None
                if varvamp_results is None
                else portable_manifest_path(varvamp_results)
            ),
        },
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
            if runtime.LOG_PATH is None
            else portable_manifest_path(runtime.LOG_PATH)
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
