# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import json

from ..core import runtime
from ..core.paths import generate_run_id, portable_project_path, sanitize_identifier, sha256_file
from ..core.runtime import executable_version

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
    threads: int,
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
        "workflow_log": None if runtime.LOG_PATH is None else portable_project_path(runtime.LOG_PATH),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
