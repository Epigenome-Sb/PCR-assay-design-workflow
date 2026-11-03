#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""
pcr_workflow.design.workflow

Generic stepwise VarVAMP assay-design workflow for nucleotide sequence datasets.

Scientific workflow
-------------------
1. Validate the input nucleotide FASTA.
2. Optionally normalize sequence orientation with MAFFT.
3. For complete circular sequences, optionally normalize starts with MARS.
4. Handle redundancy with no filtering, SeqKit exact deduplication, or CD-HIT-EST.
5. Perform the final MAFFT multiple-sequence alignment.
6. Calculate alignment QC statistics and optionally trim with trimAl.
7. Analyse alignment conservation.
8. Design SINGLE, QPCR or TILED assays with VarVAMP.
9. Save VarVAMP outputs and visual summaries.
10. Write assay_design_manifest.json for reproducible handoff to
    02_in_silico_validation.py.

Interface
---------
Scientific decisions remain stepwise and explicit in interactive mode.

Normal terminal mode is concise.
Use --verbose to display complete external commands and detailed tables.
All commands are recorded in the run-specific workflow log.

Default storage architecture
----------------------------
Input design FASTA:
    data/design/

Temporary/reconstructible files:
    work/<project>/design/<run_id>/

Permanent scientific results:
    results/<project>/design/<run_id>/

Each run is isolated by a unique run_id so previous analyses are not
overwritten. The run contains preprocessing, nested alignment QC/conservation,
VarVAMP, configuration, manifest, summary and log outputs.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import pandas as pd

from . import run_management as _run_management
from .alignment import (
    MAFFT_STRATEGIES,
    calculate_alignment_statistics,
    choose_mafft_strategy,
    choose_trimming_strategy,
    print_alignment_statistics,
    read_alignment,
    run_final_mafft,
    run_trimal,
    save_alignment_statistics,
    trimming_recommendation,
    warn_about_mafft_strategy,
)
from .conservation import (
    perform_conservation_analysis,
)
from .fasta import (
    calculate_fasta_statistics,
    count_fasta_sequences,
    print_fasta_statistics,
    validate_nucleotide_fasta,
)
from .questions import (
    ask_required_choice,
    ask_required_float,
    ask_required_non_negative_integer,
    ask_yes_no,
    validate_fraction,
)
from .redundancy import (
    choose_redundancy_strategy,
    run_cdhit_clustering,
    run_seqkit_deduplication,
)
from .reporting import (
    create_varvamp_visual_summaries,
)
from .run_management import (
    available_cpu_threads,
    automatic_thread_count,
    check_required_tools,
    compact_path,
    log_line,
    portable_manifest_path,
    resolve_run_directories,
    sanitize_project_name,
    section,
    write_project_manifest,
    write_workflow_summary,
)
from .sequence_orientation import (
    calculate_cyclic_rotation_statistics,
    choose_orientation_strategy,
    choose_sequence_topology,
    normalize_orientation_with_mafft,
    print_rotation_statistics,
    run_mars_rotation,
)
from .varvamp import (
    adjust_varvamp_after_failure,
    copy_varvamp_tree,
    ensure_mode_required_interactive_parameters,
    parse_varvamp_diagnostics,
    print_varvamp_failure,
    print_varvamp_tables,
    run_varvamp_command,
)
from .varvamp_config import (
    build_varvamp_command,
    choose_varvamp_mode,
    configure_interactive_varvamp_advanced,
    stage_varvamp_config_for_run,
    validate_varvamp_mode_parameters,
)


def format_trimal_effect(
    pretrim_stats: dict[str, float | int],
    posttrim_stats: dict[str, float | int],
    strategy: str,
    threshold: float | None,
) -> str:
    """Describe the observed before/after change without interpreting biology."""
    before_columns = int(pretrim_stats["alignment_length"])
    after_columns = int(posttrim_stats["alignment_length"])
    removed_columns = before_columns - after_columns
    removed_percent = 100.0 * removed_columns / before_columns
    before_gaps = float(pretrim_stats["mean_gap_content_percent"])
    after_gaps = float(posttrim_stats["mean_gap_content_percent"])
    before_occupied = int(pretrim_stats["fully_occupied_columns"])
    after_occupied = int(posttrim_stats["fully_occupied_columns"])
    before_occupancy = float(pretrim_stats["mean_sequence_occupancy_percent"])
    after_occupancy = float(posttrim_stats["mean_sequence_occupancy_percent"])
    final_status, _ = trimming_recommendation(posttrim_stats)

    lines = [
        "trimAl effect",
        "============",
        f"Strategy                  : {strategy}",
    ]
    if strategy == "manual":
        lines.append(f"Gap threshold (-gt)       : {threshold}")
    lines.extend(
        [
            f"Columns before            : {before_columns}",
            f"Columns after             : {after_columns}",
            f"Columns removed           : {removed_columns} ({removed_percent:.2f} %)",
            f"Mean gap content before   : {before_gaps:.2f} %",
            f"Mean gap content after    : {after_gaps:.2f} %",
            "",
            "trimAl comparison summary",
            "=========================",
            f"{'Metric':<29} {'Before':>10} {'After':>10}",
            "-" * 51,
            f"{'Alignment length':<29} {before_columns:>10} {after_columns:>10}",
            f"{'Mean gap content (%)':<29} {before_gaps:>10.2f} {after_gaps:>10.2f}",
            f"{'Fully occupied columns':<29} {before_occupied:>10} {after_occupied:>10}",
            f"{'Mean sequence occupancy (%)':<29} {before_occupancy:>10.2f} {after_occupancy:>10.2f}",
            "",
            "Trimming assessment",
            "===================",
            "trimAl completed; the QC statistics below describe its effect.",
            f"Gap content change        : {after_gaps - before_gaps:+.2f} percentage points",
            f"Final QC classification   : {final_status}",
            "Note: These statistics do not establish biological alignment validity.",
        ]
    )
    return "\n".join(lines) + "\n"


def record_design_artifacts(
    manifest_file: Path,
    *,
    project_name: str,
    input_stats_file: Path,
    results_dir: Path,
    orientation_file: Path | None,
    mars_file: Path | None,
    redundancy_file: Path | None,
    mafft_alignment_file: Path,
    mafft_qc_dir: Path,
    trimal_alignment_file: Path | None,
    trimal_qc_dir: Path,
    trimal_effect_file: Path | None,
    conservation_dir: Path,
    selected_alignment_tool: str,
    inventory_file: Path,
) -> None:
    """Extend, never replace, the existing manifest schema used by validation."""
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    directories = manifest.setdefault("result_directories", {})
    directories["conservation"] = portable_manifest_path(conservation_dir)
    directories["alignment_mafft"] = portable_manifest_path(mafft_alignment_file.parent)
    directories["alignment_trimal"] = (
        portable_manifest_path(trimal_alignment_file.parent)
        if trimal_alignment_file is not None else None
    )

    redundancy_report = (
        redundancy_file.parent
        / f"{project_name}_{redundancy_file.parent.name}_redundancy_report.txt"
        if redundancy_file is not None else None
    )
    files = {
        "input_statistics": input_stats_file,
        "orientation_fasta": orientation_file,
        "orientation_report": (
            orientation_file.parent / f"{project_name}_mafft_orientation_report.txt"
            if orientation_file is not None else None
        ),
        "mars_fasta": mars_file,
        "mars_report": (
            mars_file.parent / f"{project_name}_mars_rotation_report.txt"
            if mars_file is not None else None
        ),
        "redundancy_fasta": redundancy_file,
        "redundancy_report": redundancy_report,
        "redundancy_clusters": (
            Path(str(redundancy_file) + ".clstr")
            if redundancy_file is not None and redundancy_file.parent.name == "cdhit_est"
            else None
        ),
        "mafft_alignment": mafft_alignment_file,
        "mafft_quality_control": mafft_qc_dir / f"{project_name}_mafft_alignment_qc.txt",
        "mafft_sequence_occupancy": mafft_qc_dir / f"{project_name}_mafft_sequence_occupancy.csv",
        "trimal_alignment": trimal_alignment_file,
        "trimal_effect": trimal_effect_file,
        "trimal_quality_control": (
            trimal_qc_dir / f"{project_name}_trimal_alignment_qc.txt"
            if trimal_alignment_file is not None else None
        ),
        "trimal_sequence_occupancy": (
            trimal_qc_dir / f"{project_name}_trimal_sequence_occupancy.csv"
            if trimal_alignment_file is not None else None
        ),
        "conservation_summary": conservation_dir / f"{project_name}_conservation_summary.txt",
        "conservation_by_position": conservation_dir / f"{project_name}_conservation_by_position.csv",
        "conservation_consensus": conservation_dir / f"{project_name}_consensus.fasta",
        "output_inventory": inventory_file,
    }
    manifest["selected_alignment_tool"] = selected_alignment_tool
    manifest["design_artifacts"] = {
        key: portable_manifest_path(file_path)
        if file_path is not None and (file_path.is_file() or file_path == inventory_file)
        else None
        for key, file_path in files.items()
    }
    manifest_file.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def write_output_inventory(results_dir: Path, inventory_file: Path) -> None:
    """List permanent outputs with sizes, excluding the inventory itself."""
    entries = sorted(
        path for path in results_dir.rglob("*")
        if path.is_file() and path != inventory_file
    )
    inventory_file.write_text(
        "relative_path\tsize_bytes\n" + "".join(
            f"{path.relative_to(results_dir).as_posix()}\t{path.stat().st_size}\n"
            for path in entries
        ),
        encoding="utf-8",
    )

# The CLI and workflow order are preserved from the original script.

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Nucleotide preprocessing, MAFFT alignment, alignment QC, "
            "conservation analysis and PCR/qPCR/tiled assay design with VarVAMP."
        )
    )

    parser.add_argument("--input", type=Path, default=None, help="Input FASTA file.")
    parser.add_argument("--project-name", type=str, default=None)
    parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help=(
            "Optional run identifier. By default a timestamp such as "
            "20260906_164500 is generated automatically."
        ),
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help=(
            "Optional override for this run work directory. Default: "
            "work/<project>/design/<run_id>."
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=None,
        help=(
            "Optional override for this run results directory. Default: "
            "results/<project>/design/<run_id>."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "Number of CPU threads used by supported tools. If omitted, the "
            "workflow detects the logical CPU threads available to the process "
            "and keeps two threads free when possible."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Display complete external commands and detailed file paths. "
            "Commands are always saved to workflow.log."
        ),
    )

    parser.add_argument(
        "--orientation",
        choices=("keep", "adjust", "accurate"),
        default=None,
        help=(
            "Sequence-orientation strategy: keep, MAFFT --adjustdirection, or "
            "MAFFT --adjustdirectionaccurately."
        ),
    )
    parser.add_argument(
        "--topology",
        choices=("linear", "circular"),
        default=None,
        help="Sequence topology. MARS is run only for circular datasets.",
    )

    parser.add_argument(
        "--redundancy",
        choices=("none", "seqkit", "cdhit"),
        default=None,
        help=(
            "Redundancy strategy: keep all sequences, remove exact duplicates "
            "with SeqKit, or cluster with CD-HIT-EST."
        ),
    )
    # Backward-compatible flag from the previous workflow.
    parser.add_argument(
        "--skip-cdhit",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--identity", type=float, default=None)
    parser.add_argument(
        "--cdhit-mode", choices=("fast", "accurate"), default=None
    )
    parser.add_argument(
        "--cdhit-strand", choices=("both", "same"), default=None
    )

    parser.add_argument(
        "--mafft-strategy",
        choices=tuple(MAFFT_STRATEGIES),
        default=None,
        help="Final MAFFT alignment strategy.",
    )

    parser.add_argument(
        "--trimming",
        choices=("none", "noallgaps", "gappyout", "automated1", "manual"),
        default=None,
        help="Alignment trimming strategy after final MAFFT alignment.",
    )
    parser.add_argument(
        "--trimal-gap-threshold",
        type=float,
        default=None,
        help=(
            "trimAl -gt value for --trimming manual. This is the minimum "
            "fraction of sequences without a gap in a retained column."
        ),
    )

    parser.add_argument(
        "--min-occupancy",
        type=float,
        default=0.95,
        help="Minimum occupancy for conservation analysis (default: 0.95).",
    )
    parser.add_argument(
        "--min-major-frequency",
        type=float,
        default=0.95,
        help="Minimum major-base frequency for conservation (default: 0.95).",
    )

    parser.add_argument(
        "--varvamp-mode",
        choices=("single", "qpcr", "tiled"),
        default=None,
        help=(
            "VarVAMP mode: single (conventional PCR), qpcr, or tiled "
            "(full-genome amplicon sequencing)."
        ),
    )
    parser.add_argument("--varvamp-threshold", type=float, default=None)
    parser.add_argument("--primer-ambiguity", type=int, default=None)
    parser.add_argument("--probe-ambiguity", type=int, default=None)
    parser.add_argument("--skip-varvamp", action="store_true")

    # Optional mode-specific VarVAMP CLI parameters. When omitted, VarVAMP's
    # own documented defaults are left untouched.
    parser.add_argument("--single-opt-length", type=int, default=None)
    parser.add_argument("--single-max-length", type=int, default=None)
    parser.add_argument(
        "--single-report-n",
        type=str,
        default=None,
        help="VarVAMP single -n value; integer or 'inf'.",
    )
    parser.add_argument("--tiled-opt-length", type=int, default=None)
    parser.add_argument("--tiled-max-length", type=int, default=None)
    parser.add_argument("--tiled-overlap", type=int, default=None)
    parser.add_argument("--qpcr-test-n", type=int, default=None)
    parser.add_argument("--qpcr-deltag", type=int, default=None)
    parser.add_argument(
        "--varvamp-config",
        type=Path,
        default=None,
        help="Optional custom VarVAMP config file passed through VARVAMP_CONFIG.",
    )

    # Executable names/paths make the script easier to use in custom envs.
    parser.add_argument("--mafft-executable", default="mafft")
    parser.add_argument("--mars-executable", default="mars")
    parser.add_argument("--seqkit-executable", default="seqkit")
    parser.add_argument("--cdhit-executable", default="cd-hit-est")
    parser.add_argument("--trimal-executable", default="trimal")
    parser.add_argument("--varvamp-executable", default="varvamp")

    return parser.parse_args()



def ask_for_input_file() -> Path:
    """Interactively select the Workflow 01 design FASTA."""

    design_dir = Path("data") / "design"
    candidates: list[Path] = []
    if design_dir.is_dir():
        for pattern in ("*.fasta", "*.fa", "*.fna", "*.fas"):
            candidates.extend(design_dir.glob(pattern))
        candidates = sorted({path.resolve() for path in candidates if path.is_file()})

    print("\nDesign database")
    print("===============")
    print("Workflow 01 uses the FASTA stored in data/design/ for assay design.")

    if candidates:
        for index, path in enumerate(candidates, start=1):
            print(f"{index}. {compact_path(path)}")
        print(f"{len(candidates) + 1}. Enter another FASTA path")

        while True:
            try:
                raw = input("Select the design FASTA by number: ").strip()
            except EOFError as error:
                raise RuntimeError(
                    "No interactive input is available. Use --input."
                ) from error

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
        try:
            raw_value = input("Name or path of the design FASTA file: ").strip()
        except EOFError as error:
            raise RuntimeError(
                "No interactive input is available. Use --input."
            ) from error

        raw_value = raw_value.strip("'\"")
        if not raw_value:
            print("Please enter a filename or file path.")
            continue

        candidate = Path(raw_value).expanduser()
        possible_paths = [candidate]
        if not candidate.is_absolute() and candidate.parent == Path("."):
            possible_paths.append(design_dir / candidate)
            # Compatibility fallback for older repository layouts.
            possible_paths.append(Path("data") / candidate)

        for possible_path in possible_paths:
            resolved = possible_path.resolve()
            if resolved.is_file():
                return resolved

        print(
            "File not found. Checked paths: "
            + ", ".join(str(path.resolve()) for path in possible_paths)
        )



# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------


def main() -> int:

    args = parse_arguments()
    _run_management.VERBOSE = bool(args.verbose)

    try:
        cpu_threads = available_cpu_threads()

        if args.threads is None:
            args.threads = automatic_thread_count(cpu_threads)
            thread_selection = "automatic"
        else:
            thread_selection = "manual"

        if args.threads <= 0:
            raise ValueError("--threads must be greater than zero.")

        if args.threads > cpu_threads:
            print(
                f"Warning: --threads={args.threads} exceeds the "
                f"{cpu_threads} logical CPU threads currently available."
            )

        validate_fraction(args.min_occupancy, "--min-occupancy")
        validate_fraction(args.min_major_frequency, "--min-major-frequency")

        if args.input is None:
            input_file = ask_for_input_file()
        else:
            input_file = args.input.expanduser().resolve()

        validate_nucleotide_fasta(input_file)

        project_name = sanitize_project_name(
            args.project_name if args.project_name else input_file.stem
        )
        run_id, workdir, results_dir = resolve_run_directories(
            args,
            project_name,
        )

        # Directories name the producing step/tool; create conditional ones lazily.
        orientation_workdir = workdir / "preprocessing" / "orientation_mafft"
        mars_workdir = workdir / "preprocessing" / "topology_mars"
        mafft_workdir = workdir / "alignment" / "mafft"
        trimal_workdir = workdir / "alignment" / "trimal"
        varvamp_work_root = workdir / "varvamp"

        preprocessing_results_dir = results_dir / "preprocessing"
        orientation_results_dir = preprocessing_results_dir / "orientation_mafft"
        mars_results_dir = preprocessing_results_dir / "topology_mars"
        alignment_results_dir = results_dir / "alignment"
        mafft_results_dir = alignment_results_dir / "mafft"
        mafft_qc_results_dir = mafft_results_dir / "quality_control"
        trimal_results_dir = alignment_results_dir / "trimal"
        trimal_qc_results_dir = trimal_results_dir / "quality_control"
        config_results_dir = results_dir / "configuration" / "varvamp"

        results_dir.mkdir(parents=True, exist_ok=True)

        _run_management.WORKFLOW_LOG = results_dir / "workflow.log"
        _run_management.WORKFLOW_LOG.write_text(
            "VarVAMP assay-design workflow log\n"
            f"Project: {project_name}\n"
            f"Run ID: {run_id}\n"
            f"Design FASTA: {compact_path(input_file)}\n"
            f"CPU threads available: {cpu_threads}\n"
            f"Workflow threads: {args.threads} ({thread_selection})\n",
            encoding="utf-8",
        )

        print("╔══════════════════════════════════════════════╗")
        print("║       VarVAMP ASSAY DESIGN WORKFLOW         ║")
        print("╚══════════════════════════════════════════════╝")
        print(f"Project : {project_name}")
        print(f"Run ID  : {run_id}")
        print(f"Input   : {compact_path(input_file)}")
        print(f"Work    : {compact_path(workdir)}")
        print(f"Results : {compact_path(results_dir)}")
        print(f"CPU     : {cpu_threads} logical threads available")
        print(
            f"Threads : {args.threads} used by supported tools "
            f"({thread_selection})"
        )
        if _run_management.VERBOSE:
            print(f"Log     : {_run_management.WORKFLOW_LOG}")

        interactive = sys.stdin.isatty()

        if not interactive and args.trimming is None:
            raise RuntimeError(
                "--trimming is required in non-interactive mode."
            )

        if (
            args.trimming == "manual"
            and args.trimal_gap_threshold is None
            and not interactive
        ):
            raise RuntimeError(
                "--trimal-gap-threshold is required for --trimming manual."
            )

        # ------------------------------------------------------------------
        # Input overview
        # ------------------------------------------------------------------
        current_fasta = input_file
        initial_count = count_fasta_sequences(input_file)

        # Persist exactly what the existing FASTA summary displays.
        input_stats_buffer = StringIO()
        with redirect_stdout(input_stats_buffer):
            print_fasta_statistics(
                calculate_fasta_statistics(input_file),
                title="Input dataset statistics",
            )
        print(input_stats_buffer.getvalue(), end="")
        input_stats_file = results_dir / "input" / f"{project_name}_input_statistics.txt"
        input_stats_file.parent.mkdir(parents=True, exist_ok=True)
        input_stats_file.write_text(input_stats_buffer.getvalue(), encoding="utf-8")

        # ------------------------------------------------------------------
        # 1. Orientation normalization
        # ------------------------------------------------------------------
        print("\n1. Sequence orientation checkpoint")

        if args.orientation is None:
            if not interactive:
                raise RuntimeError(
                    "--orientation is required in non-interactive mode."
                )
            args.orientation = choose_orientation_strategy()

        if args.orientation not in {"keep", "adjust", "accurate"}:
            raise ValueError("Unsupported orientation strategy.")

        reversed_count = 0
        orientation_result: Path | None = None

        if args.orientation == "keep":
            section("Orientation result")
            print("Method                    : keep original orientation")
            print("MAFFT orientation check   : skipped")
            print(f"Sequences retained        : {initial_count}")
            print("Sequences changed         : 0")
        else:
            check_required_tools([("MAFFT", args.mafft_executable)])
            print("Running sequence orientation normalization with MAFFT.")
            orientation_workdir.mkdir(parents=True, exist_ok=True)
            current_fasta, reversed_count = normalize_orientation_with_mafft(
                current_fasta,
                orientation_workdir,
                project_name,
                args.orientation,
                args.threads,
                args.mafft_executable,
            )

            section("Orientation result")
            print(
                "Method                    : "
                + (
                    "MAFFT --adjustdirection"
                    if args.orientation == "adjust"
                    else "MAFFT --adjustdirectionaccurately"
                )
            )
            print(f"Sequences checked         : {initial_count}")
            print(
                f"Reverse-complemented      : {reversed_count} "
                f"({100.0 * reversed_count / initial_count if initial_count else 0.0:.2f} %)"
            )
            print(
                f"Orientation unchanged     : {initial_count - reversed_count}"
            )

        if args.orientation != "keep":
            orientation_results_dir.mkdir(parents=True, exist_ok=True)
            orientation_result = (
                orientation_results_dir / f"{project_name}_mafft_oriented.fasta"
            )
            shutil.copy2(current_fasta, orientation_result)
            orientation_report = (
                orientation_workdir / f"{project_name}_orientation_report.txt"
            )
            if orientation_report.is_file():
                shutil.copy2(
                    orientation_report,
                    orientation_results_dir / f"{project_name}_mafft_orientation_report.txt",
                )

        # ------------------------------------------------------------------
        # 2. Circular start-position normalization
        # ------------------------------------------------------------------
        print("\n2. Sequence topology checkpoint")

        if args.topology is None:
            if not interactive:
                raise RuntimeError(
                    "--topology is required in non-interactive mode."
                )
            args.topology = choose_sequence_topology()

        mars_result: Path | None = None
        if args.topology == "circular":
            check_required_tools([("MARS", args.mars_executable)])
            print("Running circular sequence rotation with MARS.")
            before_mars = current_fasta
            mars_workdir.mkdir(parents=True, exist_ok=True)
            current_fasta = run_mars_rotation(
                current_fasta,
                mars_workdir,
                project_name,
                args.mars_executable,
                args.threads,
            )
            mars_stats = calculate_cyclic_rotation_statistics(before_mars, current_fasta)
            mars_buffer = StringIO()
            with redirect_stdout(mars_buffer):
                print_rotation_statistics(mars_stats)
            print(mars_buffer.getvalue(), end="")
            mars_results_dir.mkdir(parents=True, exist_ok=True)
            (mars_results_dir / f"{project_name}_mars_rotation_report.txt").write_text(
                mars_buffer.getvalue(), encoding="utf-8"
            )
        elif args.topology == "linear":
            section("Topology result")
            print("Topology                  : linear")
            print("MARS                      : not applicable")
            print("Cyclic starts changed     : 0")
        else:
            raise ValueError("--topology must be 'linear' or 'circular'.")

        if args.topology == "circular":
            mars_result = mars_results_dir / f"{project_name}_mars_normalized.fasta"
            shutil.copy2(current_fasta, mars_result)

        # ------------------------------------------------------------------
        # 3. Redundancy handling
        # ------------------------------------------------------------------
        print("\n3. Sequence redundancy checkpoint")

        if args.redundancy is None:
            if args.skip_cdhit:
                args.redundancy = "none"
            elif not interactive:
                raise RuntimeError(
                    "--redundancy is required in non-interactive mode."
                )
            else:
                args.redundancy = choose_redundancy_strategy()

        if args.redundancy == "cdhit":
            if args.identity is None:
                if not interactive:
                    raise RuntimeError(
                        "--identity is required for CD-HIT-EST."
                    )
                args.identity = ask_required_float(
                    "CD-HIT-EST identity threshold",
                    "-c",
                    minimum=0.80,
                    maximum=1.0,
                )

            if not 0.80 <= args.identity <= 1.0:
                raise ValueError(
                    "--identity must be between 0.80 and 1.0."
                )

            if args.cdhit_mode is None:
                if not interactive:
                    raise RuntimeError(
                        "--cdhit-mode is required for CD-HIT-EST."
                    )
                args.cdhit_mode = ask_required_choice(
                    "CD-HIT-EST clustering mode",
                    [
                        (
                            "fast",
                            "Fast mode (-g 0): assign to the first qualifying cluster.",
                        ),
                        (
                            "accurate",
                            "Accurate mode (-g 1): select the most similar qualifying representative.",
                        ),
                    ],
                )

            if args.cdhit_strand is None:
                if not interactive:
                    raise RuntimeError(
                        "--cdhit-strand is required for CD-HIT-EST."
                    )
                args.cdhit_strand = ask_required_choice(
                    "CD-HIT-EST strand comparison",
                    [
                        (
                            "both",
                            "Compare both +/+ and +/- orientations (-r 1).",
                        ),
                        (
                            "same",
                            "Compare only the same orientation, +/+ (-r 0).",
                        ),
                    ],
                )

        redundancy_result: Path | None = None
        redundancy_dir_name = {"seqkit": "seqkit", "cdhit": "cdhit_est"}.get(args.redundancy)
        if redundancy_dir_name is None:
            redundancy_output = current_fasta
            redundancy_results_dir = None
        else:
            redundancy_workdir = workdir / "preprocessing" / "redundancy" / redundancy_dir_name
            redundancy_results_dir = preprocessing_results_dir / "redundancy" / redundancy_dir_name
            redundancy_workdir.mkdir(parents=True, exist_ok=True)
            filename = (
                f"{project_name}_seqkit_deduplicated.fasta"
                if args.redundancy == "seqkit"
                else f"{project_name}_cdhit_filtered.fasta"
            )
            redundancy_output = redundancy_workdir / filename

        before_redundancy = count_fasta_sequences(current_fasta)
        word_size: int | None = None

        if args.redundancy == "none":
            pass  # MAFFT reads the unchanged FASTA directly; no extra copy.

        elif args.redundancy == "seqkit":
            check_required_tools([("SeqKit", args.seqkit_executable)])
            print("Running exact duplicate removal with SeqKit.")
            run_seqkit_deduplication(
                current_fasta,
                redundancy_output,
                args.seqkit_executable,
            )

        elif args.redundancy == "cdhit":
            check_required_tools(
                [("CD-HIT-EST", args.cdhit_executable)]
            )
            print("Running sequence clustering with CD-HIT-EST.")
            word_size = run_cdhit_clustering(
                current_fasta,
                redundancy_output,
                identity=args.identity,
                mode=args.cdhit_mode,
                strand=args.cdhit_strand,
                threads=args.threads,
                cdhit_executable=args.cdhit_executable,
            )

        else:
            raise RuntimeError(
                f"Unknown redundancy strategy: {args.redundancy}"
            )

        after_redundancy = count_fasta_sequences(redundancy_output)
        removed_redundancy = before_redundancy - after_redundancy
        reduction = (
            100.0 * removed_redundancy / before_redundancy
            if before_redundancy
            else 0.0
        )
        retained_percentage = (
            100.0 * after_redundancy / before_redundancy
            if before_redundancy
            else 0.0
        )

        section("Redundancy reduction result")
        print(f"Strategy                  : {args.redundancy}")
        if args.redundancy == "seqkit":
            print("Matching criterion         : exact sequence identity")
        elif args.redundancy == "cdhit":
            print(f"Identity threshold (-c)    : {args.identity}")
            print(f"Word size (-n)             : {word_size}")
            print(
                "Clustering mode (-g)       : "
                f"{args.cdhit_mode} "
                f"({'1' if args.cdhit_mode == 'accurate' else '0'})"
            )
            print(
                "Strand comparison (-r)     : "
                f"{args.cdhit_strand} "
                f"({'1' if args.cdhit_strand == 'both' else '0'})"
            )

        print(f"Sequences before          : {before_redundancy}")
        print(
            f"Sequences removed         : {removed_redundancy} "
            f"({reduction:.2f} %)"
        )
        print(
            f"Sequences retained        : {after_redundancy} "
            f"({retained_percentage:.2f} %)"
        )
        if args.redundancy == "cdhit":
            print(
                f"Cluster representatives   : {after_redundancy}"
            )

        if redundancy_results_dir is not None:
            redundancy_results_dir.mkdir(parents=True, exist_ok=True)
            redundancy_result = redundancy_results_dir / redundancy_output.name
            shutil.copy2(redundancy_output, redundancy_result)
            if args.redundancy == "cdhit":
                cluster_file = Path(str(redundancy_output) + ".clstr")
                if cluster_file.is_file():
                    shutil.copy2(cluster_file, redundancy_results_dir / cluster_file.name)
            report_lines = [
                "REDUNDANCY REDUCTION REPORT",
                "=" * 55,
                f"Strategy: {args.redundancy}",
                f"Sequences before: {before_redundancy}",
                f"Sequences removed: {removed_redundancy} ({reduction:.2f} %)",
                f"Sequences retained: {after_redundancy} ({retained_percentage:.2f} %)",
            ]
            if args.redundancy == "cdhit":
                report_lines.extend([
                    f"Identity threshold (-c): {args.identity}",
                    f"Word size (-n): {word_size}",
                    f"Clustering mode (-g): {args.cdhit_mode}",
                    f"Strand comparison (-r): {args.cdhit_strand}",
                ])
            (redundancy_results_dir / f"{project_name}_{redundancy_dir_name}_redundancy_report.txt").write_text(
                "\n".join(report_lines) + "\n", encoding="utf-8"
            )

        # ------------------------------------------------------------------
        # 4. Final MAFFT multiple sequence alignment
        # ------------------------------------------------------------------
        print("\n4. Final MAFFT alignment checkpoint")

        if args.mafft_strategy is None:
            if not interactive:
                raise RuntimeError(
                    "--mafft-strategy is required in non-interactive mode."
                )
            args.mafft_strategy = choose_mafft_strategy()

        if args.mafft_strategy not in MAFFT_STRATEGIES:
            raise ValueError(
                f"Unsupported MAFFT strategy: {args.mafft_strategy}"
            )

        if interactive:
            warn_about_mafft_strategy(
                args.mafft_strategy,
                after_redundancy,
            )

        check_required_tools([("MAFFT", args.mafft_executable)])
        final_alignment = mafft_workdir / f"{project_name}_mafft_aligned.fasta"

        print(
            f"Running final MAFFT alignment: {args.mafft_strategy} "
            f"({args.threads} threads)",
            flush=True,
        )
        mafft_reported_strategy = run_final_mafft(
            redundancy_output,
            final_alignment,
            args.mafft_strategy,
            args.threads,
            args.mafft_executable,
        )
        read_alignment(final_alignment)

        mafft_results_dir.mkdir(parents=True, exist_ok=True)
        final_alignment_result = (
            mafft_results_dir / f"{project_name}_mafft_aligned.fasta"
        )
        shutil.copy2(final_alignment, final_alignment_result)
        downstream_alignment_result = final_alignment_result

        pretrim_stats, pretrim_sequences = calculate_alignment_statistics(
            final_alignment
        )

        section("MAFFT alignment result")
        print(f"Sequences aligned         : {after_redundancy}")
        print(f"Requested strategy        : {args.mafft_strategy}")
        print(f"MAFFT reported strategy   : {mafft_reported_strategy}")
        print(
            f"Alignment length          : "
            f"{int(pretrim_stats['alignment_length'])} positions"
        )
        print(
            f"Mean gap content          : "
            f"{float(pretrim_stats['mean_gap_content_percent']):.2f} %"
        )
        print(
            f"Fully occupied columns    : "
            f"{int(pretrim_stats['fully_occupied_columns'])}"
        )

        # ------------------------------------------------------------------
        # 5. Alignment QC and trimming checkpoint
        # ------------------------------------------------------------------
        print("\n5. Alignment QC and trimming checkpoint")

        print_alignment_statistics(pretrim_stats)
        save_alignment_statistics(
            pretrim_stats,
            pretrim_sequences,
            mafft_qc_results_dir,
            project_name,
            "mafft",
        )
        section("Saved MAFFT results")
        print(f"Aligned FASTA             : {compact_path(final_alignment_result)}")
        print(
            "Alignment QC              : "
            + compact_path(mafft_qc_results_dir / f"{project_name}_mafft_alignment_qc.txt")
        )
        print(
            "Sequence occupancy        : "
            + compact_path(mafft_qc_results_dir / f"{project_name}_mafft_sequence_occupancy.csv")
        )

        if args.trimming is None:
            args.trimming = choose_trimming_strategy()

        if (
            args.trimming == "manual"
            and args.trimal_gap_threshold is None
        ):
            print("\nManual trimAl gap threshold")
            print("============================")
            print(
                "The -gt value is the minimum fraction of sequences without "
                "a gap in a retained column. Examples: 0.90 = >=90% "
                "occupancy; 1.00 = no gaps allowed."
            )
            args.trimal_gap_threshold = ask_required_float(
                "trimAl gap threshold",
                "-gt",
                minimum=0.0,
                maximum=1.0,
            )

            if args.trimal_gap_threshold < 0.50:
                print(
                    "Warning: this is a permissive threshold; many highly "
                    "gapped columns may be retained."
                )
                if not ask_yes_no(
                    f"Continue with -gt {args.trimal_gap_threshold}?"
                ):
                    raise RuntimeError(
                        "Manual trimAl threshold rejected by user."
                    )

        downstream_alignment = final_alignment
        trimmed_alignment_result: Path | None = None
        trimal_effect_file: Path | None = None

        if args.trimming != "none":
            check_required_tools(
                [("trimAl", args.trimal_executable)]
            )
            trimmed_alignment = (
                trimal_workdir / f"{project_name}_trimal_trimmed.fasta"
            )

            trimal_workdir.mkdir(parents=True, exist_ok=True)
            print(f"\nRunning trimAl strategy: {args.trimming}")
            run_trimal(
                final_alignment,
                trimmed_alignment,
                args.trimming,
                args.trimal_executable,
                args.trimal_gap_threshold,
            )
            downstream_alignment = trimmed_alignment

            trimal_results_dir.mkdir(parents=True, exist_ok=True)
            trimmed_alignment_result = (
                trimal_results_dir / f"{project_name}_trimal_trimmed.fasta"
            )
            shutil.copy2(trimmed_alignment, trimmed_alignment_result)
            downstream_alignment_result = trimmed_alignment_result

            posttrim_stats, posttrim_sequences = (
                calculate_alignment_statistics(trimmed_alignment)
            )

            save_alignment_statistics(
                posttrim_stats,
                posttrim_sequences,
                trimal_qc_results_dir,
                project_name,
                "trimal",
            )

            # Display and persist the same descriptive before/after assessment.
            trimal_effect_file = trimal_results_dir / f"{project_name}_trimal_effect.txt"
            comparison = format_trimal_effect(
                pretrim_stats,
                posttrim_stats,
                args.trimming,
                args.trimal_gap_threshold,
            )
            trimal_effect_file.write_text(comparison, encoding="utf-8")
            print("\n" + comparison, end="")

            print_alignment_statistics(
                posttrim_stats,
                title="Post-trimming alignment statistics",
            )

            section("Saved trimAl results")
            print(f"Trimmed alignment         : {compact_path(trimmed_alignment_result)}")
            print(f"Trimming effect report    : {compact_path(trimal_effect_file)}")
            print(
                "Post-trimming QC          : "
                + compact_path(trimal_qc_results_dir / f"{project_name}_trimal_alignment_qc.txt")
            )
            print(
                "Sequence occupancy        : "
                + compact_path(trimal_qc_results_dir / f"{project_name}_trimal_sequence_occupancy.csv")
            )
        else:
            section("trimAl result")
            print("Trimming                  : skipped")
            print("Alignment used downstream : complete MAFFT alignment")

        # ------------------------------------------------------------------
        # 6. Conservation analysis
        # ------------------------------------------------------------------
        print("\n6. Conservation analysis")
        selected_alignment_tool = "trimal" if args.trimming != "none" else "mafft"
        conservation_results_dir = (
            trimal_results_dir if selected_alignment_tool == "trimal"
            else mafft_results_dir
        ) / "conservation"
        perform_conservation_analysis(
            downstream_alignment,
            conservation_results_dir,
            project_name,
            args.min_occupancy,
            args.min_major_frequency,
        )

        # ------------------------------------------------------------------
        # 7. VarVAMP assay design with retry-on-failure
        # ------------------------------------------------------------------
        varvamp_results: Path | None = None
        varvamp_attempts: list[dict[str, object]] = []

        if not args.skip_varvamp:
            # Configuration directories are created only if a custom config is used.
            varvamp_work_root.mkdir(parents=True, exist_ok=True)
            print("\n7. VarVAMP assay design checkpoint")

            if args.varvamp_mode is None:
                if not interactive:
                    raise RuntimeError(
                        "--varvamp-mode is required in non-interactive mode."
                    )
                args.varvamp_mode = choose_varvamp_mode()

            if args.varvamp_threshold is None:
                if not interactive:
                    raise RuntimeError(
                        "--varvamp-threshold is required."
                    )
                args.varvamp_threshold = ask_required_float(
                    "VarVAMP consensus threshold",
                    "-t",
                    minimum=0.01,
                    maximum=1.0,
                )

            if args.primer_ambiguity is None:
                if not interactive:
                    raise RuntimeError(
                        "--primer-ambiguity is required."
                    )
                args.primer_ambiguity = (
                    ask_required_non_negative_integer(
                        "Maximum ambiguous bases in each primer",
                        "-a",
                    )
                )

            ensure_mode_required_interactive_parameters(
                args,
                args.varvamp_mode,
            )

            if interactive:
                configure_interactive_varvamp_advanced(
                    args,
                    args.varvamp_mode,
                )

            check_required_tools(
                [("VarVAMP", args.varvamp_executable)]
            )

            attempt_number = 0

            while True:
                attempt_number += 1
                mode = args.varvamp_mode

                ensure_mode_required_interactive_parameters(
                    args,
                    mode,
                )
                validate_varvamp_mode_parameters(args, mode)

                varvamp_workdir = (
                    varvamp_work_root
                    / f"attempt_{attempt_number:02d}_{mode}"
                )
                if varvamp_workdir.exists():
                    raise RuntimeError(
                        "Unexpected existing VarVAMP attempt directory: "
                        f"{varvamp_workdir}. Refusing to overwrite it."
                    )

                command = build_varvamp_command(
                    args,
                    mode,
                    downstream_alignment,
                    varvamp_workdir,
                )

                command_env = os.environ.copy()

                # Every new script execution starts with VarVAMP's normal
                # configuration unless a custom config is explicitly selected.
                command_env.pop("VARVAMP_CONFIG", None)

                config_path: Path | None = None
                if args.varvamp_config is not None:
                    config_path = stage_varvamp_config_for_run(
                        args.varvamp_config,
                        config_results_dir,
                        mode,
                        attempt_number,
                    )
                    # Keep the run-local copy as the canonical provenance path.
                    args.varvamp_config = config_path
                    command_env["VARVAMP_CONFIG"] = str(config_path)

                section(f"VarVAMP attempt {attempt_number}")
                print(f"Mode                      : {mode.upper()}")
                print(
                    f"Consensus threshold (-t)  : "
                    f"{args.varvamp_threshold}"
                )
                print(
                    f"Primer ambiguity (-a)     : "
                    f"{args.primer_ambiguity}"
                )
                if mode == "qpcr":
                    print(
                        f"Probe ambiguity (-pa)     : "
                        f"{args.probe_ambiguity}"
                    )
                print(
                    "Python configuration      : "
                    + (
                        compact_path(config_path)
                        if config_path is not None
                        else "VarVAMP default"
                    )
                )

                return_code, varvamp_output = run_varvamp_command(
                    command,
                    env=command_env,
                )

                diagnostics = parse_varvamp_diagnostics(
                    varvamp_output
                )

                attempt_record: dict[str, object] = {
                    "attempt": attempt_number,
                    "mode": mode,
                    "consensus_threshold": args.varvamp_threshold,
                    "primer_ambiguity": args.primer_ambiguity,
                    "probe_ambiguity": (
                        args.probe_ambiguity
                        if mode == "qpcr"
                        else None
                    ),
                    "single_opt_length": args.single_opt_length,
                    "single_max_length": args.single_max_length,
                    "single_report_n": args.single_report_n,
                    "tiled_opt_length": args.tiled_opt_length,
                    "tiled_max_length": args.tiled_max_length,
                    "tiled_overlap": args.tiled_overlap,
                    "qpcr_test_n": args.qpcr_test_n,
                    "qpcr_deltag": args.qpcr_deltag,
                    "custom_config": (
                        portable_manifest_path(config_path)
                        if config_path is not None
                        else None
                    ),
                    "return_code": return_code,
                    "status": (
                        "success"
                        if return_code == 0
                        else "failed"
                    ),
                }

                if return_code != 0:
                    attempt_record["reason"] = diagnostics.get(
                        "reason"
                    )

                varvamp_attempts.append(attempt_record)
                log_line(
                    "VARVAMP_ATTEMPT: "
                    + json.dumps(
                        attempt_record,
                        sort_keys=True,
                    )
                )

                if return_code == 0:
                    varvamp_results = copy_varvamp_tree(
                        varvamp_workdir,
                        results_dir,
                        mode,
                    )
                    print_varvamp_tables(
                        varvamp_results,
                        mode,
                    )
                    create_varvamp_visual_summaries(
                        varvamp_results
                    )
                    break

                if not interactive:
                    raise RuntimeError(
                        "VarVAMP failed in non-interactive mode: "
                        + str(
                            diagnostics.get(
                                "reason",
                                "unknown error",
                            )
                        )
                    )

                print_varvamp_failure(
                    args,
                    mode,
                    diagnostics,
                )

                action = adjust_varvamp_after_failure(
                    args,
                    mode,
                    config_results_dir,
                )

                if action == "stop":
                    print(
                        "VarVAMP design stopped by the user. "
                        "Previous workflow steps remain saved."
                    )
                    args.skip_varvamp = True
                    break

                # "retry" and "mode_changed" both return here and rerun only
                # the VarVAMP stage. All previous workflow products are reused.
                print("\nRetrying VarVAMP only...")
                print("Orientation              : already completed")
                print("Topology normalization   : already completed")
                print("Redundancy handling      : already completed")
                print("Final MAFFT              : already completed")
                print("trimAl / alignment QC    : already completed")
                print("Conservation analysis    : already completed")

        else:
            print("\n7. VarVAMP assay design skipped")

        manifest_file = results_dir / "assay_design_manifest.json"
        summary_file = results_dir / "workflow_summary.txt"

        write_project_manifest(
            manifest_file,
            project_name=project_name,
            run_id=run_id,
            input_file=input_file,
            workdir=workdir,
            results_dir=results_dir,
            downstream_alignment=downstream_alignment_result,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
            varvamp_attempts=varvamp_attempts,
        )

        inventory_file = results_dir / "output_inventory.tsv"
        record_design_artifacts(
            manifest_file,
            project_name=project_name,
            input_stats_file=input_stats_file,
            results_dir=results_dir,
            orientation_file=orientation_result,
            mars_file=mars_result,
            redundancy_file=redundancy_result,
            mafft_alignment_file=final_alignment_result,
            mafft_qc_dir=mafft_qc_results_dir,
            trimal_alignment_file=trimmed_alignment_result,
            trimal_qc_dir=trimal_qc_results_dir,
            trimal_effect_file=trimal_effect_file,
            conservation_dir=conservation_results_dir,
            selected_alignment_tool=selected_alignment_tool,
            inventory_file=inventory_file,
        )

        write_workflow_summary(
            summary_file,
            project_name=project_name,
            run_id=run_id,
            input_file=input_file,
            downstream_alignment=downstream_alignment_result,
            workdir=workdir,
            results_dir=results_dir,
            initial_count=initial_count,
            after_redundancy=after_redundancy,
            args=args,
            varvamp_results=varvamp_results,
        )

        if varvamp_results is not None:
            shutil.copy2(
                manifest_file,
                varvamp_results
                / "assay_design_manifest.json",
            )

        write_output_inventory(results_dir, inventory_file)

        section("Workflow completed")
        print(f"Run ID                   : {run_id}")
        print(f"Initial sequences        : {initial_count}")
        print(
            f"Final MAFFT sequences    : {after_redundancy}"
        )
        print(
            f"Downstream alignment     : "
            f"{compact_path(downstream_alignment_result)}"
        )

        if varvamp_results is not None:
            print(
                f"VarVAMP mode             : "
                f"{args.varvamp_mode.upper()}"
            )
            print(
                f"VarVAMP results          : "
                f"{compact_path(varvamp_results)}"
            )
        elif args.skip_varvamp:
            print("VarVAMP                  : no final design")
        else:
            print("VarVAMP                  : skipped")

        print(
            f"Project manifest         : "
            f"{compact_path(manifest_file)}"
        )
        print(
            f"Workflow summary         : "
            f"{compact_path(summary_file)}"
        )
        print(
            f"Workflow log             : "
            f"{compact_path(_run_management.WORKFLOW_LOG)}"
        )
        print(f"Permanent run results    : {compact_path(results_dir)}")
        print(f"Reconstructible work     : {compact_path(workdir)}")
        return 0

    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        pd.errors.ParserError,
    ) as error:
        log_line("ERROR: " + str(error))
        print(f"\nERROR: {error}", file=sys.stderr)
        return 1


