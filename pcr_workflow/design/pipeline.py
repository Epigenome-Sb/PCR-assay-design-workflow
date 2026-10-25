# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

import shutil
import sys

import pandas as pd

from ..core import runtime
from ..core.paths import compact_path, resolve_run_directories, sanitize_project_name
from ..core.runtime import automatic_thread_count, available_cpu_threads, log_line, section
from .cli import ask_for_input_file, parse_arguments
from .fasta import validate_fraction, validate_nucleotide_fasta
from .manifests import write_project_manifest, write_workflow_summary
from .stages import run_alignment_and_conservation, run_preprocessing
from .varvamp_stage import run_varvamp_stage

def main() -> int:
    args = parse_arguments()
    runtime.VERBOSE = bool(args.verbose)

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

        # Run-specific working areas contain reconstructible/intermediate files.
        preprocessing_workdir = workdir / "preprocessing"
        alignment_workdir = workdir / "alignment"
        varvamp_work_root = workdir / "varvamp"

        # Run-specific result areas contain permanent scientific outputs.
        preprocessing_results_dir = results_dir / "preprocessing"
        alignment_results_dir = results_dir / "alignment"
        conservation_results_dir = results_dir / "conservation"
        config_results_dir = results_dir / "config"

        for directory in (
            preprocessing_workdir,
            alignment_workdir,
            varvamp_work_root,
            preprocessing_results_dir,
            alignment_results_dir,
            conservation_results_dir,
            config_results_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        runtime.LOG_PATH = results_dir / "workflow.log"
        runtime.LOG_PATH.write_text(
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
        if runtime.VERBOSE:
            print(f"Log     : {runtime.LOG_PATH}")

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

        _, redundancy_output, initial_count, after_redundancy = run_preprocessing(
            args, interactive, input_file, project_name,
            preprocessing_workdir, preprocessing_results_dir,
        )
        downstream_alignment, downstream_alignment_result = (
            run_alignment_and_conservation(
                args, interactive, redundancy_output, after_redundancy,
                project_name, alignment_workdir, alignment_results_dir,
                conservation_results_dir,
            )
        )
        varvamp_results, varvamp_attempts = run_varvamp_stage(
            args, interactive, downstream_alignment, varvamp_work_root,
            config_results_dir, results_dir,
        )

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
            f"{compact_path(runtime.LOG_PATH)}"
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
