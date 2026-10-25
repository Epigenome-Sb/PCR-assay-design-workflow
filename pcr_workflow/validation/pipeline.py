# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

import shutil
import sys

from ..core import runtime
from ..core.paths import compact_path
from ..core.runtime import automatic_thread_count, available_cpu_threads, log_line
from .assay_inputs import find_qpcr_oligos, prepare_single_or_tiled_files, prepare_varvamp_files
from .cli import parse_cli
from .database import ensure_validation_tools_available, format_target_database_for_mfeprimer
from .fasta import fasta_statistics
from .manifests import resolve_validation_run_directories, write_assay_manifest, write_validation_database_summary, write_validation_manifest
from .mfeprimer import run_mfeprimer_pairs
from .selection import choose_design_run
from .ui import ask_yes_no, choose_mfeprimer_parameters, choose_schemes, choose_sequence_length_filter, choose_target_database, index_target_database

def main() -> int:
    cli = parse_cli()
    runtime.VERBOSE = bool(cli.verbose)

    try:
        cpu_threads = available_cpu_threads()
        if cli.threads is None:
            cli.threads = automatic_thread_count(cpu_threads)
            thread_selection = "automatic"
        else:
            thread_selection = "manual"

        if cli.threads <= 0:
            raise ValueError("--threads must be greater than zero.")

        if cli.threads > cpu_threads:
            print(
                f"Warning: --threads={cli.threads} exceeds the "
                f"{cpu_threads} logical CPU threads currently available."
            )

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

        runtime.LOG_PATH = results_dir / "validation.log"
        runtime.LOG_PATH.write_text(
            "VarVAMP in silico validation log\n"
            f"Project: {project_name}\n"
            f"Validation run ID: {run_id}\n"
            f"Parent design run ID: {design_run_id}\n"
            f"Design manifest: {compact_path(design_manifest_path)}\n"
            f"CPU threads available: {cpu_threads}\n"
            f"Workflow threads: {cli.threads} ({thread_selection})\n"
            f"SeqKit formatting threads: {min(cli.threads, 4)}\n"
            f"BLAST threads: {cli.threads if assay_mode == 'qpcr' else 'N/A'}\n",
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
        print(f"CPU         : {cpu_threads} logical threads available")
        print(
            f"Threads     : {cli.threads} used where supported "
            f"({thread_selection})"
        )
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
            cli.threads,
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
        print(
            f"CPU threads       : {cli.threads} for BLAST/supported tools; "
            f"SeqKit formatting uses {min(cli.threads, 4)}"
        )
        print(
            "MFEprimer threads : not user-configurable in the supported "
            "MFEprimer 3.1 full-QC command"
        )

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
            cli.threads,
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
            threads=cli.threads,
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
