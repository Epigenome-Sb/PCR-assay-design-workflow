# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import shutil

from ..core.runtime import section
from .alignment import calculate_alignment_statistics, print_alignment_statistics, read_alignment, run_final_mafft, run_trimal, save_alignment_statistics
from .cli import ask_required_choice, ask_required_float, ask_yes_no
from .commands import check_required_tools
from .conservation import perform_conservation_analysis
from .constants import MAFFT_STRATEGIES
from .fasta import calculate_fasta_statistics, count_fasta_sequences, print_fasta_statistics
from .orientation import calculate_cyclic_rotation_statistics, normalize_orientation_with_mafft, print_rotation_statistics, run_mars_rotation
from .redundancy import run_cdhit_clustering, run_seqkit_deduplication
from .ui import choose_mafft_strategy, choose_orientation_strategy, choose_redundancy_strategy, choose_sequence_topology, choose_trimming_strategy, warn_about_mafft_strategy

def run_preprocessing(
    args,
    interactive: bool,
    input_file: Path,
    project_name: str,
    preprocessing_workdir: Path,
    preprocessing_results_dir: Path,
) -> tuple[Path, Path, int, int]:
    # ------------------------------------------------------------------
    # Input overview
    # ------------------------------------------------------------------
    current_fasta = input_file
    initial_count = count_fasta_sequences(input_file)

    print_fasta_statistics(
        calculate_fasta_statistics(input_file),
        title="Input dataset statistics",
    )

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

    if args.orientation == "keep":
        section("Orientation result")
        print("Method                    : keep original orientation")
        print("MAFFT orientation check   : skipped")
        print(f"Sequences retained        : {initial_count}")
        print("Sequences changed         : 0")
    else:
        check_required_tools([("MAFFT", args.mafft_executable)])
        print("Running sequence orientation normalization with MAFFT.")
        current_fasta, reversed_count = normalize_orientation_with_mafft(
            current_fasta,
            preprocessing_workdir,
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
        oriented_result = (
            preprocessing_results_dir / f"{project_name}_oriented.fasta"
        )
        shutil.copy2(current_fasta, oriented_result)
        orientation_report = (
            preprocessing_workdir / f"{project_name}_orientation_report.txt"
        )
        if orientation_report.is_file():
            shutil.copy2(
                orientation_report,
                preprocessing_results_dir / orientation_report.name,
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

    if args.topology == "circular":
        check_required_tools([("MARS", args.mars_executable)])
        print("Running circular sequence rotation with MARS.")
        before_mars = current_fasta
        current_fasta = run_mars_rotation(
            current_fasta,
            preprocessing_workdir,
            project_name,
            args.mars_executable,
            args.threads,
        )
        print_rotation_statistics(
            calculate_cyclic_rotation_statistics(
                before_mars,
                current_fasta,
            )
        )
    elif args.topology == "linear":
        section("Topology result")
        print("Topology                  : linear")
        print("MARS                      : not applicable")
        print("Cyclic starts changed     : 0")
    else:
        raise ValueError("--topology must be 'linear' or 'circular'.")

    if args.topology == "circular":
        shutil.copy2(
            current_fasta,
            preprocessing_results_dir / f"{project_name}_rotated.fasta",
        )

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

    redundancy_output = (
        preprocessing_workdir / f"{project_name}_redundancy_filtered.fasta"
    )
    before_redundancy = count_fasta_sequences(current_fasta)
    word_size: int | None = None

    if args.redundancy == "none":
        shutil.copy2(current_fasta, redundancy_output)

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

    # Keep the biologically meaningful preprocessing output for provenance.
    redundancy_result = (
        preprocessing_results_dir
        / f"{project_name}_redundancy_filtered.fasta"
    )
    shutil.copy2(redundancy_output, redundancy_result)
    if args.redundancy == "cdhit":
        cluster_file = Path(str(redundancy_output) + ".clstr")
        if cluster_file.is_file():
            shutil.copy2(
                cluster_file,
                preprocessing_results_dir / cluster_file.name,
            )
    return current_fasta, redundancy_output, initial_count, after_redundancy


def run_alignment_and_conservation(
    args,
    interactive: bool,
    redundancy_output: Path,
    after_redundancy: int,
    project_name: str,
    alignment_workdir: Path,
    alignment_results_dir: Path,
    conservation_results_dir: Path,
) -> tuple[Path, Path]:
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
    final_alignment = alignment_workdir / f"{project_name}_alignment.fasta"

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

    final_alignment_result = (
        alignment_results_dir / f"{project_name}_alignment.fasta"
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
        alignment_results_dir,
        project_name,
        "pretrim",
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

    if args.trimming != "none":
        check_required_tools(
            [("trimAl", args.trimal_executable)]
        )
        trimmed_alignment = (
            alignment_workdir / f"{project_name}_alignment_trimmed.fasta"
        )

        print(f"\nRunning trimAl strategy: {args.trimming}")
        run_trimal(
            final_alignment,
            trimmed_alignment,
            args.trimming,
            args.trimal_executable,
            args.trimal_gap_threshold,
        )
        downstream_alignment = trimmed_alignment

        trimmed_alignment_result = (
            alignment_results_dir
            / f"{project_name}_alignment_trimmed.fasta"
        )
        shutil.copy2(trimmed_alignment, trimmed_alignment_result)
        downstream_alignment_result = trimmed_alignment_result

        posttrim_stats, posttrim_sequences = (
            calculate_alignment_statistics(trimmed_alignment)
        )

        save_alignment_statistics(
            posttrim_stats,
            posttrim_sequences,
            alignment_results_dir,
            project_name,
            "posttrim",
        )

        original_length = int(pretrim_stats["alignment_length"])
        trimmed_length = int(posttrim_stats["alignment_length"])
        removed_columns = original_length - trimmed_length

        section("trimAl effect")
        print(f"Strategy                  : {args.trimming}")
        if args.trimming == "manual":
            print(
                f"Gap threshold (-gt)       : "
                f"{args.trimal_gap_threshold}"
            )
        print(f"Columns before            : {original_length}")
        print(f"Columns after             : {trimmed_length}")
        print(
            f"Columns removed           : {removed_columns} "
            f"({100.0 * removed_columns / original_length:.2f} %)"
        )
        print(
            f"Mean gap content before   : "
            f"{float(pretrim_stats['mean_gap_content_percent']):.2f} %"
        )
        print(
            f"Mean gap content after    : "
            f"{float(posttrim_stats['mean_gap_content_percent']):.2f} %"
        )

        print_alignment_statistics(
            posttrim_stats,
            title="Post-trimming alignment statistics",
        )
    else:
        section("trimAl result")
        print("Trimming                  : skipped")
        print("Alignment used downstream : complete MAFFT alignment")

    # ------------------------------------------------------------------
    # 6. Conservation analysis
    # ------------------------------------------------------------------
    print("\n6. Conservation analysis")
    perform_conservation_analysis(
        downstream_alignment,
        conservation_results_dir,
        project_name,
        args.min_occupancy,
        args.min_major_frequency,
    )
    return downstream_alignment, downstream_alignment_result
