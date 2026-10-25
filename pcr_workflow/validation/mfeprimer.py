# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import csv

from ..core import runtime
from ..core.paths import compact_path, optional_portable_project_path, portable_project_path
from .amplicons import analyse_and_extract_mfeprimer_amplicons
from .commands import run_command
from .primer_qc import parse_mfeprimer_secondary_structures, print_secondary_structure_qc, save_secondary_structure_qc
from .probe.pipeline import analyse_probe_with_blast, empty_probe_analysis
from .probe.reporting import print_probe_variant_amplicon_table
from .reporting import inspect_mfeprimer_products, print_compact_primer_table, print_general_potential_amplicon_statistics, print_potential_amplicons_by_pair
from .runtime import progress
from .ui import ask_yes_no

def _run_mfeprimer_batch(
    schemes: dict[str, dict[str, str]],
    selected_schemes: list[str],
    assay_inputs_dir: Path,
    mfeprimer_results_dir: Path,
    probe_results_dir: Path,
    summary_results_dir: Path,
    database: Path,
    parameters: dict[str, float | int] | None,
    total_target_sequences: int,
) -> dict[str, object]:
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

    return {
        "run_manifest": run_manifest, "reports": reports, "scheme_dirs": scheme_dirs,
        "previews": previews, "secondary_by_scheme": secondary_by_scheme,
    }



def _review_potential_products(
    batch: dict[str, object],
    selected_schemes: list[str],
    assay_mode: str,
    total_target_sequences: int,
    summary_results_dir: Path,
) -> tuple[bool, Path, Path]:
    reports = batch["reports"]
    previews = batch["previews"]
    secondary_by_scheme = batch["secondary_by_scheme"]

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

    return exclude_self_priming, secondary_summary_file, secondary_details_file



def _filter_amplicons(
    batch: dict[str, object],
    selected_schemes: list[str],
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
    exclude_self_priming: bool,
) -> dict[str, dict[str, object]]:
    reports = batch["reports"]
    scheme_dirs = batch["scheme_dirs"]
    previews = batch["previews"]

    # Apply the user's explicit decision only after the unfiltered MFEprimer
    # results have been displayed.
    analyses: dict[str, dict[str, object]] = {}
    for scheme in selected_schemes:
        if runtime.VERBOSE:
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

    if runtime.VERBOSE:
        print(
            "\nCoverage set: "
            + (
                "TARGET only (SELF_PRIMING excluded)."
                if exclude_self_priming
                else "TARGET + SELF_PRIMING."
            )
        )

    return analyses



def _analyse_probe_batch(
    schemes: dict[str, dict[str, str]],
    selected_schemes: list[str],
    analyses: dict[str, dict[str, object]],
    probe_results_dir: Path,
    total_target_sequences: int,
    assay_mode: str,
    threads: int,
) -> tuple[dict[str, dict[str, object]], bool]:
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
                    if runtime.VERBOSE:
                        blast_analysis = analyse_probe_with_blast(
                            scheme,
                            schemes[scheme]["PROBE"],
                            analyses[scheme]["valid_fasta"],
                            analyses[scheme]["valid_rows"],
                            probe_scheme_dir,
                            total_target_sequences,
                            threads,
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
                                threads,
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

    return probe_results, probe_validation_enabled



def _build_coverage_manifest(
    selected_schemes: list[str],
    analyses: dict[str, dict[str, object]],
    probe_results: dict[str, dict[str, object]],
    assay_mode: str,
    total_target_sequences: int,
    probe_validation_enabled: bool,
) -> list[list[object]]:
    coverage_manifest: list[list[object]] = []
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

    return coverage_manifest



def _write_validation_outputs(
    coverage_manifest: list[list[object]],
    run_manifest: list[list[str]],
    summary_results_dir: Path,
    assay_mode: str,
    total_target_sequences: int,
    exclude_self_priming: bool,
    probe_validation_enabled: bool,
    secondary_summary_file: Path,
    secondary_details_file: Path,
) -> dict[str, object]:
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
    if runtime.LOG_PATH is not None:
        print(f"Log     : {compact_path(runtime.LOG_PATH)}")
    if not runtime.VERBOSE:
        print("Details : saved to files; use --verbose for expanded terminal output.")

    return {
        "exclude_self_priming": exclude_self_priming,
        "probe_validation_enabled": probe_validation_enabled,
        "coverage_results": coverage_manifest_file,
        "mfeprimer_runs": mfeprimer_runs_file,
        "secondary_summary": secondary_summary_file,
        "secondary_details": secondary_details_file,
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
    threads: int,
) -> dict[str, object]:
    """Run primer QC, apply the global product filter, and summarize coverage."""
    batch = _run_mfeprimer_batch(
        schemes, selected_schemes, assay_inputs_dir, mfeprimer_results_dir,
        probe_results_dir, summary_results_dir, database, parameters,
        total_target_sequences,
    )
    exclude_self_priming, secondary_summary_file, secondary_details_file = (
        _review_potential_products(
            batch, selected_schemes, assay_mode, total_target_sequences,
            summary_results_dir,
        )
    )
    analyses = _filter_amplicons(
        batch, selected_schemes, total_target_sequences, min_target_size,
        max_target_size, exclude_self_priming,
    )
    probe_results, probe_validation_enabled = _analyse_probe_batch(
        schemes, selected_schemes, analyses, probe_results_dir,
        total_target_sequences, assay_mode, threads,
    )
    coverage_manifest = _build_coverage_manifest(
        selected_schemes, analyses, probe_results, assay_mode,
        total_target_sequences, probe_validation_enabled,
    )
    return _write_validation_outputs(
        coverage_manifest, batch["run_manifest"], summary_results_dir,
        assay_mode, total_target_sequences, exclude_self_priming,
        probe_validation_enabled, secondary_summary_file, secondary_details_file,
    )

