# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv
import re

from ...core import runtime
from .matching import extract_full_length_mismatches, format_probe_blast_alignment_block

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

def write_probe_variant_outputs(
    stats_rows: list[dict[str, object]],
    matrix: dict[str, dict[str, str]],
    variants: list[dict[str, str]],
    pcr_positive_hitids: list[str],
    probe_length: int,
    scheme_dir: Path,
    category_rank: dict[str, int],
) -> tuple[Path, Path]:
    """Write per-variant TSVs and the optional expanded terminal table."""
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
    if runtime.VERBOSE and len(stats_rows) > 1:
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

