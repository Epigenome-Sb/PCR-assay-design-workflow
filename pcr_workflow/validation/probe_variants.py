
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Probe IUPAC-variant statistics and terminal summaries for Workflow 02."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from . import run_management as _run_management
from .fasta import read_fasta


def summarize_probe_variants(
    query_file: Path,
    raw_hits: list[dict[str, object]],
    valid_rows: list[dict[str, object]],
    probe_length: int,
    scheme_dir: Path,
    category_rank: dict[str, int],
) -> dict[str, object]:
    """
    Analyse every concrete IUPAC-expanded probe variant independently.

    Statistics are calculated at unique HitID level:
      - each concrete probe variant is evaluated against every PCR-positive HitID;
      - if one HitID has multiple valid MFEprimer amplicons, the best result for
        that concrete variant is retained once;
      - variants are NOT merged in this table.

    A separate best-of-all-variants summary is still produced elsewhere and is
    used for the final qPCR coverage.
    """
    query_records = read_fasta(query_file)

    variants = [
        {
            "query_id": header,
            "sequence": sequence,
        }
        for header, sequence in query_records
    ]

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    # Stable first-occurrence order of unique biological PCR-positive targets.
    pcr_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    # Hits have already been classified as full-length / partial by the caller.
    hits_by_variant_amp: dict[
        tuple[str, int],
        list[dict[str, object]],
    ] = {}

    for hit in raw_hits:
        amp_id = hit.get("AmpID")
        if amp_id is None:
            continue

        amp_id = int(amp_id)
        if amp_id not in valid_by_ampid:
            continue

        query_id = str(hit["qseqid"])
        hits_by_variant_amp.setdefault(
            (query_id, amp_id),
            [],
        ).append(hit)

    def choose_best_amp_result(
        query_id: str,
        amp_id: int,
    ) -> dict[str, object]:
        amp_hits = hits_by_variant_amp.get(
            (query_id, amp_id),
            [],
        )

        full_hits = [
            hit
            for hit in amp_hits
            if hit.get("FullLengthUngapped") == "YES"
        ]

        if full_hits:
            full_hits.sort(
                key=lambda hit: (
                    int(hit["mismatch"]),
                    -float(hit["pident"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = full_hits[0]
            return {
                "category": str(best["BlastCategory"]),
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        if amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            return {
                "category": "PARTIAL_ONLY",
                "mismatches": int(best["mismatch"]),
                "alignment_length": int(best["length"]),
                "bitscore": float(best["bitscore"]),
                "amp_id": amp_id,
            }

        return {
            "category": "NO_HIT",
            "mismatches": 999999,
            "alignment_length": 0,
            "bitscore": 0.0,
            "amp_id": amp_id,
        }

    def result_rank(result: dict[str, object]) -> tuple:
        return (
            category_rank[str(result["category"])],
            int(result["mismatches"]),
            -int(result["alignment_length"]),
            -float(result["bitscore"]),
            int(result["amp_id"]),
        )

    # Map each HitID to all its valid amplicons.
    ampids_by_hitid: dict[str, list[int]] = {}
    for row in valid_rows:
        ampids_by_hitid.setdefault(
            str(row["HitID"]),
            [],
        ).append(int(row["AmpID"]))

    matrix: dict[str, dict[str, str]] = {
        hit_id: {}
        for hit_id in pcr_positive_hitids
    }

    stats_rows: list[dict[str, object]] = []

    for variant in variants:
        query_id = str(variant["query_id"])
        sequence = str(variant["sequence"])

        # A. Statistics by unique PCR-positive HitID.
        #    These remain the basis of biological coverage so one biological
        #    target is not counted several times when MFEprimer reports more
        #    than one retained amplicon for the same HitID.
        counts = {
            category: 0
            for category in category_rank
        }

        for hit_id in pcr_positive_hitids:
            amp_results = [
                choose_best_amp_result(
                    query_id,
                    amp_id,
                )
                for amp_id in ampids_by_hitid[hit_id]
            ]

            best = min(
                amp_results,
                key=result_rank,
            )
            category = str(best["category"])

            counts[category] += 1
            matrix[hit_id][query_id] = category

        full_length = (
            counts["0_MISMATCH"]
            + counts["1_MISMATCH"]
            + counts["2_MISMATCHES"]
            + counts["GE3_MISMATCHES"]
        )

        denominator = len(pcr_positive_hitids)

        # B. Statistics by retained TARGET amplicon.
        #    This is displayed immediately after probe BLAST so the user can
        #    see, for every concrete probe version, how many of the retained
        #    amplicons are exact, 1-mismatch, 2-mismatch, etc.
        amplicon_counts = {
            category: 0
            for category in category_rank
        }

        for row in valid_rows:
            amp_id = int(row["AmpID"])
            best_amp = choose_best_amp_result(query_id, amp_id)
            amplicon_counts[str(best_amp["category"])] += 1

        amplicon_full_length = (
            amplicon_counts["0_MISMATCH"]
            + amplicon_counts["1_MISMATCH"]
            + amplicon_counts["2_MISMATCHES"]
            + amplicon_counts["GE3_MISMATCHES"]
        )

        stats_rows.append(
            {
                "QueryVariant": query_id,
                "ConcreteSequence": sequence,

                # Retained TARGET amplicon-level counts
                "TargetAmplicons": len(valid_rows),
                "AmpliconExact0Mismatch": amplicon_counts["0_MISMATCH"],
                "AmpliconOneMismatch": amplicon_counts["1_MISMATCH"],
                "AmpliconTwoMismatches": amplicon_counts["2_MISMATCHES"],
                "AmpliconGE3Mismatches": amplicon_counts["GE3_MISMATCHES"],
                "AmpliconPartialOnly": amplicon_counts["PARTIAL_ONLY"],
                "AmpliconNoHit": amplicon_counts["NO_HIT"],
                "AmpliconFullLengthSite": amplicon_full_length,

                # Unique-HitID counts used for coverage
                "PCRPositiveHitIDs": denominator,
                "Exact0Mismatch": counts["0_MISMATCH"],
                "OneMismatch": counts["1_MISMATCH"],
                "TwoMismatches": counts["2_MISMATCHES"],
                "GE3Mismatches": counts["GE3_MISMATCHES"],
                "PartialOnly": counts["PARTIAL_ONLY"],
                "NoHit": counts["NO_HIT"],
                "FullLengthSite": full_length,
                "ExactPercentPCRPositive": (
                    100.0 * counts["0_MISMATCH"] / denominator
                    if denominator
                    else 0.0
                ),
                "LE1MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "LE2MismatchPercentPCRPositive": (
                    100.0
                    * (
                        counts["0_MISMATCH"]
                        + counts["1_MISMATCH"]
                        + counts["2_MISMATCHES"]
                    )
                    / denominator
                    if denominator
                    else 0.0
                ),
                "FullLengthPercentPCRPositive": (
                    100.0 * full_length / denominator
                    if denominator
                    else 0.0
                ),
            }
        )

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
    if _run_management.VERBOSE and len(stats_rows) > 1:
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
