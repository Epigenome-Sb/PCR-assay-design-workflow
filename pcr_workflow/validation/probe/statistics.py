# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv

from ...core import runtime
from ..fasta import read_fasta
from .reporting import write_probe_variant_outputs

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

    stats_file, matrix_file = write_probe_variant_outputs(
        stats_rows, matrix, variants, pcr_positive_hitids, probe_length,
        scheme_dir, category_rank,
    )
    return {
        "stats_rows": stats_rows,
        "stats_file": stats_file,
        "matrix_file": matrix_file,
    }
