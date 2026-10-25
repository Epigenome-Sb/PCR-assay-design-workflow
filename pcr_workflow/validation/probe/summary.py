# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv

from ...core import runtime
from ...core.paths import compact_path
from ...core.runtime import section
from .blast import read_probe_blast_raw
from .matching import classify_full_length_probe_hit
from .reporting import write_probe_blast_alignment_reports
from .statistics import summarize_probe_variants

def _build_probe_summary_rows(
    scheme: str,
    probe: str,
    raw_file: Path,
    query_file: Path,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
):
    probe_length = len(probe)
    raw_hits = read_probe_blast_raw(raw_file)

    valid_by_ampid = {
        int(row["AmpID"]): row
        for row in valid_rows
    }

    hits_by_ampid: dict[int, list[dict[str, object]]] = {}

    for hit in raw_hits:
        amp_id = hit["AmpID"]

        if amp_id is None:
            continue

        amp_id = int(amp_id)

        if amp_id not in valid_by_ampid:
            continue

        full_length, category = classify_full_length_probe_hit(
            hit,
            probe_length,
        )

        hit["FullLengthUngapped"] = (
            "YES" if full_length else "NO"
        )
        hit["BlastCategory"] = category

        hits_by_ampid.setdefault(amp_id, []).append(hit)

    category_rank = {
        "0_MISMATCH": 0,
        "1_MISMATCH": 1,
        "2_MISMATCHES": 2,
        "GE3_MISMATCHES": 3,
        "PARTIAL_ONLY": 4,
        "NO_HIT": 5,
    }

    variant_analysis = summarize_probe_variants(
        query_file,
        raw_hits,
        valid_rows,
        probe_length,
        scheme_dir,
        category_rank,
    )

    amplicon_summary: list[dict[str, object]] = []

    for amp_id in sorted(valid_by_ampid):
        source = valid_by_ampid[amp_id]
        amp_hits = hits_by_ampid.get(amp_id, [])

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
            category = str(best["BlastCategory"])
        elif amp_hits:
            amp_hits.sort(
                key=lambda hit: (
                    -int(hit["length"]),
                    int(hit["mismatch"]),
                    -float(hit["bitscore"]),
                    float(hit["evalue"]),
                )
            )
            best = amp_hits[0]
            category = "PARTIAL_ONLY"
        else:
            best = None
            category = "NO_HIT"

        if best is None:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": "NO",
                "QueryVariant": "",
                "PercentIdentity": "",
                "AlignmentLength": "",
                "Mismatches": "",
                "GapOpenings": "",
                "QueryStart": "",
                "QueryEnd": "",
                "SubjectStart": "",
                "SubjectEnd": "",
                "Evalue": "",
                "BitScore": "",
                "QueryAligned": "",
                "SubjectAligned": "",
                "RawHSPCount": 0,
            }
        else:
            row = {
                "AmpID": amp_id,
                "HitID": str(source["HitID"]),
                "AmpliconSize": int(source["Size"]),
                "BlastCategory": category,
                "FullLengthUngapped": str(
                    best.get("FullLengthUngapped", "NO")
                ),
                "QueryVariant": str(best["qseqid"]),
                "PercentIdentity": float(best["pident"]),
                "AlignmentLength": int(best["length"]),
                "Mismatches": int(best["mismatch"]),
                "GapOpenings": int(best["gapopen"]),
                "QueryStart": int(best["qstart"]),
                "QueryEnd": int(best["qend"]),
                "SubjectStart": int(best["sstart"]),
                "SubjectEnd": int(best["send"]),
                "Evalue": float(best["evalue"]),
                "BitScore": float(best["bitscore"]),
                "QueryAligned": str(best["qseq"]),
                "SubjectAligned": str(best["sseq"]),
                "RawHSPCount": len(amp_hits),
            }

        amplicon_summary.append(row)

    # One best amplicon result per unique target HitID.
    hitid_best: dict[str, dict[str, object]] = {}

    def result_rank(row: dict[str, object]) -> tuple:
        category = str(row["BlastCategory"])

        mismatch_value = (
            int(row["Mismatches"])
            if row["Mismatches"] != ""
            else 999999
        )

        alignment_length = (
            int(row["AlignmentLength"])
            if row["AlignmentLength"] != ""
            else 0
        )

        bitscore = (
            float(row["BitScore"])
            if row["BitScore"] != ""
            else 0.0
        )

        return (
            category_rank[category],
            mismatch_value,
            -alignment_length,
            -bitscore,
            int(row["AmpID"]),
        )

    for row in amplicon_summary:
        hit_id = str(row["HitID"])

        if (
            hit_id not in hitid_best
            or result_rank(row) < result_rank(hitid_best[hit_id])
        ):
            hitid_best[hit_id] = row

    hitid_summary = list(hitid_best.values())
    hitid_summary.sort(
        key=lambda row: str(row["HitID"])
    )

    (
        alignment_amplicon_file,
        alignment_hitid_file,
        mismatch_positions_file,
    ) = write_probe_blast_alignment_reports(
        amplicon_summary,
        hitid_summary,
        scheme_dir,
    )

    return (
        probe_length, raw_hits, category_rank, variant_analysis,
        amplicon_summary, hitid_summary, alignment_amplicon_file,
        alignment_hitid_file, mismatch_positions_file,
    )


def summarize_probe_blast(
    scheme: str,
    probe: str,
    raw_file: Path,
    query_file: Path,
    concrete_query_count: int,
    valid_rows: list[dict[str, object]],
    scheme_dir: Path,
    total_target_sequences: int,
) -> dict[str, object]:
    """
    Produce BLAST-only summaries.

    There is intentionally:
    - no IUPAC optimization;
    - no ambiguity proposal;
    - no PASS/FAIL probe decision.
    """
    (
        probe_length, raw_hits, category_rank, variant_analysis,
        amplicon_summary, hitid_summary, alignment_amplicon_file,
        alignment_hitid_file, mismatch_positions_file,
    ) = _build_probe_summary_rows(
        scheme, probe, raw_file, query_file, valid_rows, scheme_dir
    )

    counts = {
        category: 0
        for category in category_rank
    }

    for row in hitid_summary:
        counts[str(row["BlastCategory"])] += 1

    amplicon_counts = {
        category: 0
        for category in category_rank
    }

    for row in amplicon_summary:
        amplicon_counts[str(row["BlastCategory"])] += 1

    full_length_amplicons = (
        amplicon_counts["0_MISMATCH"]
        + amplicon_counts["1_MISMATCH"]
        + amplicon_counts["2_MISMATCHES"]
        + amplicon_counts["GE3_MISMATCHES"]
    )

    exact_probe_amplicons = amplicon_counts["0_MISMATCH"]

    full_length_hitids = (
        counts["0_MISMATCH"]
        + counts["1_MISMATCH"]
        + counts["2_MISMATCHES"]
        + counts["GE3_MISMATCHES"]
    )

    exact_probe_hitids = counts["0_MISMATCH"]
    pcr_positive_hitids = len(hitid_summary)
    total_valid_amplicons = len(amplicon_summary)

    # -------------------- write processed outputs --------------------
    columns = [
        "AmpID",
        "HitID",
        "AmpliconSize",
        "BlastCategory",
        "FullLengthUngapped",
        "QueryVariant",
        "PercentIdentity",
        "AlignmentLength",
        "Mismatches",
        "GapOpenings",
        "QueryStart",
        "QueryEnd",
        "SubjectStart",
        "SubjectEnd",
        "Evalue",
        "BitScore",
        "QueryAligned",
        "SubjectAligned",
        "RawHSPCount",
    ]

    amplicon_file = scheme_dir / "probe_blast_amplicon_summary.tsv"
    hitid_file = scheme_dir / "probe_blast_hitid_summary.tsv"
    summary_file = scheme_dir / "probe_blast_summary.txt"

    def write_table(
        path: Path,
        rows: list[dict[str, object]],
    ) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_table(amplicon_file, amplicon_summary)
    write_table(hitid_file, hitid_summary)

    summary_lines = [
        f"BLAST probe summary — {scheme}",
        "=" * (22 + len(scheme)),
        "",
        f"Original VarVAMP probe: {probe}",
        f"Probe length: {probe_length} nt",
        f"Concrete BLAST query variants: {concrete_query_count}",
        f"Valid MFEprimer target amplicons: {len(valid_rows)}",
        f"Unique PCR-positive HitIDs: {pcr_positive_hitids}",
        f"Raw BLAST HSPs: {len(raw_hits)}",
        "",
        "AMPlicon-level probe statistics:",
        f"Valid target amplicons analysed: {total_valid_amplicons}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * exact_probe_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_amplicons} / {total_valid_amplicons} "
            f"= {(100.0 * full_length_amplicons / total_valid_amplicons if total_valid_amplicons else 0.0):.2f} %"
        ),
        f"0 mismatch: {amplicon_counts['0_MISMATCH']}",
        f"1 mismatch: {amplicon_counts['1_MISMATCH']}",
        f"2 mismatches: {amplicon_counts['2_MISMATCHES']}",
        f">=3 mismatches: {amplicon_counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {amplicon_counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {amplicon_counts['NO_HIT']}",
        "",
        "HITID-level probe statistics (unique PCR-positive targets):",
        f"Unique PCR-positive HitIDs analysed: {pcr_positive_hitids}",
        (
            "Exact probe sequence (0 mismatch): "
            f"{exact_probe_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * exact_probe_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length probe alignments (any mismatch count): "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= {(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        f"0 mismatch: {counts['0_MISMATCH']}",
        f"1 mismatch: {counts['1_MISMATCH']}",
        f"2 mismatches: {counts['2_MISMATCHES']}",
        f">=3 mismatches: {counts['GE3_MISMATCHES']}",
        f"Partial/gapped hit only: {counts['PARTIAL_ONLY']}",
        f"No BLAST hit: {counts['NO_HIT']}",
        "",
        (
            "Full-length ungapped BLAST matches among PCR-positive HitIDs: "
            f"{full_length_hitids} / {pcr_positive_hitids} "
            f"= "
            f"{(100.0 * full_length_hitids / pcr_positive_hitids if pcr_positive_hitids else 0.0):.2f} %"
        ),
        (
            "Full-length ungapped probe matches relative to the complete "
            f"validation database: {full_length_hitids} / "
            f"{total_target_sequences} = "
            f"{(100.0 * full_length_hitids / total_target_sequences if total_target_sequences else 0.0):.2f} %"
        ),
        "",
        "Full-length ungapped definition:",
        f"- query starts at 1 and ends at {probe_length};",
        f"- alignment length = {probe_length} nt;",
        "- zero gap openings.",
        "",
        "This version reports BLAST results only.",
        "No probe optimization or IUPAC substitution is proposed.",
        f"Raw BLAST output: {compact_path(raw_file)}",
        f"BLAST query FASTA: {compact_path(query_file)}",
        f"Visual alignments by amplicon: {alignment_amplicon_file}",
        f"Visual alignments by HitID: {alignment_hitid_file}",
        f"Mismatch positions TSV: {mismatch_positions_file}",
        f"Per-variant statistics TSV: {variant_analysis['stats_file']}",
        f"HitID x probe-variant matrix TSV: {variant_analysis['matrix_file']}",
    ]

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"Probe results — {scheme}")

    def _pct(value: int, total: int) -> float:
        return 100.0 * value / total if total else 0.0

    print("Biological coverage — unique PCR-positive targets")
    print(f"  PCR-positive HitIDs   {pcr_positive_hitids}")
    print(
        f"  Exact probe           {counts['0_MISMATCH']:>4}  "
        f"{_pct(counts['0_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  1 mismatch            {counts['1_MISMATCH']:>4}  "
        f"{_pct(counts['1_MISMATCH'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  2 mismatches          {counts['2_MISMATCHES']:>4}  "
        f"{_pct(counts['2_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  >=3 mismatches        {counts['GE3_MISMATCHES']:>4}  "
        f"{_pct(counts['GE3_MISMATCHES'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Partial/gapped        {counts['PARTIAL_ONLY']:>4}  "
        f"{_pct(counts['PARTIAL_ONLY'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  No hit                {counts['NO_HIT']:>4}  "
        f"{_pct(counts['NO_HIT'], pcr_positive_hitids):6.2f} %"
    )
    print(
        f"  Full-length site      {full_length_hitids:>4}  "
        f"{_pct(full_length_hitids, pcr_positive_hitids):6.2f} %"
    )

    print("\nTechnical amplicon statistics")
    print(f"  Valid amplicons       {total_valid_amplicons}")
    print(f"  Exact probe           {amplicon_counts['0_MISMATCH']}")
    print(f"  1 mismatch            {amplicon_counts['1_MISMATCH']}")
    print(f"  2 mismatches          {amplicon_counts['2_MISMATCHES']}")
    print(f"  >=3 mismatches        {amplicon_counts['GE3_MISMATCHES']}")
    print(f"  Partial/gapped        {amplicon_counts['PARTIAL_ONLY']}")
    print(f"  No hit                {amplicon_counts['NO_HIT']}")

    if runtime.VERBOSE:
        print(f"\nRaw BLAST               {raw_file}")
        print(f"Amplicon summary         {amplicon_file}")
        print(f"HitID summary            {hitid_file}")
        print(f"Alignment by amplicon    {alignment_amplicon_file}")
        print(f"Alignment by HitID       {alignment_hitid_file}")
        print(f"Mismatch positions       {mismatch_positions_file}")
        print(f"Variant statistics       {variant_analysis['stats_file']}")
        print(f"HitID x variant matrix   {variant_analysis['matrix_file']}")
        print(f"BLAST summary            {summary_file}")

    return {
        "raw_hsps": len(raw_hits),
        "amplicon_total": total_valid_amplicons,
        "amplicon_exact_0": exact_probe_amplicons,
        "amplicon_full_length": full_length_amplicons,
        "amplicon_1_mismatch": amplicon_counts["1_MISMATCH"],
        "amplicon_2_mismatch": amplicon_counts["2_MISMATCHES"],
        "amplicon_ge3_mismatch": amplicon_counts["GE3_MISMATCHES"],
        "exact_0": counts["0_MISMATCH"],
        "mismatch_1": counts["1_MISMATCH"],
        "mismatch_2": counts["2_MISMATCHES"],
        "mismatch_ge3": counts["GE3_MISMATCHES"],
        "partial_only": counts["PARTIAL_ONLY"],
        "no_hit": counts["NO_HIT"],
        "full_length_hitids": full_length_hitids,
        "pcr_positive_hitids": pcr_positive_hitids,
        "amplicon_summary": amplicon_file,
        "hitid_summary": hitid_file,
        "blast_summary": summary_file,
        "alignment_amplicon_file": alignment_amplicon_file,
        "alignment_hitid_file": alignment_hitid_file,
        "mismatch_positions_file": mismatch_positions_file,
        "variant_statistics_file": variant_analysis["stats_file"],
        "variant_hitid_matrix_file": variant_analysis["matrix_file"],
    }
