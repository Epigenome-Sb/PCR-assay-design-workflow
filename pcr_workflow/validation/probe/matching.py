# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

def classify_full_length_probe_hit(
    hit: dict[str, object],
    probe_length: int,
) -> tuple[bool, str]:
    """
    Define the transparent BLAST classes used by this workflow.

    Full-length ungapped means:
    - alignment starts at query base 1;
    - alignment ends at the complete probe length;
    - aligned length equals probe length;
    - zero gap openings.

    No biological PASS/FAIL label is imposed.
    """
    full_length = (
        int(hit["qstart"]) == 1
        and int(hit["qend"]) == probe_length
        and int(hit["length"]) == probe_length
        and int(hit["gapopen"]) == 0
    )

    if not full_length:
        return False, "PARTIAL_OR_GAPPED"

    mismatches = int(hit["mismatch"])

    if mismatches == 0:
        return True, "0_MISMATCH"

    if mismatches == 1:
        return True, "1_MISMATCH"

    if mismatches == 2:
        return True, "2_MISMATCHES"

    return True, "GE3_MISMATCHES"
def build_blast_match_line(
    query_aligned: str,
    subject_aligned: str,
) -> str:
    """
    Create a BLAST-like visual match line.

    |  exact nucleotide identity
    .  aligned non-identical bases
       gap
    """
    symbols: list[str] = []

    for query_base, subject_base in zip(
        query_aligned.upper(),
        subject_aligned.upper(),
    ):
        if query_base == "-" or subject_base == "-":
            symbols.append(" ")
        elif query_base == subject_base:
            symbols.append("|")
        else:
            symbols.append(".")

    return "".join(symbols)
def extract_full_length_mismatches(
    row: dict[str, object],
) -> list[dict[str, object]]:
    """
    Extract mismatch positions from a full-length ungapped BLAST alignment.

    ProbePosition is 1-based relative to the concrete BLAST query variant.
    SubjectCoordinate respects the BLAST subject orientation.
    """
    if str(row["FullLengthUngapped"]) != "YES":
        return []

    query = str(row["QueryAligned"]).upper()
    subject = str(row["SubjectAligned"]).upper()

    if len(query) != len(subject):
        return []

    subject_start = int(row["SubjectStart"])
    subject_end = int(row["SubjectEnd"])
    subject_step = 1 if subject_end >= subject_start else -1
    subject_coord = subject_start

    mismatch_rows: list[dict[str, object]] = []

    for probe_position, (query_base, subject_base) in enumerate(
        zip(query, subject),
        start=1,
    ):
        if query_base == "-" or subject_base == "-":
            if subject_base != "-":
                subject_coord += subject_step
            continue

        current_subject_coord = subject_coord

        if query_base != subject_base:
            mismatch_rows.append(
                {
                    "AmpID": int(row["AmpID"]),
                    "HitID": str(row["HitID"]),
                    "BlastCategory": str(row["BlastCategory"]),
                    "QueryVariant": str(row["QueryVariant"]),
                    "ProbePosition": probe_position,
                    "ProbeBase": query_base,
                    "TargetBase": subject_base,
                    "SubjectCoordinate": current_subject_coord,
                    "SubjectStrand": (
                        "+"
                        if subject_step == 1
                        else "-"
                    ),
                }
            )

        subject_coord += subject_step

    return mismatch_rows
def format_probe_blast_alignment_block(
    row: dict[str, object],
) -> str:
    """Create one human-readable BLAST alignment block."""
    category = str(row["BlastCategory"])
    hit_id = str(row["HitID"])
    amp_id = int(row["AmpID"])
    amp_size = int(row["AmpliconSize"])

    lines = [
        "=" * 78,
        f"HitID: {hit_id} | AmpID: {amp_id} | Amplicon size: {amp_size} bp",
        f"Category: {category}",
    ]

    if category == "NO_HIT":
        lines.extend(
            [
                "No BLAST HSP was reported for this amplicon.",
                "",
            ]
        )
        return "\n".join(lines)

    qseq = str(row["QueryAligned"])
    sseq = str(row["SubjectAligned"])
    match_line = build_blast_match_line(qseq, sseq)

    sstart = int(row["SubjectStart"])
    send = int(row["SubjectEnd"])
    strand = "+" if send >= sstart else "-"

    lines.extend(
        [
            (
                f"Identity: {row['PercentIdentity']}% | "
                f"Alignment: {row['AlignmentLength']} nt | "
                f"Mismatches: {row['Mismatches']} | "
                f"Gap openings: {row['GapOpenings']}"
            ),
            (
                f"Query variant: {row['QueryVariant']} | "
                f"Query coordinates: {row['QueryStart']}..{row['QueryEnd']}"
            ),
            (
                f"Subject coordinates: {sstart}..{send} | "
                f"Subject strand: {strand}"
            ),
            f"E-value: {row['Evalue']} | Bit score: {row['BitScore']}",
            "",
            f"Probe  {str(row['QueryStart']).rjust(4)}  {qseq}  {row['QueryEnd']}",
            f"             {match_line}",
            f"Target {str(sstart).rjust(4)}  {sseq}  {send}",
        ]
    )

    mismatch_rows = extract_full_length_mismatches(row)

    if str(row["FullLengthUngapped"]) == "YES":
        if mismatch_rows:
            mismatch_text = ", ".join(
                (
                    f"pos{item['ProbePosition']} "
                    f"{item['ProbeBase']}>{item['TargetBase']} "
                    f"(target:{item['SubjectCoordinate']})"
                )
                for item in mismatch_rows
            )
            lines.extend(
                [
                    "",
                    f"Mismatch positions: {mismatch_text}",
                ]
            )
        else:
            lines.extend(
                [
                    "",
                    "Mismatch positions: none",
                ]
            )
    else:
        lines.extend(
            [
                "",
                (
                    "Note: this is not a full-length ungapped alignment; "
                    "mismatch-position interpretation is not used."
                ),
            ]
        )

    lines.append("")
    return "\n".join(lines)
