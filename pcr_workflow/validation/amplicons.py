# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv
import re

from ..core import runtime
from ..core.paths import compact_path
from ..core.runtime import section

def normalize_primer_identifier(value: str) -> str:
    """
    Normalize an identifier for robust comparisons.

    Examples:
        varVAMP_0      -> varvamp0
        VARVAMP-0      -> varvamp0
        varVAMP_0_LEFT -> varvamp0left
    """
    return re.sub(
        r"[^A-Za-z0-9]+",
        "",
        value,
    ).lower()
def parse_primer_identity(
    primer_name: str,
) -> tuple[str | None, str | None]:
    """
    Extract the VarVAMP scheme and LEFT/RIGHT role from a MFEprimer
    primer identifier.

    This intentionally tolerates harmless decorations added around the
    FASTA identifier, for example:
        varVAMP_0_LEFT
        varVAMP_0_LEFT(+)
        varVAMP_0-LEFT
        varVAMP_0_LEFT_1

    LEFT/RIGHT must still appear as a real token, not as part of another
    ordinary word.
    """
    raw = primer_name.strip().lstrip(">")

    role_match = re.search(
        r"(^|[^A-Za-z0-9])(LEFT|RIGHT)(?=$|[^A-Za-z0-9])",
        raw,
        flags=re.IGNORECASE,
    )

    if not role_match:
        return None, None

    role = role_match.group(2).upper()

    # Everything before the LEFT/RIGHT token is treated as the scheme
    # component. Remove punctuation/separators surrounding the token.
    scheme = raw[:role_match.start(2)]
    scheme = re.sub(
        r"[^A-Za-z0-9]+$",
        "",
        scheme,
    ).strip()

    if not scheme:
        return None, role

    return scheme, role
def same_scheme(
    observed: str | None,
    expected: str,
) -> bool:
    """Compare scheme identifiers while ignoring separators and case."""
    if observed is None:
        return False

    return (
        normalize_primer_identifier(observed)
        == normalize_primer_identifier(expected)
    )
def classify_mfeprimer_product(
    primer_1: str,
    primer_2: str,
    expected_scheme: str,
) -> tuple[str, str, str]:
    """
    Classify a predicted MFEprimer product.

    TARGET:
        one LEFT and one RIGHT primer from the selected VarVAMP scheme.

    SELF_PRIMING:
        LEFT+LEFT or RIGHT+RIGHT from the selected scheme.

    CROSS_SCHEME:
        LEFT/RIGHT roles are recognized, but the primers do not belong to
        the same scheme.

    OTHER:
        at least one raw primer identifier cannot be interpreted safely.
    """
    scheme_1, role_1 = parse_primer_identity(primer_1)
    scheme_2, role_2 = parse_primer_identity(primer_2)

    if role_1 is None or role_2 is None:
        return (
            "OTHER",
            role_1 or "UNKNOWN",
            role_2 or "UNKNOWN",
        )

    # If roles were recognized but either scheme prefix was not recoverable,
    # do not guess.
    if scheme_1 is None or scheme_2 is None:
        return "OTHER", role_1, role_2

    if normalize_primer_identifier(
        scheme_1
    ) != normalize_primer_identifier(
        scheme_2
    ):
        return "CROSS_SCHEME", role_1, role_2

    if not same_scheme(
        scheme_1,
        expected_scheme,
    ):
        return "OTHER", role_1, role_2

    roles = {role_1, role_2}

    if roles == {"LEFT", "RIGHT"}:
        return "TARGET", role_1, role_2

    if (
        role_1 == role_2
        and role_1 in {"LEFT", "RIGHT"}
    ):
        return "SELF_PRIMING", role_1, role_2

    return "OTHER", role_1, role_2
def parse_mfeprimer3_report(report_file: Path) -> dict[str, object]:
    """
    Parse a legacy MFEprimer 3.x specificity text report.

    The parser deliberately avoids fixed-width column positions.

    It uses two independent structures:
    1. the AmpID/HitID/Size summary table;
    2. FASTA-like blocks beginning with:
           >Amp_N PRIMER1 + PRIMER2 ==> HitID

    The sequence block is considered authoritative for extraction, while
    the summary table is used for consistency checks.
    """
    if not report_file.is_file():
        raise FileNotFoundError(
            f"MFEprimer report not found: {report_file}"
        )

    lines = report_file.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    expected_count: int | None = None

    for line in lines:
        match = re.search(
            r"Descriptions\s+of\s+\[\s*(\d+)\s*\]\s+potential\s+amplicons",
            line,
            flags=re.IGNORECASE,
        )
        if match:
            expected_count = int(match.group(1))
            break

    # Parse the summary table without relying on fixed-width spacing.
    summary_by_ampid: dict[int, dict[str, object]] = {}
    in_table = False

    for line in lines:
        stripped = line.strip()

        if re.match(
            r"^AmpID\s+HitID\s+Size\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            in_table = True
            continue

        if in_table and re.match(
            r"^Amplicon\s+details\b",
            stripped,
            flags=re.IGNORECASE,
        ):
            break

        if not in_table:
            continue

        row = re.match(
            r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+",
            line,
        )
        if not row:
            continue

        amp_id = int(row.group(1))
        hit_id = row.group(2)
        size = int(row.group(3))

        if amp_id in summary_by_ampid:
            raise RuntimeError(
                f"Duplicate AmpID {amp_id} in MFEprimer summary table."
            )

        summary_by_ampid[amp_id] = {
            "hit_id": hit_id,
            "size": size,
        }

    # Parse FASTA-like amplicon blocks.
    header_re = re.compile(
        r"^>Amp_(\d+)\s+(\S+)\s+\+\s+(\S+)\s+==>\s+(\S+)"
    )
    sequence_re = re.compile(
        r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$"
    )

    parsed_blocks: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    sequence_parts: list[str] = []

    def flush_current() -> None:
        nonlocal current, sequence_parts

        if current is None:
            return

        sequence = "".join(sequence_parts).upper()

        if not sequence:
            raise RuntimeError(
                f"No amplicon sequence found for Amp_{current['amp_id']}."
            )

        current["sequence"] = sequence
        current["sequence_length"] = len(sequence)
        parsed_blocks.append(current)

        current = None
        sequence_parts = []

    for raw_line in lines:
        stripped = raw_line.strip()

        header_match = header_re.match(stripped)

        if header_match:
            flush_current()

            current = {
                "amp_id": int(header_match.group(1)),
                "primer_1": header_match.group(2),
                "primer_2": header_match.group(3),
                "hit_id": header_match.group(4),
            }
            sequence_parts = []
            continue

        if current is None:
            continue

        if stripped and sequence_re.fullmatch(stripped):
            sequence_parts.append(stripped)
            continue

        # A blank line or a non-sequence line terminates the FASTA block.
        if sequence_parts:
            flush_current()

    flush_current()

    # If no amplicons were reported, allow a clean empty result.
    if expected_count == 0 and not parsed_blocks:
        return {
            "expected_count": 0,
            "summary_by_ampid": summary_by_ampid,
            "amplicons": [],
        }

    if not parsed_blocks:
        raise RuntimeError(
            "No '>Amp_N ... ==> HitID' amplicon sequence blocks were found "
            "in the MFEprimer report. This parser expects the legacy "
            "MFEprimer 3.x text report."
        )

    amp_ids = [int(item["amp_id"]) for item in parsed_blocks]

    if len(amp_ids) != len(set(amp_ids)):
        raise RuntimeError(
            "Duplicate AmpID values were found in the MFEprimer sequence blocks."
        )

    if (
        expected_count is not None
        and len(parsed_blocks) != expected_count
    ):
        raise RuntimeError(
            "MFEprimer report consistency check failed: "
            f"report announces {expected_count} potential amplicons, "
            f"but {len(parsed_blocks)} sequence blocks were parsed."
        )

    # Cross-check HitID and Size between the table and sequence blocks.
    for item in parsed_blocks:
        amp_id = int(item["amp_id"])
        table_item = summary_by_ampid.get(amp_id)

        if table_item is None:
            raise RuntimeError(
                f"Amp_{amp_id} is present in amplicon details but missing "
                "from the AmpID/HitID/Size table."
            )

        if str(table_item["hit_id"]) != str(item["hit_id"]):
            raise RuntimeError(
                f"HitID mismatch for Amp_{amp_id}: "
                f"table={table_item['hit_id']} vs "
                f"details={item['hit_id']}."
            )

        reported_size = int(table_item["size"])
        sequence_length = int(item["sequence_length"])

        if reported_size != sequence_length:
            raise RuntimeError(
                f"Size mismatch for Amp_{amp_id}: "
                f"table reports {reported_size} bp but extracted sequence "
                f"contains {sequence_length} nt."
            )

        item["size"] = reported_size

    return {
        "expected_count": expected_count,
        "summary_by_ampid": summary_by_ampid,
        "amplicons": parsed_blocks,
    }
def analyse_and_extract_mfeprimer_amplicons(
    report_file: Path,
    scheme: str,
    scheme_dir: Path,
    total_target_sequences: int,
    min_target_size: int | None,
    max_target_size: int | None,
    exclude_self_priming: bool = True,
) -> dict[str, object]:
    """
    Classify MFEprimer products, extract valid target amplicons and calculate
    unique primer-pair target coverage.
    """
    parsed = parse_mfeprimer3_report(report_file)
    amplicons = list(parsed["amplicons"])

    all_rows: list[dict[str, object]] = []
    valid_rows: list[dict[str, object]] = []
    non_target_rows: list[dict[str, object]] = []

    for item in amplicons:
        classification, role_1, role_2 = classify_mfeprimer_product(
            str(item["primer_1"]),
            str(item["primer_2"]),
            scheme,
        )

        size = int(item["size"])

        if min_target_size is None or max_target_size is None:
            size_pass = True
        else:
            size_pass = min_target_size <= size <= max_target_size

        allowed_classes = (
            {"TARGET"}
            if exclude_self_priming
            else {"TARGET", "SELF_PRIMING"}
        )
        coverage_eligible = (
            classification in allowed_classes
            and size_pass
        )

        row = {
            "AmpID": int(item["amp_id"]),
            "HitID": str(item["hit_id"]),
            "Primer1": str(item["primer_1"]),
            "Primer2": str(item["primer_2"]),
            "Role1": role_1,
            "Role2": role_2,
            "ProductClass": classification,
            "Size": size,
            "SizePass": "YES" if size_pass else "NO",
            "CoverageEligible": "YES" if coverage_eligible else "NO",
            "Sequence": str(item["sequence"]).upper(),
        }

        all_rows.append(row)

        if coverage_eligible:
            valid_rows.append(row)

        if classification != "TARGET":
            non_target_rows.append(row)

    # Stable first-occurrence order, without double-counting a HitID.
    unique_positive_hitids = list(
        dict.fromkeys(
            str(row["HitID"])
            for row in valid_rows
        )
    )

    valid_amplicons = len(valid_rows)
    unique_positive_targets = len(unique_positive_hitids)
    duplicate_target_products = (
        valid_amplicons - unique_positive_targets
    )

    if total_target_sequences <= 0:
        raise ValueError(
            "Total target sequence count must be greater than zero."
        )

    coverage_pct = (
        100.0
        * unique_positive_targets
        / total_target_sequences
    )

    # Output files.
    all_tsv = scheme_dir / "amplicons_all.tsv"
    valid_tsv = scheme_dir / "amplicons_valid.tsv"
    valid_fasta = scheme_dir / "amplicons_valid.fasta"
    hitids_file = scheme_dir / "positive_hitids.txt"
    non_target_tsv = scheme_dir / "non_target_products.tsv"
    summary_file = scheme_dir / "coverage_summary.txt"

    tsv_columns = [
        "AmpID",
        "HitID",
        "Primer1",
        "Primer2",
        "Role1",
        "Role2",
        "ProductClass",
        "Size",
        "SizePass",
        "CoverageEligible",
    ]

    def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=tsv_columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    write_rows(all_tsv, all_rows)
    write_rows(valid_tsv, valid_rows)
    write_rows(non_target_tsv, non_target_rows)

    with valid_fasta.open("w", encoding="utf-8") as handle:
        for row in valid_rows:
            pair = f"{row['Role1']}+{row['Role2']}"
            handle.write(
                f">Amp_{row['AmpID']}|HitID={row['HitID']}|"
                f"Scheme={scheme}|Pair={pair}|Size={row['Size']}\n"
            )
            sequence = str(row["Sequence"])
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    hitids_file.write_text(
        "".join(f"{hit_id}\n" for hit_id in unique_positive_hitids),
        encoding="utf-8",
    )

    class_counts: dict[str, int] = {}
    for row in all_rows:
        key = str(row["ProductClass"])
        class_counts[key] = class_counts.get(key, 0) + 1

    unrecognized_examples: list[str] = []
    seen_unrecognized: set[tuple[str, str, str]] = set()

    for row in all_rows:
        product_class = str(row["ProductClass"])

        if product_class not in {"OTHER", "CROSS_SCHEME"}:
            continue

        key = (
            str(row["Primer1"]),
            str(row["Primer2"]),
            product_class,
        )

        if key in seen_unrecognized:
            continue

        seen_unrecognized.add(key)
        unrecognized_examples.append(
            f"{product_class}: "
            f"{row['Primer1']} + {row['Primer2']}"
        )

        if len(unrecognized_examples) >= 10:
            break

    if min_target_size is None or max_target_size is None:
        size_filter_text = "none"
    else:
        size_filter_text = f"{min_target_size}-{max_target_size} bp"

    summary_lines = [
        f"MFEprimer target-coverage summary — {scheme}",
        "=" * (38 + len(scheme)),
        "",
        f"MFEprimer report: {compact_path(report_file)}",
        f"Total validation sequences: {total_target_sequences}",
        f"Potential amplicons parsed: {len(all_rows)}",
        f"TARGET products (LEFT+RIGHT or RIGHT+LEFT): "
        f"{class_counts.get('TARGET', 0)}",
        f"SELF_PRIMING products (LEFT+LEFT or RIGHT+RIGHT): "
        f"{class_counts.get('SELF_PRIMING', 0)}",
        f"CROSS_SCHEME products: {class_counts.get('CROSS_SCHEME', 0)}",
        f"OTHER products: {class_counts.get('OTHER', 0)}",
        f"Target amplicon size filter: {size_filter_text}",
        f"Retained amplicons after filtering: {valid_amplicons}",
        f"Unique positive HitIDs: {unique_positive_targets}",
        f"Additional valid amplicons on already-positive HitIDs: "
        f"{duplicate_target_products}",
        (
            "Primer-pair target coverage: "
            f"{unique_positive_targets} / {total_target_sequences} "
            f"= {coverage_pct:.2f} %"
        ),
        "",
        "Coverage-set definition:",
        (
            "- TARGET only (LEFT+RIGHT / RIGHT+LEFT);"
            if exclude_self_priming
            else "- TARGET plus SELF_PRIMING products;"
        ),
        "- optional target-amplicon size filter passed.",
    ]

    if unrecognized_examples:
        summary_lines.extend(
            [
                "",
                "Examples of products not recognized as TARGET/SELF_PRIMING:",
                *[
                    f"- {example}"
                    for example in unrecognized_examples
                ],
            ]
        )

    summary_file.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    section(f"MFEprimer results — {scheme}")
    print(f"Validation sequences       {total_target_sequences}")
    print(f"Potential products         {len(all_rows)}")
    print(f"Retained coverage products {valid_amplicons}")
    print(f"Self-priming products      {class_counts.get('SELF_PRIMING', 0)}")
    print(f"Unique PCR-positive HitIDs {unique_positive_targets}")
    print(f"PCR coverage               {coverage_pct:.2f} %")
    if class_counts.get("CROSS_SCHEME", 0):
        print(f"Cross-scheme products      {class_counts.get('CROSS_SCHEME', 0)}")
    if class_counts.get("OTHER", 0):
        print(f"Unclassified products      {class_counts.get('OTHER', 0)}")
    if unrecognized_examples and runtime.VERBOSE:
        print("\nUnrecognized primer-name examples:")
        for example in unrecognized_examples[:5]:
            print(f"- {example}")
    if runtime.VERBOSE:
        print(f"Valid amplicons TSV        {valid_tsv}")
        print(f"Valid amplicons FASTA      {valid_fasta}")
        print(f"Positive HitIDs            {hitids_file}")
        print(f"Non-target products        {non_target_tsv}")
        print(f"Coverage summary           {summary_file}")

    return {
        "scheme": scheme,
        "potential_amplicons": len(all_rows),
        "target_products": class_counts.get("TARGET", 0),
        "self_priming_products": class_counts.get("SELF_PRIMING", 0),
        "valid_target_amplicons": valid_amplicons,
        "unique_positive_hitids": unique_positive_targets,
        "coverage_pct": coverage_pct,
        "valid_tsv": valid_tsv,
        "valid_fasta": valid_fasta,
        "positive_hitids": hitids_file,
        "summary": summary_file,
        "valid_rows": valid_rows,
        "self_priming_excluded": exclude_self_priming,
    }
