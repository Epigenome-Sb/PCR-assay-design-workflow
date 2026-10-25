# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path
import csv
import re

def parse_mfeprimer_secondary_structures(report_file: Path) -> dict[str, object]:
    """
    Parse Hairpin and Dimer sections from a legacy MFEprimer 3.x text report.

    Reported counts are authoritative. Detailed structure blocks are parsed
    when their textual format is recognized. If MFEprimer reports structures
    but their detail format is not recognized, the counts are still preserved
    and the raw report remains the authoritative source.
    """
    if not report_file.is_file():
        raise FileNotFoundError(f"MFEprimer report not found: {report_file}")

    lines = report_file.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count_patterns = {
        "hairpin": re.compile(r"^\s*Hairpin\s+List\s*\(\s*(\d+)\s*\)\s*$", re.IGNORECASE),
        "dimer": re.compile(r"^\s*Dimer\s+List\s*\(\s*(\d+)\s*\)\s*$", re.IGNORECASE),
    }

    counts: dict[str, int | None] = {"hairpin": None, "dimer": None}
    section_start: dict[str, int | None] = {"hairpin": None, "dimer": None}

    for index, line in enumerate(lines):
        for kind, pattern in count_patterns.items():
            match = pattern.match(line)
            if match:
                counts[kind] = int(match.group(1))
                section_start[kind] = index + 1

    # MFEprimer 3.x reports normally contain both sections. Keep "not found"
    # distinct from a true zero so the workflow does not silently invent data.
    hairpin_count = counts["hairpin"]
    dimer_count = counts["dimer"]

    def section_lines(kind: str) -> list[str]:
        start = section_start[kind]
        if start is None:
            return []

        end = len(lines)
        if kind == "hairpin" and section_start["dimer"] is not None:
            end = int(section_start["dimer"]) - 1
        elif kind == "dimer":
            for idx in range(start, len(lines)):
                if re.search(
                    r"Descriptions\s+of\s+\[\s*\d+\s*\]\s+potential\s+amplicons",
                    lines[idx],
                    flags=re.IGNORECASE,
                ):
                    end = idx
                    break
        return lines[start:end]

    def parse_numeric(block: list[str], patterns: list[re.Pattern[str]]) -> float | None:
        joined = "\n".join(block)
        for pattern in patterns:
            match = pattern.search(joined)
            if match:
                try:
                    return float(match.group(1))
                except ValueError:
                    return None
        return None

    score_patterns = [
        re.compile(r"\bScore\s*[:=]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE),
    ]
    tm_patterns = [
        re.compile(r"\bTm\s*[:=]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE),
        re.compile(
            r"\bmelting\s+temperature\s*[:=]\s*(-?\d+(?:\.\d+)?)",
            re.IGNORECASE,
        ),
    ]
    dg_patterns = [
        re.compile(
            r"(?:Delta\s*G|Dg|ΔG)\s*[:=]\s*(-?\d+(?:\.\d+)?)",
            re.IGNORECASE,
        ),
    ]

    def parse_items(kind: str) -> list[dict[str, object]]:
        block_lines = section_lines(kind)
        if not block_lines:
            return []

        label = "Hairpin" if kind == "hairpin" else "Dimer"
        item_re = re.compile(
            rf"^\s*{label}\s+(\d+)\s*:\s*(.*?)\s*$",
            re.IGNORECASE,
        )

        starts: list[tuple[int, re.Match[str]]] = []
        for idx, raw in enumerate(block_lines):
            match = item_re.match(raw)
            if match:
                starts.append((idx, match))

        items: list[dict[str, object]] = []
        for pos, (start_idx, match) in enumerate(starts):
            end_idx = starts[pos + 1][0] if pos + 1 < len(starts) else len(block_lines)
            item_block = block_lines[start_idx:end_idx]
            subject = match.group(2).strip()

            oligo_1 = subject
            oligo_2 = ""
            if kind == "dimer":
                split = re.split(r"\s+(?:x|×)\s+", subject, maxsplit=1, flags=re.IGNORECASE)
                if len(split) == 2:
                    oligo_1, oligo_2 = split[0].strip(), split[1].strip()

            items.append(
                {
                    "index": int(match.group(1)),
                    "subject": subject,
                    "oligo_1": oligo_1,
                    "oligo_2": oligo_2,
                    "score": parse_numeric(item_block, score_patterns),
                    "tm_c": parse_numeric(item_block, tm_patterns),
                    "delta_g": parse_numeric(item_block, dg_patterns),
                    "raw_block": "\n".join(item_block).strip(),
                }
            )
        return items

    hairpins = parse_items("hairpin")
    dimers = parse_items("dimer")

    return {
        "hairpin_section_found": hairpin_count is not None,
        "dimer_section_found": dimer_count is not None,
        "hairpin_count": hairpin_count,
        "dimer_count": dimer_count,
        "hairpins": hairpins,
        "dimers": dimers,
    }
def print_secondary_structure_qc(
    secondary_by_scheme: dict[str, dict[str, object]],
    selected_schemes: list[str],
    assay_mode: str,
) -> None:
    """Display MFEprimer hairpin/dimer data only when those report sections exist."""
    available = [
        scheme
        for scheme in selected_schemes
        if secondary_by_scheme[scheme]["hairpin_count"] is not None
        or secondary_by_scheme[scheme]["dimer_count"] is not None
    ]

    print("\nSecondary-structure QC — MFEprimer")
    print("----------------------------------")

    if not available:
        print(
            "Hairpin/dimer sections are not present in these MFEprimer specificity "
            "reports, so no secondary-structure conclusion is drawn here."
        )
        if assay_mode == "qpcr":
            print(
                "Only LEFT/RIGHT primers are present in the MFEprimer QC input. "
                "Probe hairpins/dimers are therefore not assessed here; the probe is "
                "evaluated separately by BLAST in the probe-validation step."
            )
        return

    print(f"{'Scheme':<16} {'Hairpins':>10} {'Dimers':>10}")
    print("-" * 38)

    for scheme in selected_schemes:
        result = secondary_by_scheme[scheme]
        hairpins = (
            str(result["hairpin_count"])
            if result["hairpin_count"] is not None
            else "N/A"
        )
        dimers = (
            str(result["dimer_count"])
            if result["dimer_count"] is not None
            else "N/A"
        )
        print(f"{scheme:<16} {hairpins:>10} {dimers:>10}")

    print(
        "\nCounts mean structures reported by MFEprimer under the selected "
        "analysis parameters; zero does not prove that a structure can never form."
    )

    positive_found = False
    for scheme in selected_schemes:
        result = secondary_by_scheme[scheme]
        reported_positive = (
            (result["hairpin_count"] or 0) > 0
            or (result["dimer_count"] or 0) > 0
        )
        if not reported_positive:
            continue

        positive_found = True
        print(f"\nDetected structures — {scheme}")
        print("-" * (23 + len(scheme)))

        details = list(result["hairpins"]) + list(result["dimers"])
        if not details:
            print(
                "MFEprimer reported one or more structures, but the detailed "
                "block could not be parsed. See the raw MFEprimer report."
            )
            continue

        print(
            f"{'Type':<8} {'Oligo(s)':<31} {'Score':>7} "
            f"{'Tm°C':>8} {'ΔG':>9}"
        )
        print("-" * 68)

        for kind, items in (
            ("Hairpin", result["hairpins"]),
            ("Dimer", result["dimers"]),
        ):
            for item in items:
                if kind == "Dimer" and item["oligo_2"]:
                    oligos = f"{item['oligo_1']} × {item['oligo_2']}"
                else:
                    oligos = str(item["oligo_1"])

                score = (
                    f"{float(item['score']):.2f}"
                    if item["score"] is not None
                    else "N/A"
                )
                tm_c = (
                    f"{float(item['tm_c']):.2f}"
                    if item["tm_c"] is not None
                    else "N/A"
                )
                delta_g = (
                    f"{float(item['delta_g']):.2f}"
                    if item["delta_g"] is not None
                    else "N/A"
                )
                print(
                    f"{kind:<8} {oligos[:31]:<31} "
                    f"{score:>7} {tm_c:>8} {delta_g:>9}"
                )

    if not positive_found:
        print("\nNo hairpins or dimers were reported in the available QC sections.")
def save_secondary_structure_qc(
    coverage_root: Path,
    secondary_by_scheme: dict[str, dict[str, object]],
    selected_schemes: list[str],
) -> tuple[Path, Path]:
    """Save MFEprimer hairpin/dimer summary and any parsed detailed structures."""
    summary_file = coverage_root / "secondary_structure_summary.tsv"
    details_file = coverage_root / "secondary_structure_details.tsv"

    with summary_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "hairpin_section_found",
                "hairpins_reported",
                "hairpin_details_parsed",
                "dimer_section_found",
                "dimers_reported",
                "dimer_details_parsed",
            ]
        )
        for scheme in selected_schemes:
            result = secondary_by_scheme[scheme]
            writer.writerow(
                [
                    scheme,
                    result["hairpin_section_found"],
                    "" if result["hairpin_count"] is None else result["hairpin_count"],
                    len(result["hairpins"]),
                    result["dimer_section_found"],
                    "" if result["dimer_count"] is None else result["dimer_count"],
                    len(result["dimers"]),
                ]
            )

    with details_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scheme",
                "structure_type",
                "structure_index",
                "oligo_1",
                "oligo_2",
                "score",
                "tm_c",
                "delta_g",
                "raw_block",
            ]
        )
        for scheme in selected_schemes:
            result = secondary_by_scheme[scheme]
            for kind, items in (
                ("HAIRPIN", result["hairpins"]),
                ("DIMER", result["dimers"]),
            ):
                for item in items:
                    writer.writerow(
                        [
                            scheme,
                            kind,
                            item["index"],
                            item["oligo_1"],
                            item["oligo_2"],
                            "" if item["score"] is None else item["score"],
                            "" if item["tm_c"] is None else item["tm_c"],
                            "" if item["delta_g"] is None else item["delta_g"],
                            item["raw_block"],
                        ]
                    )

    return summary_file, details_file
def parse_mfeprimer_primer_properties(report_file: Path) -> list[dict[str, object]]:
    """Parse the compact primer-property rows from a legacy MFEprimer report."""
    row_re = re.compile(
        r"^\s*(\S+)\s+([ACGTRYSWKMBDHVNacgtryswkmbdhvn]+)\s+"
        r"(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+(-?[0-9.]+)\s+"
        r"(\d+)\s+(\d+)\s*$"
    )
    rows: list[dict[str, object]] = []
    for raw in report_file.read_text(encoding="utf-8", errors="replace").splitlines():
        match = row_re.match(raw)
        if not match:
            continue
        rows.append({
            "id": match.group(1),
            "sequence": match.group(2).upper(),
            "length": int(match.group(3)),
            "gc": float(match.group(4)),
            "tm": float(match.group(5)),
            "dg": float(match.group(6)),
            "bind_plus": int(match.group(7)),
            "bind_minus": int(match.group(8)),
        })
    return rows
def compact_primer_label(primer_id: str, scheme: str) -> str:
    """Remove a repeated scheme prefix while preserving LEFT/RIGHT variant labels."""
    label = primer_id
    prefix_re = re.compile(rf"^{re.escape(scheme)}[_-]?", flags=re.IGNORECASE)
    label = prefix_re.sub("", label)
    return label or primer_id
