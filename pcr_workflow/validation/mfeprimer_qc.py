# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""MFEprimer hairpin/dimer QC, primer properties and product overviews."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from . import run_management as _run_management
from .mfeprimer import classify_mfeprimer_product, parse_mfeprimer3_report


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
    coverage_root.mkdir(parents=True, exist_ok=True)
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



def print_compact_primer_table(report_file: Path, scheme: str) -> None:
    """Display a terminal-width friendly MFEprimer primer table."""
    rows = parse_mfeprimer_primer_properties(report_file)
    if not rows:
        if _run_management.VERBOSE:
            print(f"No primer-property rows parsed from {report_file.name}.")
        return

    print(f"\nConcrete primer variants — {scheme}")
    print("-" * (27 + len(scheme)))
    print(
        f"{'Primer':<12} {'Sequence':<30} {'Len':>3} "
        f"{'GC%':>6} {'Tm°C':>7} {'ΔG':>8} {'Bind+':>7} {'Bind-':>7}"
    )
    print("-" * 86)
    for row in rows:
        label = compact_primer_label(str(row['id']), scheme)
        print(
            f"{label:<12} {str(row['sequence']):<30} {int(row['length']):>3} "
            f"{float(row['gc']):>6.2f} {float(row['tm']):>7.2f} "
            f"{float(row['dg']):>8.2f} {int(row['bind_plus']):>7} "
            f"{int(row['bind_minus']):>7}"
        )



def inspect_mfeprimer_products(report_file: Path, scheme: str) -> dict[str, object]:
    """
    Inspect ALL MFEprimer potential products before any workflow filtering.

    The returned structure keeps:
    - total potential amplicons;
    - product classes;
    - unique HitIDs represented in each class;
    - exact concrete primer-pair combinations.
    """
    parsed = parse_mfeprimer3_report(report_file)

    classes = ("TARGET", "SELF_PRIMING", "CROSS_SCHEME", "OTHER")
    class_counts = {key: 0 for key in classes}
    class_hitids: dict[str, set[str]] = {key: set() for key in classes}

    pair_counts = {
        key: 0
        for key in ("LEFT×RIGHT", "RIGHT×LEFT", "LEFT×LEFT", "RIGHT×RIGHT")
    }

    all_hitids: set[str] = set()
    concrete: dict[tuple[str, str, str], dict[str, object]] = {}

    for item in parsed["amplicons"]:
        primer_1 = str(item["primer_1"])
        primer_2 = str(item["primer_2"])
        hit_id = str(item["hit_id"])

        product_class, role_1, role_2 = classify_mfeprimer_product(
            primer_1,
            primer_2,
            scheme,
        )

        class_counts[product_class] = class_counts.get(product_class, 0) + 1
        class_hitids.setdefault(product_class, set()).add(hit_id)
        all_hitids.add(hit_id)

        role_pair = f"{role_1}×{role_2}"
        if role_pair in pair_counts:
            pair_counts[role_pair] += 1

        label_1 = compact_primer_label(primer_1, scheme)
        label_2 = compact_primer_label(primer_2, scheme)
        key = (label_1, label_2, product_class)

        if key not in concrete:
            concrete[key] = {
                "primer_1": label_1,
                "primer_2": label_2,
                "class": product_class,
                "amplicons": 0,
                "hitids": set(),
            }

        concrete[key]["amplicons"] = int(concrete[key]["amplicons"]) + 1
        concrete[key]["hitids"].add(hit_id)

    concrete_rows: list[dict[str, object]] = []
    for item in concrete.values():
        concrete_rows.append(
            {
                "primer_1": item["primer_1"],
                "primer_2": item["primer_2"],
                "class": item["class"],
                "amplicons": int(item["amplicons"]),
                "unique_hitids": len(item["hitids"]),
            }
        )

    concrete_rows.sort(
        key=lambda row: (
            0 if row["class"] == "TARGET" else
            1 if row["class"] == "SELF_PRIMING" else 2,
            str(row["primer_1"]),
            str(row["primer_2"]),
        )
    )

    return {
        "potential": len(parsed["amplicons"]),
        "all_unique_hitids": len(all_hitids),
        "class_counts": class_counts,
        "class_unique_hitids": {
            key: len(class_hitids.get(key, set()))
            for key in classes
        },
        "pair_counts": pair_counts,
        "concrete_pairs": concrete_rows,
    }



def print_general_potential_amplicon_statistics(
    previews: dict[str, dict[str, object]],
    selected_schemes: list[str],
    total_target_sequences: int,
) -> None:
    """
    Display the complete MFEprimer product picture BEFORE any workflow filter.
    """
    print("\nMFEprimer potential-product overview — BEFORE filtering")
    print("-------------------------------------------------------")
    print(
        "These are all potential products reported by MFEprimer. "
        "No TARGET/SELF_PRIMING exclusion has been applied yet."
    )

    for scheme in selected_schemes:
        preview = previews[scheme]
        counts = preview["class_counts"]
        hitids = preview["class_unique_hitids"]

        any_hitids = int(preview["all_unique_hitids"])
        any_pct = (
            100.0 * any_hitids / total_target_sequences
            if total_target_sequences else 0.0
        )

        print(f"\n{scheme}")
        print("=" * len(scheme))
        print(f"Validation database                 : {total_target_sequences} sequences")
        print(f"Potential amplicons reported        : {int(preview['potential'])}")
        print(
            f"Sequences with >=1 potential product: "
            f"{any_hitids}/{total_target_sequences} ({any_pct:.2f} %)"
        )

        print("\nProduct classes")
        print("---------------")
        print(
            f"{'Class':<18} {'Amplicons':>10} "
            f"{'Unique HitIDs':>14} {'% database':>12}"
        )
        print("-" * 58)

        for product_class, label in (
            ("TARGET", "TARGET"),
            ("SELF_PRIMING", "SELF_PRIMING"),
            ("CROSS_SCHEME", "CROSS_SCHEME"),
            ("OTHER", "OTHER"),
        ):
            n_amp = int(counts.get(product_class, 0))
            n_hits = int(hitids.get(product_class, 0))
            pct = (
                100.0 * n_hits / total_target_sequences
                if total_target_sequences else 0.0
            )
            print(
                f"{label:<18} {n_amp:>10} "
                f"{n_hits:>14} {pct:>11.2f}%"
            )

        print(
            "\nTARGET = LEFT×RIGHT or RIGHT×LEFT. "
            "SELF_PRIMING = LEFT×LEFT or RIGHT×RIGHT."
        )
        print(
            "Note: class-specific unique HitID counts are not additive; the same "
            "sequence may generate products in more than one class."
        )



def print_potential_amplicons_by_pair(
    previews: dict[str, dict[str, object]],
    selected_schemes: list[str],
) -> None:
    """Show potential MFEprimer amplicons at the exact primer-pair level."""
    print("\nPotential amplicons by primer pair")
    print("----------------------------------")
    print(
        "Counts below are raw MFEprimer potential products before the workflow "
        "classification/filtering step."
    )

    for scheme in selected_schemes:
        preview = previews[scheme]
        rows = list(preview.get("concrete_pairs", []))

        print(f"\n{scheme}")
        print("-" * len(scheme))
        print(
            f"{'Primer pair':<31} {'Class':<14} "
            f"{'Amplicons':>10} {'HitIDs':>8}"
        )
        print("-" * 67)

        if not rows:
            print("(no potential amplicons)")
            continue

        for row in rows:
            pair_label = f"{row['primer_1']} × {row['primer_2']}"
            print(
                f"{pair_label:<31} {str(row['class']):<14} "
                f"{int(row['amplicons']):>10} "
                f"{int(row['unique_hitids']):>8}"
            )

        print(
            f"{'TOTAL':<31} {'':<14} "
            f"{int(preview['potential']):>10}"
        )

