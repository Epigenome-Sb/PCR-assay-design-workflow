# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from pathlib import Path

from ..core import runtime
from .amplicons import classify_mfeprimer_product, parse_mfeprimer3_report
from .primer_qc import compact_primer_label, parse_mfeprimer_primer_properties

def print_compact_primer_table(report_file: Path, scheme: str) -> None:
    """Display a terminal-width friendly MFEprimer primer table."""
    rows = parse_mfeprimer_primer_properties(report_file)
    if not rows:
        if runtime.VERBOSE:
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
