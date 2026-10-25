# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from ..core.runtime import section


import matplotlib
matplotlib.use("Agg")

from PIL import Image
from io import BytesIO
from pathlib import Path
import fitz
import matplotlib.pyplot as plt
import pandas as pd
import shutil

from ..core.paths import compact_path

def copy_varvamp_tree(varvamp_dir: Path, results_dir: Path, mode: str) -> Path:
    """Copy the complete successful VarVAMP tree into this isolated run."""

    destination = results_dir / "varvamp"
    if destination.exists():
        raise RuntimeError(
            "The VarVAMP result directory already exists inside this run: "
            f"{destination}. Refusing to overwrite an existing result."
        )
    shutil.copytree(varvamp_dir, destination)
    (destination / "assay_mode.txt").write_text(
        mode.upper() + "\n",
        encoding="utf-8",
    )
    return destination
def plot_bed_segments(
    bed_file: Path,
    output_file: Path,
    title: str,
    y_label: str,
    line_width: float,
) -> None:
    """Create a simple overview plot from a BED file."""

    if not bed_file.is_file():
        return

    try:
        data = pd.read_csv(bed_file, sep="\t", header=None, comment="#")
    except pd.errors.EmptyDataError:
        return

    if data.empty or data.shape[1] < 3:
        return

    plt.figure(figsize=(12, 4))
    for index, row in data.iterrows():
        plt.plot(
            [row.iloc[1], row.iloc[2]],
            [index, index],
            linewidth=line_width,
        )
    plt.xlabel("Alignment / consensus position")
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
def save_pdf_first_page_as_png(
    pdf_file: Path,
    output_file: Path,
    *,
    zoom: int = 4,
) -> None:
    """Convert the first page of a PDF to PNG."""

    if not pdf_file.is_file():
        return
    with fitz.open(pdf_file) as pdf:
        if len(pdf) == 0:
            return
        page = pdf[0]
        pixmap = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    image = Image.open(BytesIO(pixmap.tobytes("png"))).copy()
    image.save(output_file)
def print_varvamp_tables(varvamp_results_dir: Path, mode: str) -> None:
    """
    Display only the principal VarVAMP result files in the terminal.

    Large FASTA alignments, consensus sequences, BED files, logs, and other
    technical outputs remain saved in the VarVAMP result directory but are not
    printed to the terminal.
    """

    principal_files_by_mode = {
        "qpcr": [
            "oligos.fasta",
            "qpcr_design.tsv",
            "qpcr_primers.tsv",
        ],
        "single": [
            "primer.tsv",
            "primer_to_amplicon_assignments.tabular",
        ],
        "tiled": [
            "primer.tsv",
            "primer_to_amplicon_assignments.tabular",
        ],
    }

    principal_files = principal_files_by_mode.get(mode, [])

    section(f"VarVAMP {mode.upper()} principal results")

    displayed = 0

    for filename in principal_files:
        matches = sorted(
            path
            for path in varvamp_results_dir.rglob(filename)
            if path.is_file()
        )

        for path in matches:
            relative = path.relative_to(varvamp_results_dir)
            title = str(relative)

            print(f"\n{title}")
            print("-" * len(title))

            try:
                content = path.read_text(
                    encoding="utf-8",
                    errors="replace",
                ).rstrip()
            except OSError as error:
                print(f"[Could not read file: {error}]")
                continue

            if content:
                print(content)
            else:
                print("[empty file]")

            displayed += 1

    if displayed == 0:
        print("No principal VarVAMP result files were found.")

    print(
        f"\nComplete VarVAMP result directory: "
        f"{compact_path(varvamp_results_dir)}"
    )
def create_varvamp_visual_summaries(varvamp_results_dir: Path) -> None:
    """Create optional PNG summaries from standard VarVAMP outputs."""

    plot_bed_segments(
        varvamp_results_dir / "amplicons.bed",
        varvamp_results_dir / "amplicons_overview.png",
        "Amplicon positions",
        "Amplicon",
        6,
    )
    plot_bed_segments(
        varvamp_results_dir / "primers.bed",
        varvamp_results_dir / "primers_overview.png",
        "Primer / probe positions",
        "Oligonucleotide",
        4,
    )
    save_pdf_first_page_as_png(
        varvamp_results_dir / "amplicon_plot.pdf",
        varvamp_results_dir / "amplicon_plot.png",
    )
    save_pdf_first_page_as_png(
        varvamp_results_dir / "per_base_mismatches.pdf",
        varvamp_results_dir / "per_base_mismatches.png",
    )
