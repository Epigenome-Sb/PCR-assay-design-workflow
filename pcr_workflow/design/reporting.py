# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham
"""Generate visual summaries from VarVAMP results."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

import matplotlib

# Allow execution on servers without a display.
matplotlib.use("Agg")

import fitz
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image


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
        data = pd.read_csv(
            bed_file,
            sep="\t",
            header=None,
            comment="#",
        )
    except pd.errors.EmptyDataError:
        return

    if data.empty or data.shape[1] < 3:
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)

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

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with fitz.open(pdf_file) as pdf:
        if len(pdf) == 0:
            return
        page = pdf[0]
        pixmap = page.get_pixmap(
            matrix=fitz.Matrix(zoom, zoom),
            alpha=False,
        )

    image = Image.open(BytesIO(pixmap.tobytes("png"))).copy()
    image.save(output_file)
    image.close()


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