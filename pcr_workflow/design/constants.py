# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

VALID_BASES = {"A", "C", "G", "T"}
GAP_CHARACTERS = {"-", "."}
ALLOWED_NUCLEOTIDES = set("ACGTURYSWKMBDHVN-.")
MAFFT_STRATEGIES: dict[str, dict[str, object]] = {
    "auto": {
        "description": (
            "Auto — MAFFT automatically selects an appropriate strategy "
            "according to dataset size."
        ),
        "arguments": ["--auto"],
    },
    "fft-ns-1": {
        "description": (
            "FFT-NS-1 — very fast progressive alignment; useful for very "
            "large datasets."
        ),
        "arguments": ["--retree", "1", "--maxiterate", "0"],
    },
    "fft-ns-2": {
        "description": (
            "FFT-NS-2 — fast progressive alignment with two guide-tree "
            "calculations."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "0"],
    },
    "fft-ns-i-2": {
        "description": "FFT-NS-i (2 cycles) — fast iterative refinement.",
        "arguments": ["--retree", "2", "--maxiterate", "2"],
    },
    "fft-ns-i-1000": {
        "description": (
            "FFT-NS-i (up to 1000 cycles) — intensive iterative refinement "
            "while remaining more scalable than pairwise strategies."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "1000"],
    },
    "nw-ns-2": {
        "description": "NW-NS-2 — progressive alignment without FFT approximation.",
        "arguments": ["--retree", "2", "--maxiterate", "0", "--nofft"],
    },
    "nw-ns-i-2": {
        "description": "NW-NS-i (2 cycles) — iterative refinement without FFT.",
        "arguments": ["--retree", "2", "--maxiterate", "2", "--nofft"],
    },
    "nw-ns-i-1000": {
        "description": (
            "NW-NS-i (up to 1000 cycles) — intensive iterative refinement "
            "without FFT."
        ),
        "arguments": ["--retree", "2", "--maxiterate", "1000", "--nofft"],
    },
    "l-ins-i": {
        "description": (
            "L-INS-i — high accuracy for locally alignable regions with "
            "variable flanks; generally for relatively small datasets."
        ),
        "arguments": ["--localpair", "--maxiterate", "1000"],
    },
    "g-ins-i": {
        "description": (
            "G-INS-i — high accuracy for globally alignable sequences of "
            "similar length; generally for relatively small datasets."
        ),
        "arguments": ["--globalpair", "--maxiterate", "1000"],
    },
    "e-ins-i": {
        "description": (
            "E-INS-i — high accuracy when conserved regions are separated by "
            "large difficult-to-align regions."
        ),
        "arguments": ["--ep", "0", "--genafpair", "--maxiterate", "1000"],
    },
    "parttree": {
        "description": (
            "NW-NS-PartTree-1 — scalable strategy intended for extremely "
            "large datasets."
        ),
        "arguments": [
            "--retree",
            "1",
            "--maxiterate",
            "0",
            "--nofft",
            "--parttree",
        ],
    },
}
VARVAMP_ADVANCED_CONFIG_DEFAULTS: dict[str, dict[str, object]] = {
    "single": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
    },
    "tiled": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
    },
    "qpcr": {
        "PRIMER_SIZES": (18, 25, 21),
        "PRIMER_TMP": (56, 63, 60),
        "PRIMER_GC_RANGE": (30, 75, 50),
        "QPROBE_SIZES": (18, 30, 25),
        "QPROBE_TMP": (60, 72, 67),
        "QPROBE_GC_RANGE": (35, 85, 60),
        "QPRIMER_DIFF": 3,
        "QPROBE_TEMP_DIFF": (2, 12),
        "QPROBE_DISTANCE": (4, 25),
        "QAMPLICON_LENGTH": (70, 200),
        "QAMPLICON_GC": (40, 60),
    },
}
VARVAMP_ADVANCED_CONFIG_ORDER: dict[str, tuple[str, ...]] = {
    "single": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
    ),
    "tiled": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
    ),
    "qpcr": (
        "PRIMER_SIZES",
        "PRIMER_TMP",
        "PRIMER_GC_RANGE",
        "QPROBE_SIZES",
        "QPROBE_TMP",
        "QPROBE_GC_RANGE",
        "QPRIMER_DIFF",
        "QPROBE_TEMP_DIFF",
        "QPROBE_DISTANCE",
        "QAMPLICON_LENGTH",
        "QAMPLICON_GC",
    ),
}
