# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

IUPAC_DNA = set("ACGTRYSWKMBDHVN")
PROBE_IUPAC_TO_BASES = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}
