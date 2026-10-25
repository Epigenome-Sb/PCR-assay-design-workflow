# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Aymane Faham

from __future__ import annotations

from ..constants import PROBE_IUPAC_TO_BASES

def expand_iupac_probe(
    sequence: str,
    max_variants: int = 256,
) -> list[str]:
    """
    Expand an already-degenerate VarVAMP probe into concrete A/C/G/T queries.

    This is NOT probe optimization. It only allows BLAST to evaluate an
    existing IUPAC probe because short-query BLAST works best with concrete
    nucleotide queries.
    """
    sequence = sequence.upper()

    invalid = sorted(set(sequence) - set(PROBE_IUPAC_TO_BASES))
    if invalid:
        raise ValueError(
            "Unsupported probe character(s): " + ", ".join(invalid)
        )

    degeneracy = 1
    for symbol in sequence:
        degeneracy *= len(PROBE_IUPAC_TO_BASES[symbol])

    if degeneracy > max_variants:
        raise RuntimeError(
            f"Probe degeneracy is {degeneracy}, which would require more "
            f"than {max_variants} concrete BLAST queries. "
            "This BLAST-only version stops instead of silently truncating "
            "the probe variants."
        )

    variants = [""]

    for symbol in sequence:
        variants = [
            prefix + base
            for prefix in variants
            for base in PROBE_IUPAC_TO_BASES[symbol]
        ]

    return variants
