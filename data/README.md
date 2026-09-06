# Input data

This directory contains user-provided nucleotide FASTA datasets.

## design/

FASTA datasets used by `01_varvamp_assay_design.py` for sequence preprocessing,
multiple-sequence alignment, conservation analysis, and VarVAMP assay design.

## validation/

FASTA datasets used by `02_in_silico_validation.py` for independent or broader
in silico validation with MFEprimer and BLAST+.

The validation dataset should ideally be larger and/or independent from the
design dataset when possible.

Large FASTA files are not intended to be committed to Git..