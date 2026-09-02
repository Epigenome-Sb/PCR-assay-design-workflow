# Results

Workflow 01 creates one result subdirectory per project:

```text
results/<project_name>/
```

This directory contains the main assay-design outputs produced by:

```text
01_varvamp_assay_design.py
```

The in silico validation workflow:

```text
02_in_silico_validation.py
```

stores its generated files separately under:

```text
validation_inputs/
```

---

## Workflow 01 outputs

Depending on the selected preprocessing, alignment, trimming, and VarVAMP options, a project result directory may contain:

- the final multiple sequence alignment;
- pre-trimming and post-trimming alignment QC reports;
- the conservative consensus sequence;
- per-position conservation metrics;
- a conservation summary;
- the workflow execution log;
- a human-readable workflow summary;
- a machine-readable assay-design manifest;
- VarVAMP outputs for `SINGLE`, `QPCR`, or `TILED` mode;
- tabular assay-design files;
- BED coordinate files;
- PDF visualizations;
- PNG visualizations.

Typical project-level files include:

```text
<project_name>_alignment.fasta
<project_name>_alignment_trimmed.fasta
<project_name>_pretrim_alignment_qc.txt
<project_name>_posttrim_alignment_qc.txt
<project_name>_consensus.fasta
<project_name>_conservation_by_position.csv
<project_name>_conservation_summary.txt
workflow.log
workflow_summary.txt
assay_design_manifest.json
```

Not every file is guaranteed to be generated. Output availability depends on the selected workflow steps and whether trimming or VarVAMP design is performed.

---

## VarVAMP result directories

Depending on the selected assay mode, Workflow 01 generates one of:

```text
results/<project_name>/varvamp_single/
results/<project_name>/varvamp_qpcr/
results/<project_name>/varvamp_tiled/
```

The exact files depend on the VarVAMP mode and whether candidate assays satisfy the selected filters.

### SINGLE

Typical files may include:

```text
primer.tsv
primer_to_amplicon_assignments.tabular
primers.bed
amplicons.bed
```

### QPCR

Typical files may include:

```text
qpcr_primers.tsv
qpcr_design.tsv
primers.bed
amplicons.bed
amplicon_plot.pdf
per_base_mismatches.pdf
amplicons_overview.png
primers_overview.png
ambiguous_consensus.fasta
varvamp_log.txt
```

### TILED

Typical files may include:

```text
primer.tsv
primer_to_amplicon_assignments.tabular
primers.bed
amplicons.bed
```

Additional VarVAMP-generated files may also be present.

---

## Assay-design manifest

Workflow 01 generates:

```text
assay_design_manifest.json
```

This file records machine-readable metadata describing the assay-design run.

It may include:

- project name;
- input FASTA path;
- sequence counts;
- orientation strategy;
- sequence topology;
- redundancy strategy;
- CD-HIT-EST parameters, when applicable;
- MAFFT strategy;
- trimming strategy;
- conservation thresholds;
- VarVAMP execution status;
- VarVAMP assay mode;
- VarVAMP parameters;
- result directory;
- software-version information;
- workflow log path.

A copy is also placed in the selected VarVAMP result directory to facilitate downstream validation.

---

## Workflow log

The file:

```text
workflow.log
```

records external commands and execution information.

This is particularly useful because normal terminal output is intentionally kept compact.

When the workflow is run with:

```bash
--verbose
```

more detailed information is also displayed in the terminal.

---

## Workflow summary

The file:

```text
workflow_summary.txt
```

provides a concise human-readable description of the main processing and assay-design choices used for the project.

Together with:

```text
assay_design_manifest.json
workflow.log
```

it helps preserve the reproducibility of the analysis.

---

# Parameter-dependent results

Workflow 01 results depend directly on the scientific parameters selected during execution.

For reproducibility, the relevant parameters should be retained for every analysis.

## Sequence orientation

Record whether sequences were:

```text
kept unchanged
```

or reoriented using:

```text
MAFFT --adjustdirection
MAFFT --adjustdirectionaccurately
```

Orientation handling can influence downstream alignment and assay coordinates.

---

## Sequence topology

Record whether the dataset was treated as:

```text
linear
```

or:

```text
circular
```

For circular datasets, also record whether MARS cyclic start-position normalization was performed.

---

## Redundancy handling

Record whether redundancy reduction used:

```text
none
SeqKit exact deduplication
CD-HIT-EST clustering
```

When CD-HIT-EST is used, record:

- identity threshold (`-c`);
- clustering mode (`-g`);
- strand comparison (`-r`);
- automatically selected word size (`-n`).

Different redundancy settings may change the number and composition of sequences used for alignment and assay design.

---

## MAFFT

Record the selected multiple sequence alignment strategy.

Supported strategies include:

```text
Auto
FFT-NS-1
FFT-NS-2
FFT-NS-i (2 cycles)
FFT-NS-i (up to 1000 cycles)
NW-NS-2
NW-NS-i (2 cycles)
NW-NS-i (up to 1000 cycles)
L-INS-i
G-INS-i
E-INS-i
NW-NS-PartTree-1
```

Alignment strategy may influence gap placement, conservation estimates, consensus generation, and downstream VarVAMP designs.

---

## Alignment trimming

Record whether trimAl was applied.

Supported strategies include:

```text
-noallgaps
-gappyout
-automated1
-gt <value>
```

or no trimming.

When a manual `-gt` threshold is used, record the selected value.

---

## Conservation analysis

Record:

- minimum occupancy threshold;
- minimum major-base frequency threshold.

These parameters determine which alignment positions contribute to the conservative consensus.

---

## VarVAMP

Record the selected assay mode:

```text
single
qpcr
tiled
```

and the common design parameters:

- consensus threshold (`-t`);
- maximum primer ambiguity (`-a`);
- number of threads.

For QPCR, also record:

- maximum probe ambiguity (`-pa`).

Mode-specific parameters should also be retained when used.

Examples include:

```text
-ol
-ml
-n
-o
-d
```

depending on the selected VarVAMP mode.

---

# Workflow 02 validation results

Workflow 02 does not write its main results under `results/`.

Validation files are generated separately under:

```text
validation_inputs/
```

Typical validation outputs may include:

```text
validation.log
manifest.tsv
mfeprimer_runs.tsv
coverage_results.tsv
amplicons_valid.tsv
amplicons_valid.fasta
positive_hitids.txt
coverage_summary.txt
```

For QPCR assays, additional files may include:

```text
probe_blast_raw.tsv
probe_blast_amplicon_summary.tsv
probe_blast_hitid_summary.tsv
probe_blast_mismatch_positions.tsv
probe_blast_summary.txt
probe_variant_statistics.tsv
probe_variant_hitid_matrix.tsv
```

These files document primer-pair target coverage and qPCR probe matching.

---

# Reproducibility metadata

When preserving or publishing results, retain:

- project name;
- input dataset name or version;
- source database;
- sequence retrieval date;
- accession numbers;
- inclusion and exclusion criteria;
- sequence completeness criteria;
- orientation strategy;
- topology;
- MARS usage, when applicable;
- redundancy strategy;
- CD-HIT-EST parameters, when applicable;
- MAFFT strategy;
- trimAl strategy;
- conservation thresholds;
- VarVAMP mode;
- VarVAMP parameters;
- software versions;
- workflow version;
- execution date;
- validation dataset;
- MFEprimer parameters;
- BLAST parameters for QPCR validation;
- any manual filtering or post-processing.

Workflow 01 automatically stores many of these details in:

```text
assay_design_manifest.json
workflow_summary.txt
workflow.log
```

---

# Interpretation

Generated outputs are computational predictions.

They do not, on their own, demonstrate:

- analytical specificity;
- analytical sensitivity;
- diagnostic sensitivity;
- diagnostic specificity;
- amplification efficiency;
- experimental primer performance;
- experimental probe performance;
- absence of off-target amplification;
- experimental validity;
- clinical validity.

Primer and probe mismatch categories from Workflow 02 are descriptive computational results and should not automatically be interpreted as successful or failed experimental detection.

Candidate assays should undergo appropriate independent in silico evaluation and experimental validation before diagnostic or research use.

---

# Git tracking

Generated project result directories are normally excluded from Git through:

```gitignore
results/*
!results/README.md
```

Workflow 02 generated files are normally excluded through:

```gitignore
validation_inputs/
```

This prevents large or dataset-specific outputs from being committed accidentally while preserving this documentation file.

Results intended for publication or public reproducibility can instead be archived separately, for example in:

- a GitHub release;
- supplementary material;
- Zenodo;
- another scientific data repository.