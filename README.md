# VarVAMP Assay Design and In Silico Validation Workflow

A reproducible and interactive command-line workflow for nucleotide sequence preprocessing, multiple sequence alignment, conservation analysis, VarVAMP assay design, and in silico validation.

The repository is organized into two complementary workflows:

1. `01_varvamp_assay_design.py` — sequence preprocessing, alignment, conservation analysis, and assay design with VarVAMP.
2. `02_in_silico_validation.py` — target-coverage and probe validation of VarVAMP designs with MFEprimer and BLAST+.

The workflow supports the three VarVAMP assay modes:

- **SINGLE** — conventional PCR / individual amplicons;
- **QPCR** — primer pairs with an internal qPCR probe;
- **TILED** — overlapping amplicons for tiled sequencing.

---

## Scope

This repository is designed for nucleotide FASTA datasets, including viral genomes, genes, genomic regions, bacterial genomes, and other comparable nucleotide sequence collections.

It is not intended for protein sequences.

The project performs computational assay design and in silico validation. Computational predictions alone do not establish diagnostic sensitivity, specificity, amplification efficiency, or clinical validity. Selected assays should ultimately be validated experimentally.

---

## Workflow overview

```text
Input nucleotide FASTA
        │
        ▼
┌──────────────────────────────────────────────┐
│ 01_varvamp_assay_design.py                  │
└──────────────────────────────────────────────┘
        │
        ├── FASTA validation
        ├── Sequence orientation
        │      ├── keep
        │      ├── MAFFT --adjustdirection
        │      └── MAFFT --adjustdirectionaccurately
        │
        ├── Sequence topology
        │      ├── linear
        │      └── circular → MARS rotation
        │
        ├── Redundancy handling
        │      ├── none
        │      ├── SeqKit exact deduplication
        │      └── CD-HIT-EST clustering
        │
        ├── Final MAFFT multiple sequence alignment
        ├── Alignment quality control
        ├── Optional trimAl trimming
        ├── Conservation analysis
        └── VarVAMP assay design
               ├── SINGLE
               ├── QPCR
               └── TILED
                        │
                        ▼
             assay_design_manifest.json
                        │
                        ▼
┌──────────────────────────────────────────────┐
│ 02_in_silico_validation.py                  │
└──────────────────────────────────────────────┘
        │
        ├── Detect VarVAMP assay mode
        ├── Prepare target validation database
        ├── MFEprimer primer-pair validation
        ├── Target amplicon extraction
        └── QPCR only:
               ├── expand IUPAC-degenerate probe
               ├── blastn-short
               ├── per-probe-variant statistics
               ├── HitID × probe-variant matrix
               └── global PCR + probe coverage
```

Scientific choices are collected step by step during interactive execution. The workflow does not silently impose fixed choices for orientation, topology, redundancy reduction, MAFFT, trimAl, or VarVAMP parameters.

---

## Repository structure

```text
.
├── 01_varvamp_assay_design.py
├── 02_in_silico_validation.py
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── .gitignore
├── CONTRIBUTING.md
├── CHANGELOG.md
│
├── data/
│   ├── README.md
│   └── accession_numbers.tsv
│
└── results/
    └── README.md
```

The `work/` directory is created automatically during assay design and is not intended to be tracked by Git.

Validation outputs are generated under `validation_inputs/`.

---

## Requirements

### Conda / pip dependencies

The provided `environment.yml` installs:

- Python 3.11;
- MAFFT;
- CD-HIT-EST;
- SeqKit;
- trimAl;
- BLAST+;
- Biopython;
- pandas;
- matplotlib;
- Pillow;
- PyMuPDF;
- VarVAMP 1.3.2.

### Additional external programs

Two programs are currently installed separately from the Conda environment:

**MARS** — used by Workflow 01 for cyclic start-position normalization of complete circular nucleotide sequences.

**MFEprimer** — used by Workflow 02 for primer-pair specificity and target-coverage analysis.

The expected executable names are:

```text
mars
mfeprimer
```

They must be available in the system `PATH`, unless custom executable paths are supplied where supported.

---

## Installation with Conda

From the repository root:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate varvamp-qpcr-workflow
```

Verify the main Conda-installed programs:

```bash
python --version
mafft --version
cd-hit-est -h
seqkit version
trimal --version
blastn -version
varvamp --help
```

Verify the separately installed programs:

```bash
mars --help
mfeprimer --help
```

---

## Input data

The primary input is a nucleotide FASTA file containing at least one sequence.

Accepted characters include standard DNA/RNA nucleotides and IUPAC ambiguity codes. FASTA record IDs must be unique so that sequences can be tracked through orientation and rotation steps.

The FASTA file may be stored under:

```text
data/my_sequences.fasta
```

or supplied from any valid filesystem location.

For reproducible studies, users should document accession numbers, source database, retrieval date, search strategy, inclusion/exclusion criteria, genotype or lineage information, sequence completeness, and analysed genomic region.

---

# Workflow 01 — Assay design

## Interactive execution

```bash
python3 01_varvamp_assay_design.py
```

For detailed external commands and file paths:

```bash
python3 01_varvamp_assay_design.py --verbose
```

In normal mode, external commands are recorded in:

```text
results/<project>/workflow.log
```

### Step 1 — Sequence orientation

The workflow allows three strategies:

```text
1. Keep the original orientation
2. MAFFT --adjustdirection
3. MAFFT --adjustdirectionaccurately
```

### Step 2 — Sequence topology

The workflow asks whether sequences are `linear` or `circular`. For complete circular datasets, MARS can normalize cyclic start positions before the final alignment.

### Step 3 — Redundancy handling

Supported strategies:

```text
none
SeqKit exact deduplication
CD-HIT-EST clustering
```

For CD-HIT-EST, the user can control identity threshold, clustering mode, and strand comparison. A compatible word size is selected automatically.

### Step 4 — Final MAFFT alignment

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

### Step 5 — Alignment QC and trimAl

The workflow calculates alignment length, gap content, column gap statistics, sequence occupancy, and low-occupancy sequence counts. A transparent heuristic assessment is shown before trimming.

Available trimAl strategies are:

```text
-noallgaps
-gappyout
-automated1
-gt VALUE
no trimming
```

### Conservation analysis

For each alignment position, the workflow calculates A/C/G/T counts, gaps, ambiguous characters, occupancy, major base, major-base frequency, Shannon entropy, strict conservation, and threshold-based conservation.

Default thresholds:

```text
minimum occupancy       = 0.95
minimum major frequency = 0.95
```

They can be changed with:

```bash
--min-occupancy
--min-major-frequency
```

### VarVAMP assay design

Workflow 01 supports all three VarVAMP modes.

**SINGLE** — conventional PCR / individual amplicon design.

**QPCR** — primer-pair plus internal-probe design.

Main parameters include:

```text
-t     VarVAMP consensus threshold
-a     maximum ambiguous positions per primer
-pa    maximum ambiguous positions in the probe
```

**TILED** — overlapping tiled amplicon design.

---

## Assay-design manifest

At the end of Workflow 01, the script generates:

```text
results/<project>/assay_design_manifest.json
```

and copies it into the selected VarVAMP result directory.

The manifest records the project name, input FASTA, sequence counts, orientation, topology, redundancy strategy, MAFFT strategy, trimming strategy, conservation thresholds, VarVAMP mode and parameters, result directory, software versions, and workflow log path.

This file provides a reproducible handoff between assay design and in silico validation.

---

# Workflow 02 — In silico validation

Run:

```bash
python3 02_in_silico_validation.py
```

Detailed mode:

```bash
python3 02_in_silico_validation.py --verbose
```

Workflow 02 detects the available VarVAMP design and adapts validation to the selected assay mode.

## Target validation database

The script asks the user to choose a local nucleotide FASTA database. SeqKit is used to produce a one-sequence-line FASTA copy compatible with the MFEprimer validation workflow:

```bash
seqkit seq -w 0
```

The script verifies that sequence number, sequence length, and ambiguous-base content are preserved.

## MFEprimer validation

MFEprimer evaluates primer pairs against the selected target database. The workflow reports potential products, valid LEFT/RIGHT products, self-priming products, unique PCR-positive HitIDs, and primer-pair target coverage.

Biological target coverage is calculated using unique HitIDs rather than simply counting predicted amplicons.

## QPCR probe validation

For QPCR designs, Workflow 02 additionally validates the internal probe with local BLAST+ using `blastn-short` against valid MFEprimer amplicons.

### IUPAC-degenerate probes

If a VarVAMP probe contains ambiguous IUPAC positions, Workflow 02 expands it into every concrete A/C/G/T-compatible version before BLAST.

For example:

```text
Y = C/T
K = G/T
```

A probe containing one `Y` and one `K` produces four concrete variants.

No new ambiguous bases are introduced; only ambiguity already present in the original VarVAMP probe is expanded.

### Per-probe-variant analysis

Every concrete probe version is analysed independently at unique-PCR-positive-HitID level. The workflow reports:

```text
Exact match
1 mismatch
2 mismatches
>=3 mismatches
Partial/gapped alignment
No hit
Full-length probe site
```

These rows must not be naively added together because one target can match more than one probe variant.

The workflow also calculates the best result across all concrete variants for each biological target.

### HitID × probe-variant matrix

For degenerate probes, Workflow 02 writes:

```text
probe_variant_hitid_matrix.tsv
probe_variant_statistics.tsv
```

The matrix records the category obtained by every concrete probe variant against every PCR-positive HitID.

### Probe mismatch categories

```text
0_MISMATCH
1_MISMATCH
2_MISMATCHES
GE3_MISMATCHES
PARTIAL_ONLY
NO_HIT
```

A full-length probe alignment requires complete query coverage without a gap.

Mismatch categories are descriptive computational results. Experimental effects depend on mismatch position, nucleotide substitution, probe chemistry, melting temperature, and reaction conditions.

---

## Final qPCR coverage

The final comparison uses the **complete validation database** as denominator.

For a validation database containing `N` sequences:

```text
PCR
= unique PCR-positive targets / N

PCR + probe exact
= PCR-positive targets with an exact probe match / N

PCR + probe <=1 MM
= PCR-positive targets whose best probe variant has 0 or 1 mismatch / N

PCR + probe <=2 MM
= PCR-positive targets whose best probe variant has 0, 1, or 2 mismatches / N

PCR + full probe site
= PCR-positive targets with a full-length ungapped probe alignment / N
```

Thus, `PCR` represents primer-pair coverage only, whereas each `PCR + probe` column represents the complete primer-plus-probe system.

---

## Workflow 02 outputs

Typical validation outputs include:

```text
validation_inputs/
├── manifest.tsv
├── validation.log
├── primer_pairs/
├── probes/
├── schemes/
├── databases/
│
└── mfeprimer_coverage/
    └── <database>/
        ├── mfeprimer_runs.tsv
        ├── coverage_results.tsv
        │
        └── <scheme>/
            ├── amplicons_valid.tsv
            ├── amplicons_valid.fasta
            ├── positive_hitids.txt
            ├── non_target_products.tsv
            ├── coverage_summary.txt
            ├── probe_blast_raw.tsv
            ├── probe_blast_amplicon_summary.tsv
            ├── probe_blast_hitid_summary.tsv
            ├── probe_blast_alignments_by_amplicon.txt
            ├── probe_blast_alignments_by_hitid.txt
            ├── probe_blast_mismatch_positions.tsv
            ├── probe_blast_summary.txt
            ├── probe_variant_statistics.tsv
            └── probe_variant_hitid_matrix.tsv
```

Probe-specific files are generated only for QPCR designs.

---

## Non-interactive execution of Workflow 01

Example for a QPCR design:

```bash
python3 01_varvamp_assay_design.py \
  --input data/my_sequences.fasta \
  --orientation accurate \
  --topology circular \
  --redundancy cdhit \
  --identity 0.95 \
  --cdhit-mode accurate \
  --cdhit-strand both \
  --mafft-strategy auto \
  --trimming none \
  --varvamp-mode qpcr \
  --varvamp-threshold 0.95 \
  --primer-ambiguity 2 \
  --probe-ambiguity 2 \
  --threads 8
```

Display available options with:

```bash
python3 01_varvamp_assay_design.py --help
python3 02_in_silico_validation.py --help
```

---

## Output organization of Workflow 01

For an input named:

```text
my_sequences.fasta
```

the default project name is `my_sequences`.

Intermediate files are written to:

```text
work/my_sequences/
```

Project results are written to:

```text
results/my_sequences/
```

Typical outputs include:

```text
my_sequences_alignment.fasta
my_sequences_alignment_trimmed.fasta
my_sequences_pretrim_alignment_qc.txt
my_sequences_posttrim_alignment_qc.txt
my_sequences_conservation_by_position.csv
my_sequences_consensus.fasta
my_sequences_conservation_summary.txt
workflow.log
workflow_summary.txt
assay_design_manifest.json

varvamp_single/
varvamp_qpcr/
or
varvamp_tiled/
```

---

## Reproducibility

For reproducible analysis, retain or report:

1. sequence database and retrieval date;
2. accession numbers;
3. inclusion and exclusion criteria;
4. sequence completeness;
5. orientation strategy;
6. sequence topology;
7. MARS use, when applicable;
8. redundancy strategy;
9. CD-HIT-EST parameters, when applicable;
10. MAFFT strategy;
11. alignment QC;
12. trimAl strategy and threshold, when applicable;
13. conservation-analysis thresholds;
14. VarVAMP mode;
15. VarVAMP consensus threshold;
16. primer ambiguity limit;
17. probe ambiguity limit for QPCR;
18. mode-specific VarVAMP parameters;
19. software versions;
20. validation FASTA database;
21. MFEprimer parameters;
22. target-amplicon extraction filter;
23. BLAST probe-analysis parameters;
24. primer coverage on the complete validation dataset;
25. probe-match statistics;
26. concrete IUPAC probe variants;
27. experimental validation of the selected assay.

Workflow 01 automatically records many of these settings in:

```text
assay_design_manifest.json
workflow_summary.txt
workflow.log
```

Workflow 02 writes corresponding validation tables and logs.

---

## Important interpretation notes

**Cluster representatives:** if clustering is used for design, validation should ideally also be assessed against the original or another appropriate validation dataset, not only cluster representatives.

**Circular genomes:** arbitrary FASTA start positions can affect linearized analyses. MARS normalization helps standardize cyclic starts before alignment and assay design.

**Primer/probe mismatches:** computational mismatch counts are descriptive. Experimental impact depends on mismatch position, chemistry, melting temperature, reaction conditions, and other factors.

**Internal versus independent validation:** validation against the same dataset used for design measures internal target coverage. Independent external sequence datasets provide stronger evidence of assay generalizability.

---

## Third-party software

This repository does not redistribute third-party software.

The workflow uses or can use:

- CD-HIT/CD-HIT-EST;
- MAFFT;
- SeqKit;
- MARS;
- trimAl;
- VarVAMP;
- MFEprimer;
- NCBI BLAST+.

Users should cite the original software publications when reporting analyses produced using these tools.

---

## Citation

Citation metadata for this repository are provided in:

```text
CITATION.cff
```

---

## License

The workflow code in this repository is licensed under the **GNU General Public License v3.0 or later**.

See:

```text
LICENSE
```

Input sequence datasets and generated scientific results may be subject to separate terms depending on source databases, institutional policies, collaborators, or data-use agreements.