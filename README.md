# PCR Assay Design and In Silico Validation Workflow

[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](environment.yml)
[![Docker](https://img.shields.io/badge/Docker-supported-blue.svg)](Dockerfile)

A reproducible command-line workflow for **PCR/qPCR assay design from variable nucleotide sequence datasets** and **in silico evaluation of primer-pair and probe coverage**.

The repository connects sequence preprocessing, multiple-sequence alignment, conservation analysis, assay design with **VarVAMP**, primer-pair quality control with **MFEprimer**, and qPCR probe evaluation with **BLAST+** in a run-specific, provenance-aware workflow.

> **Scope.** This project performs computational assay design and in silico sequence-level validation. It does not replace experimental optimization or wet-lab validation of amplification efficiency, analytical sensitivity, analytical specificity, or diagnostic performance.

---

## Overview

The repository contains two complementary workflows.

| Workflow | Purpose |
|---|---|
| `01_varvamp_assay_design.py` | Preprocess nucleotide sequences, generate a final multiple-sequence alignment, quantify conservation, and design assays with VarVAMP |
| `02_in_silico_validation.py` | Evaluate selected primer pairs with MFEprimer and, for QPCR assays, evaluate the existing VarVAMP probe against retained PCR products with BLAST+ |

Workflow 01 supports all three VarVAMP assay modes:

- **SINGLE** — individual PCR amplicons
- **QPCR** — primer pairs with an internal qPCR probe
- **TILED** — overlapping amplicons for tiled sequencing

The two workflows are connected by `assay_design_manifest.json`, allowing a completed design run to be selected explicitly during downstream validation.

---

## Scientific workflow

```text
Design FASTA
    │
    ▼
┌────────────────────────────────────┐
│ Workflow 01 — Assay design         │
└────────────────────────────────────┘
    │
    ├─ FASTA validation
    ├─ Sequence orientation
    ├─ Circular-sequence normalization (optional)
    ├─ Redundancy handling
    ├─ MAFFT multiple-sequence alignment
    ├─ Alignment QC and optional trimAl trimming
    ├─ Conservation analysis
    └─ VarVAMP assay design
           │
           ▼
   assay_design_manifest.json
           │
           │  reproducible handoff
           ▼
Validation FASTA
    │
    ▼
┌────────────────────────────────────┐
│ Workflow 02 — In silico validation │
└────────────────────────────────────┘
    │
    ├─ Validation-database preparation
    ├─ MFEprimer primer-pair QC
    ├─ Predicted-product classification
    ├─ PCR target-coverage analysis
    └─ QPCR only
         ├─ IUPAC probe expansion
         ├─ blastn-short against retained PCR products
         ├─ mismatch and alignment analysis
         └─ combined PCR + probe coverage
```

The workflow deliberately keeps major scientific choices explicit instead of hiding them behind fixed defaults.

---

## Design and validation datasets

The workflow distinguishes two sequence collections with different scientific roles:

```text
data/
├── design/
│   └── <design_dataset>.fasta
└── validation/
    └── <validation_dataset>.fasta
```

### Design dataset

`data/design/` contains the sequences used to construct the assay. Workflow 01 uses this dataset for preprocessing, alignment, conservation analysis, and VarVAMP design.

### Validation dataset

`data/validation/` contains the sequences used to evaluate the resulting assay. Workflow 02 uses this dataset for predicted PCR-product analysis, target coverage, and optional qPCR probe evaluation.

For a generalizability assessment, the validation dataset should ideally be **independent of**, and often larger or more diverse than, the design dataset.

The reported coverage metrics always refer to the supplied validation dataset after any user-selected length filtering. They should not be interpreted as organism-level analytical specificity unless the validation database was deliberately constructed to answer that question.

---

## Key features

- Nucleotide FASTA validation with IUPAC ambiguity support
- Optional orientation normalization with MAFFT
- Optional cyclic start-position normalization of complete circular sequences with MARS
- Redundancy handling with:
  - no filtering
  - SeqKit exact deduplication
  - CD-HIT-EST clustering
- Multiple MAFFT alignment strategies, from scalable progressive methods to high-accuracy iterative methods
- Alignment QC and optional trimAl trimming
- Per-position conservation analysis
- VarVAMP `SINGLE`, `QPCR`, and `TILED` assay design
- Run-specific result directories that avoid silent overwriting
- Machine-readable design and validation manifests
- MFEprimer 3.x primer-pair QC and predicted-product extraction
- Product classification into `TARGET`, `SELF_PRIMING`, `CROSS_SCHEME`, and `OTHER`
- PCR target coverage calculated from unique biological sequence identifiers (`HitID`)
- QPCR probe evaluation only after primer-pair filtering
- Expansion of existing IUPAC-degenerate probes into concrete A/C/G/T variants
- `blastn-short` probe analysis with exact, mismatch, partial, and no-hit categories
- Human-readable probe alignment reports and mismatch-position tables
- Docker-based reproducible execution
- Portable project-relative provenance paths for files stored inside the repository

---

## Repository structure

```text
PCR-assay-design-workflow/
├── 01_varvamp_assay_design.py
├── 02_in_silico_validation.py
├── Dockerfile
├── compose.yaml
├── .dockerignore
├── environment.yml
├── README.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── CHANGELOG.md
├── .gitignore
│
├── data/
│   ├── README.md
│   ├── design/
│   └── validation/
│
├── work/
└── results/
```

The storage areas have distinct purposes:

| Directory | Role |
|---|---|
| `data/` | User-provided biological input datasets |
| `work/` | Reconstructible intermediate and technical files |
| `results/` | Scientific outputs, summaries, logs, and provenance records |

Each execution receives an isolated `run_id`, normally a timestamp such as `20260909_011500`.

```text
work/<project>/design/<run_id>/
results/<project>/design/<run_id>/

work/<project>/validation/<run_id>/
results/<project>/validation/<run_id>/
```

---

## Installation

### Docker — recommended

**Docker installation is self-contained.** Install Docker with the Compose plugin
(Docker Desktop includes it). No local Conda, Python, MAFFT, CD-HIT-EST, SeqKit,
trimAl, BLAST+, VarVAMP, MARS, MFEprimer, or Python libraries are required.
Conda itself is included in the image, and the environment
`pcr-assay-design-workflow` is created during the build.

From the repository root:

```bash
docker compose build
docker compose run --rm design
docker compose run --rm validation
```

Both interactive services use the same image and explicitly launch the script
through `conda run --no-capture-output -n pcr-assay-design-workflow`. No host-side
`conda activate` is needed. Run them from an interactive terminal.

The repository is bind-mounted at `/app`, including local `data/`, `work/`, and
`results/`: inputs and outputs remain accessible on the host after the container
exits. Script edits are immediately visible; rebuild after changing `Dockerfile`
or `environment.yml`. Files created by the default container user may be owned
by root on Linux. Store inputs within the repository or add an explicit volume
for external datasets; host absolute paths are not automatically available.

The image includes the following runtime dependencies:

| Dependency | Version / installation |
|---|---|
| Linux and Conda | Miniconda image `26.7.1`, pinned to its amd64 manifest digest in `Dockerfile` |
| Python | Conda, 3.11.11 |
| MAFFT; CD-HIT (`cd-hit-est`) | Conda, 7.505; 4.8.1 |
| SeqKit; trimAl | Conda, 2.9.0; 1.5.0 |
| BLAST+ (`blastn`, `makeblastdb`) | Conda, 2.16.0 |
| Biopython; pandas; NumPy | Conda, 1.85; 2.2.3; 1.26.4 |
| matplotlib (`matplotlib-base`); Pillow; PyMuPDF (`fitz`) | Conda, 3.9.4; 11.1.0; 1.25.3 |
| pip; setuptools (`pkg_resources` for seqfold) | Conda, 25.0.1; 75.8.0 |
| VarVAMP; primer3-py; seqfold | pip, 1.3.2; 2.0.3; 0.7.18 |
| MARS | Compiled from commit `cbf8f594e96db8f757f1d84edf773406b4c72701`, including bundled SeqAn, SDSL/divsufsort and Edlib sources |
| MFEprimer | Official Linux amd64 binary, 3.1.0 |
| C++/OpenMP runtime and CA certificates | apt: `libstdc++6`, `libgomp1`, `ca-certificates` |

MARS and MFEprimer downloads are SHA-256 checked. MARS compilation uses
`build-essential`, `cmake`, `git`, `curl`, `unzip`, and `ca-certificates` in a
separate build stage; these build dependencies are not copied into the runtime.
The build checks executable availability, versions, Python imports, pip
consistency, and both scripts' `--help` without running a biological analysis.
Workflow 02 directly imports only the Python standard library; `makeblastdb`
is included with BLAST+ although the current script uses `blastn -subject`.

Direct application versions are pinned. Apt packages, Conda build variants and
transitive dependencies are **not fully locked**, so a fresh solve is not
promised to be bit-for-bit identical. Exact installed package inventories are
saved inside the image in `/opt/workflow-manifests/` (`conda-explicit.txt`,
`pip-freeze.txt`, `dpkg-packages.txt`). Retain the built image by digest for
identical reruns. Building requires Internet access; running uses dependencies
already inside the image.

Rebuild from scratch and perform non-destructive startup checks:

```bash
docker compose config --quiet
docker compose build --no-cache
docker compose run --rm design --help
docker compose run --rm validation --help
```

For a terminal/input check without starting an analysis:

```bash
docker compose run --rm --entrypoint conda design run --no-capture-output -n pcr-assay-design-workflow python3 -c 'import sys; assert sys.stdin.isatty() and sys.stdout.isatty(); print(input("Terminal check: type OK then Enter: "))'
```

**Architecture:** Compose explicitly selects `linux/amd64` because MFEprimer
3.1.0 is provided as an x86_64 binary. ARM hosts require Docker's amd64 emulation;
this is not a native ARM image, and emulation may be slower. For a direct build,
use `docker build --platform linux/amd64 -t pcr-assay-design-workflow:latest .`.

### Conda — native execution

Create the environment:

```bash
conda env create -f environment.yml
conda activate pcr-assay-design-workflow
```

The Conda environment installs the main dependencies, including Python 3.11, MAFFT, CD-HIT-EST, SeqKit, trimAl, BLAST+, Biopython, pandas, matplotlib, Pillow, PyMuPDF, and VarVAMP.

For native execution, **MARS** and **MFEprimer** must additionally be installed and available in `PATH`.

---

## Quick start

### 1. Add the design dataset

```text
data/design/my_design_dataset.fasta
```

Run:

```bash
python3 01_varvamp_assay_design.py
```

or with Docker:

```bash
docker compose run --rm design
```

The workflow interactively asks for the major preprocessing, alignment, trimming, and VarVAMP design choices.

### 2. Add the validation dataset

```text
data/validation/my_validation_dataset.fasta
```

Run:

```bash
python3 02_in_silico_validation.py
```

or with Docker:

```bash
docker compose run --rm validation
```

Workflow 02 lists completed Workflow 01 design runs, reads the selected assay-design manifest, and then evaluates the chosen assay against the validation dataset.

Use `--verbose` with either script to display complete external commands and detailed paths. Commands are also retained in the run-specific log files.

---

# Workflow 01 — Assay design

## 1. Input validation

Workflow 01 accepts nucleotide FASTA data containing standard nucleotides and IUPAC ambiguity codes.

FASTA record identifiers should be unique so that sequences can be tracked consistently through preprocessing and downstream analyses.

## 2. Sequence orientation

Available strategies are:

- keep the original orientation
- MAFFT `--adjustdirection`
- MAFFT `--adjustdirectionaccurately`

Orientation normalization is useful when sequences representing the same biological region are not consistently provided in the same strand orientation.

## 3. Sequence topology

Datasets can be treated as:

- **linear**
- **circular**

For complete circular sequences, MARS can cyclically shift sequence starts before the final multiple-sequence alignment. This changes the linearization start position, not the biological nucleotide content.

## 4. Redundancy handling

Three strategies are available:

- no redundancy reduction
- exact duplicate removal with SeqKit
- similarity-based clustering with CD-HIT-EST

For CD-HIT-EST, the workflow exposes identity threshold, clustering mode, and strand comparison as explicit scientific parameters.

## 5. Final multiple-sequence alignment

Workflow 01 provides a range of MAFFT strategies, including:

- Auto
- FFT-NS-1
- FFT-NS-2
- FFT-NS-i
- NW-NS-2
- NW-NS-i
- L-INS-i
- G-INS-i
- E-INS-i
- NW-NS-PartTree-1

The selected strategy is recorded for provenance.

### Computational resources

Workflow 01 supports multithreaded execution for tools that expose a thread option.

In the current version, if `--threads` is omitted, the script detects the logical CPU threads available to the process and keeps two threads free when possible. A fixed value can still be requested explicitly:

```bash
python3 01_varvamp_assay_design.py --threads 8
```

The selected thread count is displayed at startup and recorded with the run.

## 6. Alignment QC and optional trimming

The workflow calculates alignment-level quality metrics, including:

- alignment length
- gap content
- highly gapped columns
- fully occupied columns
- per-sequence occupancy
- low-occupancy sequence counts

Optional trimAl strategies include:

- `-noallgaps`
- `-gappyout`
- `-automated1`
- manual gap-threshold filtering
- no trimming

## 7. Conservation analysis

For each alignment position, Workflow 01 calculates:

- A/C/G/T counts
- gap count
- ambiguous-character count
- occupancy
- major base
- major-base frequency
- Shannon entropy
- strict conservation
- threshold-based conservation

Default thresholds are:

```text
minimum occupancy       = 0.95
minimum major frequency = 0.95
```

They can be modified with:

```bash
--min-occupancy
--min-major-frequency
```

## 8. VarVAMP assay design

Workflow 01 supports all VarVAMP modes:

| Mode | Intended design |
|---|---|
| `SINGLE` | Individual PCR amplicons |
| `QPCR` | Primer pair plus internal qPCR probe |
| `TILED` | Overlapping tiled amplicons |

The workflow exposes the VarVAMP consensus threshold and ambiguity limits, as well as mode-specific options.

If a VarVAMP attempt fails during an interactive run, design parameters can be adjusted and **only the VarVAMP stage is retried**; the preceding preprocessing and alignment stages do not need to be recomputed.

---

# Workflow 02 — In silico validation

Workflow 02 consumes:

1. a completed Workflow 01 `assay_design_manifest.json`
2. a validation FASTA dataset

Its purpose is to measure sequence-level assay coverage in the supplied validation collection while retaining the exact design-to-validation provenance relationship.

## 1. Validation-database preparation

The workflow summarizes the selected validation FASTA and optionally allows sequence-length filtering.

The retained dataset is then reformatted for MFEprimer using SeqKit without changing sequence content or FASTA identifiers.

MFEprimer indexes and other reconstructible database files remain in the run-specific `work/` directory.

## 2. Primer-pair QC with MFEprimer

For each selected assay, Workflow 02 runs the LEFT and RIGHT primers through MFEprimer.

The current parser expects the **legacy MFEprimer 3.x text report**, including the `Hairpin List` and `Dimer List` sections.

The qPCR probe is **not** included in this MFEprimer primer hairpin/dimer analysis.

## 3. Predicted-product classification

MFEprimer-predicted products are classified as:

| Class | Interpretation |
|---|---|
| `TARGET` | LEFT × RIGHT or RIGHT × LEFT from the selected scheme |
| `SELF_PRIMING` | LEFT × LEFT or RIGHT × RIGHT |
| `CROSS_SCHEME` | Recognized primers belonging to different schemes |
| `OTHER` | Product that cannot be safely assigned to the selected assay |

The user can choose whether downstream coverage uses:

```text
TARGET only
```

or:

```text
TARGET + SELF_PRIMING
```

`CROSS_SCHEME` and `OTHER` products are never treated as valid products for the selected assay.

## 4. PCR target coverage

PCR coverage is calculated from **unique retained HitIDs**, not raw amplicon counts:

```text
PCR coverage
=
unique retained PCR-positive HitIDs
/
number of sequences in the final retained validation database
```

This prevents a biological sequence from being counted multiple times when MFEprimer predicts more than one retained product for the same HitID.

## 5. QPCR probe validation

For QPCR designs, the existing VarVAMP probe can be evaluated after primer-pair filtering.

Probe analysis is therefore conditional on PCR positivity:

```text
validation database
      │
      ▼
MFEprimer primer-pair analysis
      │
      ▼
retained PCR products
      │
      ▼
probe BLAST
```

### IUPAC-degenerate probes

If the original VarVAMP probe contains ambiguous IUPAC positions, Workflow 02 expands that **existing probe** into all compatible concrete A/C/G/T sequences.

The workflow does not redesign the probe and does not introduce new ambiguity.

### Probe BLAST categories

Concrete probe variants are searched against retained PCR products with `blastn-short`.

Results are classified as:

```text
0_MISMATCH
1_MISMATCH
2_MISMATCHES
GE3_MISMATCHES
PARTIAL_ONLY
NO_HIT
```

A **full-length ungapped site** requires the complete probe query to align without gaps. It may still contain mismatches and should therefore not be interpreted as proof of efficient probe hybridization.

### Probe outputs

For each QPCR assay, Workflow 02 reports both amplicon-level and unique-HitID-level results and writes, among other files:

```text
probe_blast_raw.tsv
probe_blast_amplicon_summary.tsv
probe_blast_hitid_summary.tsv
probe_blast_alignments_by_amplicon.txt
probe_blast_alignments_by_hitid.txt
probe_blast_mismatch_positions.tsv
probe_variant_statistics.tsv
probe_variant_hitid_matrix.tsv
probe_blast_summary.txt
```

The human-readable alignment reports make it possible to inspect the probe/target alignments and mismatch positions directly.

## 6. Combined PCR + probe coverage

For a final validation database containing `N` retained sequences:

```text
PCR
= unique PCR-positive HitIDs / N

PCR + probe exact
= PCR-positive HitIDs with a best concrete probe match of 0 mismatches / N

PCR + probe <=1 mismatch
= PCR-positive HitIDs with a best concrete probe match of <=1 mismatch / N

PCR + probe <=2 mismatches
= PCR-positive HitIDs with a best concrete probe match of <=2 mismatches / N

PCR + full probe site
= PCR-positive HitIDs with a full-length ungapped probe alignment / N
```

When several concrete IUPAC-derived probe variants match the same biological sequence, the target is counted only once using its best probe result.

---

## Main outputs

### Workflow 01

The principal scientific results are stored under:

```text
results/<project>/design/<run_id>/
```

They include:

- final multiple-sequence alignment
- alignment QC outputs
- conservation tables and visualizations
- selected VarVAMP output tree
- `assay_design_manifest.json`
- `workflow_summary.txt`
- `workflow.log`

### Workflow 02

The principal validation results are stored under:

```text
results/<project>/validation/<run_id>/
```

The main cross-assay summary is:

```text
summary/coverage_results.tsv
```

Additional result groups contain:

- MFEprimer reports and retained amplicons
- positive HitID lists
- non-target product tables
- primer secondary-structure summaries
- qPCR probe BLAST results and alignments
- mismatch-position tables
- per-probe-variant statistics
- `validation_manifest.json`
- `validation.log`

<details>
<summary><strong>Example validation result structure</strong></summary>

```text
results/<project>/validation/<run_id>/
├── inputs/
│   ├── assay_design_manifest.json
│   ├── assay_manifest.tsv
│   └── validation_database.tsv
│
├── mfeprimer/
│   └── <scheme>/
│       ├── <raw MFEprimer report>
│       ├── amplicons_all.tsv
│       ├── amplicons_valid.tsv
│       ├── amplicons_valid.fasta
│       ├── positive_hitids.txt
│       ├── non_target_products.tsv
│       └── coverage_summary.txt
│
├── probe_blast/
│   └── <scheme>/
│       ├── probe_blast_queries.fasta
│       ├── probe_blast_raw.tsv
│       ├── probe_blast_amplicon_summary.tsv
│       ├── probe_blast_hitid_summary.tsv
│       ├── probe_blast_alignments_by_amplicon.txt
│       ├── probe_blast_alignments_by_hitid.txt
│       ├── probe_blast_mismatch_positions.tsv
│       ├── probe_blast_summary.txt
│       ├── probe_variant_statistics.tsv
│       └── probe_variant_hitid_matrix.tsv
│
├── summary/
│   ├── mfeprimer_runs.tsv
│   ├── coverage_results.tsv
│   ├── secondary_structure_summary.tsv
│   └── secondary_structure_details.tsv
│
├── validation_manifest.json
└── validation.log
```

</details>

---

## Reproducibility and provenance

Each run is isolated by project, workflow stage, and `run_id`.

The workflow records provenance through:

```text
assay_design_manifest.json
validation_manifest.json
workflow_summary.txt
workflow.log
validation.log
```

The manifests are designed to preserve information such as:

- source FASTA paths
- sequence counts
- dataset SHA-256 checksums
- preprocessing decisions
- alignment strategy
- conservation thresholds
- VarVAMP mode and parameters
- selected assays
- validation database and filtering
- MFEprimer settings
- PCR coverage denominator
- probe-validation status
- software versions
- parent design run and validation run identifiers

Paths for files located inside the project are stored relative to the repository root whenever possible, improving portability between native execution, Docker, and relocated project directories.

For publication-quality analyses, retain the original accession list and document the source database, retrieval date, inclusion/exclusion criteria, lineage/genotype composition, sequence completeness criteria, and all user-selected scientific parameters.

---

## Scientific interpretation and limitations

### In silico coverage is not wet-lab performance

Sequence-level coverage does not directly establish:

- amplification efficiency
- limit of detection
- analytical sensitivity
- analytical specificity
- robustness to reaction conditions
- probe fluorescence performance
- diagnostic sensitivity or specificity

Experimental validation remains required.

### Validation-database composition matters

Coverage estimates are conditional on the sequences included in the validation dataset. Database bias, incomplete lineage representation, duplicated sequences, partial genomes, and temporal sampling can all affect the resulting percentages.

### Circular genomes

For circular molecules, arbitrary FASTA start positions can distort a conventional linear alignment. MARS can reduce this problem by normalizing cyclic start positions before MAFFT.

### Redundancy reduction

Clustering or exact deduplication can be useful during assay design, but validation should normally be performed against the original or another biologically representative dataset rather than only against cluster representatives.

### Probe mismatch counts

A mismatch count is descriptive. Its experimental effect depends on position, neighboring sequence context, probe chemistry, melting temperature, and reaction conditions.

### MFEprimer compatibility

Workflow 02 currently depends on the structure of the **MFEprimer 3.x text report**. Recognition of newer index formats does not imply full MFEprimer 4.x report compatibility.

---

## Software used

The workflow integrates established bioinformatics tools:

| Tool | Role |
|---|---|
| [VarVAMP](https://github.com/jonas-fuchs/varVAMP) | Degenerate primer/probe design for variable sequence alignments |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Multiple-sequence alignment and orientation normalization |
| [MARS](https://github.com/lorrainea/MARS) | Circular-sequence start-position normalization |
| [SeqKit](https://github.com/shenwei356/seqkit) | FASTA manipulation and exact deduplication |
| [CD-HIT](https://github.com/weizhongli/cdhit) | Similarity-based nucleotide clustering with CD-HIT-EST |
| [trimAl](https://github.com/inab/trimal) | Multiple-sequence-alignment trimming |
| [MFEprimer](https://github.com/quwubin/MFEprimer-3.0) | Primer quality control and predicted PCR-product analysis |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) | qPCR probe sequence matching |

Users should cite the original publications for the tools used in a reported analysis.

---

## Selected software references

- Fuchs J, Kleine J, Schemmerer M, et al. **varVAMP: degenerate primer design for tiled full genome sequencing and qPCR.** *Nature Communications*. 2025;16:5067.
- Katoh K, Standley DM. **MAFFT multiple sequence alignment software version 7: improvements in performance and usability.** *Molecular Biology and Evolution*. 2013;30:772–780. doi:10.1093/molbev/mst010.
- Ayad LAK, Pissis SP. **MARS: improving multiple circular sequence alignment using refined sequences.** *BMC Genomics*. 2017;18:86.
- Li W, Godzik A. **Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences.** *Bioinformatics*. 2006;22:1658–1659. doi:10.1093/bioinformatics/btl158.
- Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. **trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses.** *Bioinformatics*. 2009;25:1972–1973. doi:10.1093/bioinformatics/btp348.
- Wang K, Li H, Xu Y, et al. **MFEprimer-3.0: quality control for PCR primers.** *Nucleic Acids Research*. 2019;47:W610–W613. doi:10.1093/nar/gkz351.
- Shen W, Sipos B, Zhao L. **SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing.** *iMeta*. 2024;e191. doi:10.1002/imt2.191.

---

## Citation

Repository citation metadata are provided in:

```text
CITATION.cff
```

When publishing analyses generated with this workflow, please also cite the relevant third-party software used in the analysis.

---

## Contributing

Contributions, bug reports, reproducibility improvements, and scientifically motivated feature requests are welcome.

See:

```text
CONTRIBUTING.md
```

---

## License

The workflow source code is distributed under the **GNU General Public License v3.0 or later**.

See:

```text
LICENSE
```

Input datasets and generated scientific results may be subject to separate terms imposed by source databases, institutions, collaborators, or data-use agreements.