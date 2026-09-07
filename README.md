# PCR Assay Design and In Silico Validation Workflow

A reproducible and interactive command-line workflow for nucleotide sequence preprocessing, multiple-sequence alignment, conservation analysis, VarVAMP assay design, and in silico validation.

The repository contains two complementary workflows:

1. `01_varvamp_assay_design.py` — sequence preprocessing, alignment, conservation analysis, and assay design with VarVAMP.

2. `02_in_silico_validation.py` — primer-pair target coverage and, for qPCR designs, probe validation with MFEprimer and BLAST+.

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

## Two-dataset design

The workflows distinguish two biological sequence datasets with different roles.

```text
data/
├── design/
│   └── <design_dataset>.fasta
└── validation/
    └── <validation_dataset>.fasta
```

### Design dataset

`data/design/` contains the FASTA used by Workflow 01 for:

- preprocessing;

- multiple-sequence alignment;

- conservation analysis;

- VarVAMP assay design.

### Validation dataset

`data/validation/` contains the FASTA used by Workflow 02 for:

- MFEprimer primer-pair validation;

- PCR target-coverage calculation;

- qPCR probe BLAST analysis.

The validation dataset may be larger and more diverse than the design dataset. An independent validation dataset provides stronger evidence of assay generalizability than validation against the same sequences used for design.

---

## Workflow overview

```text
data/design/<design_dataset>.fasta
        │
        ▼
┌──────────────────────────────────────────────┐
│ 01_varvamp_assay_design.py                   │
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
        ├── Final MAFFT multiple-sequence alignment
        ├── Alignment quality control
        ├── Optional trimAl trimming
        ├── Conservation analysis
        └── VarVAMP assay design
               ├── SINGLE
               ├── QPCR
               └── TILED
                        │
                        ▼
 results/<project>/design/<run_id>/
        assay_design_manifest.json
                        │
                        │ reproducible handoff
                        ▼
data/validation/<validation_dataset>.fasta
                        │
                        ▼
┌──────────────────────────────────────────────┐
│ 02_in_silico_validation.py                   │
└──────────────────────────────────────────────┘
        │
        ├── Select Workflow 01 design manifest
        ├── Read VarVAMP assay mode and outputs
        ├── Select validation FASTA
        ├── Optional sequence-length filtering
        ├── MFEprimer FASTA formatting
        ├── MFEprimer database indexing
        ├── Full MFEprimer primer QC
        ├── Potential-product classification
        ├── TARGET / SELF_PRIMING filtering decision
        ├── PCR target-coverage calculation
        └── QPCR only:
               ├── expand existing IUPAC probe
               ├── blastn-short on retained PCR products
               ├── per-probe-variant statistics
               ├── HitID × probe-variant matrix
               └── global PCR + probe coverage
```

Scientific choices remain explicit during interactive execution.

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
│   ├── design/
│   └── validation/
│
├── work/
└── results/
```

The three storage areas have distinct roles:

```text
data/      user-provided biological input datasets
work/      generated, reconstructible intermediate files
results/   permanent scientific outputs and provenance records
```

`work/` and `results/` may initially be empty. Their project/run subdirectories are created automatically by the workflows.

The old `validation_inputs/` directory is no longer used.

---

## Run-specific storage

Every workflow execution receives a `run_id`.

By default, it is a local timestamp such as:

```text
20260906_164500
```

Workflow 01 writes by default to:

```text
work/<project>/design/<run_id>/
results/<project>/design/<run_id>/
```

Workflow 02 writes by default to:

```text
work/<project>/validation/<run_id>/
results/<project>/validation/<run_id>/
```

This isolates successive runs and prevents normal timestamp-based executions from overwriting previous analyses.

If a generated run ID already exists, a suffix such as `_02`, `_03`, etc. is added.

---

## Requirements

### Conda / pip dependencies

The provided `environment.yml` installs the main workflow dependencies, including:

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

The Conda environment name used by the repository is:

```text
pcr-assay-design-workflow
```

The environment name must remain consistent between `environment.yml`, `Dockerfile`, and `compose.yaml`.

### Additional external programs

Two additional programs are required:

**MARS** — used by Workflow 01 for cyclic start-position normalization of complete circular nucleotide sequences.

**MFEprimer** — used by Workflow 02 for primer-pair QC, specificity analysis, and target-coverage analysis.

For a native Conda installation, these programs must be installed separately and available in `PATH`.

For the Docker installation provided by this repository, MARS and MFEprimer are installed automatically during image construction.

The expected executable names are:

```text
mars
mfeprimer
```

### MFEprimer compatibility note

Workflow 02 currently parses the **legacy MFEprimer 3.x text report** for amplicon extraction, product classification, and the `Hairpin List` / `Dimer List` sections.

The workflow can recognize both of these index layouts:

```text
<database>.primerqc.bin
```

or:

```text
<database>.primerqc
<database>.primerqc.fai
```

Recognition of an index layout does **not** imply full MFEprimer 4.x report compatibility. If only a `.spec.tsv` report is produced, the current parser stops instead of silently interpreting an unsupported format.

The Dockerfile currently installs MFEprimer 3.1.0 specifically to preserve compatibility with this parser.

---

## Installation with Docker

Docker is the simplest way to reproduce the complete software environment because the image contains the Conda environment plus MARS and MFEprimer.

From the repository root, build the image:

```bash
docker compose build
```

The two Compose services use the same image:

```text
pcr-assay-design-workflow:latest
```

Run Workflow 01:

```bash
docker compose run --rm design
```

Run Workflow 02:

```bash
docker compose run --rm validation
```

The Compose configuration starts commands inside the `pcr-assay-design-workflow` Conda environment automatically. No manual `conda activate` is required for these commands.

You can also run an explicit command:

```bash
docker compose run --rm design python3 01_varvamp_assay_design.py
docker compose run --rm validation python3 02_in_silico_validation.py
```

To open an interactive shell through the same environment:

```bash
docker compose run --rm validation bash
```

Then, for example:

```bash
python3 --version
python3 02_in_silico_validation.py
```

The repository is bind-mounted at runtime:

```yaml
volumes:
  - .:/app
```

Therefore edits made to the Python scripts in the local repository are visible immediately inside the container. Rebuilding is normally required after changing `Dockerfile` or `environment.yml`, but not after ordinary edits to the Python source files.

### Docker input-data policy

The default Docker configuration exposes the repository itself, including:

```text
data/design/
data/validation/
work/
results/
```

For normal Docker usage, place design FASTA files under `data/design/` and validation FASTA files under `data/validation/`.

Arbitrary host paths outside the repository, such as `/home/<user>/another_directory/`, are not mounted by default. This is intentional and keeps the default execution model reproducible and limited to the project directory.

### Current Docker architecture note

The current Dockerfile downloads the MFEprimer 3.1.0 `linux-amd64` binary. The provided build is therefore currently intended for `amd64`/`x86_64` Linux-compatible Docker execution.

---

## Installation with Conda

For native execution outside Docker, create the environment from the repository root:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate pcr-assay-design-workflow
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

MARS and MFEprimer are not installed by `environment.yml`; install them separately for native execution and verify that they are available in `PATH`:

```bash
mars --help
mfeprimer --help
```

---

## Input data

### Workflow 01 input

Place the assay-design FASTA under:

```text
data/design/
```

For example:

```text
data/design/hdv_design.fasta
```

FASTA record IDs must be unique so sequences can be tracked through orientation and circular-rotation steps.

Accepted sequence characters include standard DNA/RNA nucleotides and IUPAC ambiguity codes.

### Workflow 02 input

Place the validation FASTA under:

```text
data/validation/
```

For example:

```text
data/validation/hdv_validation.fasta
```

Workflow 02 keeps the design dataset and validation dataset conceptually separate.

For reproducible studies, document accession numbers, source database, retrieval date, search strategy, inclusion/exclusion criteria, genotype or lineage information, sequence completeness, and the analysed genomic region.

---

# Workflow 01 — Assay design

## Interactive execution

```bash
python3 01_varvamp_assay_design.py
```

Workflow 01 lists suitable FASTA files found under `data/design/`.

For detailed external commands and paths:

```bash
python3 01_varvamp_assay_design.py --verbose
```

The run-specific log is written to:

```text
results/<project>/design/<run_id>/workflow.log
```

### Project name

By default, the project name is derived from the input FASTA filename.

For example:

```text
hdv_design.fasta
```

produces the project name:

```text
hdv_design
```

To use a shorter shared project namespace for design and validation, specify it explicitly:

```bash
--project-name hdv
```

---

## Step 1 — Sequence orientation

The workflow allows three strategies:

```text
1. Keep the original orientation
2. MAFFT --adjustdirection
3. MAFFT --adjustdirectionaccurately
```

Orientation normalization statistics are reported when MAFFT orientation correction is used.

---

## Step 2 — Sequence topology

The workflow asks whether sequences are `linear` or `circular`.

For complete circular datasets, MARS can normalize cyclic start positions before the final alignment.

MARS changes the cyclic start position of a sequence; it does not change the biological nucleotide content.

---

## Step 3 — Redundancy handling

Supported strategies are:

```text
none
SeqKit exact deduplication
CD-HIT-EST clustering
```

For CD-HIT-EST, the user can control:

- identity threshold;

- clustering mode;

- strand comparison.

The supported CD-HIT-EST identity range is:

```text
0.80 to 1.00
```

A compatible word size is selected automatically.

---

## Step 4 — Final MAFFT alignment

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

The workflow records the requested strategy and, when detectable, the strategy reported by MAFFT.

---

## Step 5 — Alignment QC and trimAl

The workflow calculates:

- alignment length;

- mean gap content;

- columns with high gap fractions;

- fully occupied columns;

- sequence occupancy;

- low-occupancy sequence counts.

A transparent heuristic assessment is displayed before trimming.

Available trimAl strategies are:

```text
-noallgaps
-gappyout
-automated1
-gt VALUE
no trimming
```

---

## Conservation analysis

For each alignment position, Workflow 01 calculates:

- A/C/G/T counts;

- gap count;

- ambiguous-character count;

- occupancy;

- major base;

- major-base frequency;

- Shannon entropy;

- strict conservation;

- threshold-based conservation.

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

---

## VarVAMP assay design

Workflow 01 supports all three VarVAMP modes.

### SINGLE

Conventional PCR / individual amplicon design.

### QPCR

Primer-pair plus internal-probe design.

Main parameters include:

```text
-t     VarVAMP consensus threshold
-a     maximum ambiguous positions per primer
-pa    maximum ambiguous positions in the probe
```

### TILED

Overlapping tiled-amplicon design.

If VarVAMP fails to produce a scheme during an interactive run, the workflow can adjust VarVAMP parameters and retry **only the VarVAMP stage**. Earlier preprocessing, alignment, trimming, and conservation steps are reused.

Custom VarVAMP configuration files created during a run are stored with the run for provenance.

---

## Workflow 01 output organization

For a project called `hdv`, a run may create:

```text
work/
└── hdv/
    └── design/
        └── 20260906_164500/
            ├── preprocessing/
            ├── alignment/
            └── varvamp/
                ├── attempt_01_qpcr/
                └── ...
results/
└── hdv/
    └── design/
        └── 20260906_164500/
            ├── preprocessing/
            ├── alignment/
            ├── conservation/
            ├── config/
            ├── varvamp/
            │   ├── assay_mode.txt
            │   ├── assay_design_manifest.json
            │   └── <VarVAMP outputs>
            ├── assay_design_manifest.json
            ├── workflow_summary.txt
            └── workflow.log
```

**###** `**work/**`

Contains reconstructible technical/intermediate files, including orientation files, MARS inputs, redundancy-processing files, MAFFT working alignments, and VarVAMP attempt directories.

**###** `**results/**`

Contains permanent scientific outputs and provenance records.

Typical result groups are:

```text
preprocessing/
alignment/
conservation/
config/
varvamp/
```

---

## Assay-design manifest

At the end of Workflow 01, the script generates:

```text
results/<project>/design/<run_id>/assay_design_manifest.json
```

A copy is also placed inside the successful `varvamp/` result directory.

The manifest uses schema version 2 and records information including:

- workflow stage;

- project name;

- run ID;

- design FASTA path;

- SHA-256 of the design FASTA;

- work and results directories;

- sequence counts;

- orientation strategy;

- topology;

- redundancy strategy;

- MAFFT strategy;

- trimAl strategy;

- conservation thresholds;

- VarVAMP mode and parameters;

- VarVAMP attempts;

- VarVAMP result directory;

- software versions;

- workflow log.

Workflow 02 uses the run-level manifest as the reproducible handoff from assay design to validation.

### Manifest path portability

For files located inside the repository, Workflow 01 stores manifest paths relative to the repository root, for example:

```text
results/hdv/design/20260906_164500/varvamp
work/hdv/design/20260906_164500
data/design/hdv_design.fasta
```

This avoids embedding execution-specific prefixes such as `/app/...` from Docker or `/home/<user>/...` from a native installation.

A truly external input path remains absolute because it cannot be represented safely relative to the project.

The manifest records the path policy with:

```json
"path_base": "repository_root",
"path_policy": "repository-relative when inside the project; absolute only for external paths"
```

Workflow 02 also contains a compatibility fallback for older Workflow 01 manifests that stored absolute Docker paths such as `/app/results/...`. If the equivalent project-relative path exists in the current repository, it is remapped automatically.


---

## Non-interactive execution of Workflow 01

Example for a qPCR design:

```bash
python3 01_varvamp_assay_design.py \
  --input data/design/hdv_design.fasta \
  --project-name hdv \
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

An optional run ID can be supplied:

```bash
--run-id hdv_qpcr_test
```

Custom exact work/results directories can also be supplied with:

```bash
--workdir
--results
```

---

# Workflow 02 — In silico validation

## Interactive execution

Run:

```bash
python3 02_in_silico_validation.py
```

Detailed mode:

```bash
python3 02_in_silico_validation.py --verbose
```

Workflow 02 first searches for run-level manifests under:

```text
results/<project>/design/<run_id>/assay_design_manifest.json
```

It displays the available completed design runs and reads the selected VarVAMP mode and result directory from the manifest.

Workflow 02 then selects the validation FASTA, normally from:

```text
data/validation/
```

---

## Explicit design and validation inputs

A specific Workflow 01 run and validation database can be supplied directly:

```bash
python3 02_in_silico_validation.py \
  --design-manifest results/hdv/design/<design_run_id>/assay_design_manifest.json \
  --validation-db data/validation/hdv_validation.fasta
```

These options preselect the two main input sources. The remaining scientific validation decisions remain interactive.

An optional validation run ID can also be supplied:

```bash
--run-id validation_test
```

---

## Validation database preparation

Workflow 02 first summarizes the selected validation FASTA.

The user can then choose to:

```text
1. Keep all sequences
2. Keep sequences above a relative fraction of the maximum sequence length
3. Apply a custom minimum sequence length
```

The relative-length option is an operational sequence-length criterion. It does not prove biological completeness.

If filtering is applied, the retained FASTA is written under the validation run's `work/` directory.

---

## MFEprimer database formatting

Only the final retained validation set is formatted for MFEprimer.

SeqKit is used to convert the FASTA to one sequence line per record:

```bash
seqkit seq -w 0
```

Workflow 02 verifies that:

- complete FASTA headers are preserved;

- sequence order is preserved;

- sequence content is preserved;

- IUPAC characters are preserved.

The formatted FASTA and MFEprimer index files are technical/reconstructible files and therefore remain under:

```text
work/<project>/validation/<run_id>/database/
```

If the formatted FASTA content changes inside a reused run directory, recognized adjacent MFEprimer indexes are invalidated rather than silently reused against changed content.

---

## MFEprimer validation

Workflow 02 runs the full MFEprimer command for each selected primer pair:

```text
mfeprimer -i <primers.fasta> -d <database.fasta> -o <output>
```

The workflow expects the legacy text report to contain both:

```text
Hairpin List
Dimer List
```

If these expected QC sections are missing, the workflow stops before downstream product filtering.

For qPCR assays, this MFEprimer secondary-structure QC concerns the **LEFT and RIGHT primers** supplied to MFEprimer. The qPCR probe is not included in this hairpin/dimer check.

---

## Potential-product classification

Before downstream filtering, Workflow 02 displays all potential products reported by MFEprimer.

Products are classified as:

```text
TARGET
SELF_PRIMING
CROSS_SCHEME
OTHER
```

Definitions:

```text
TARGET
LEFT×RIGHT or RIGHT×LEFT from the selected scheme
SELF_PRIMING
LEFT×LEFT or RIGHT×RIGHT from the selected scheme
CROSS_SCHEME
recognized primer roles from different schemes
OTHER
a product that cannot be safely interpreted as belonging to the selected assay
```

The workflow then asks:

```text
Apply the TARGET-only filter before coverage and probe analysis? [y/n]:
```

If `yes`, only `TARGET` products are retained.

If `no`, `TARGET + SELF_PRIMING` products are retained.

`CROSS_SCHEME` and `OTHER` are never treated as valid products for the selected assay.

Class-specific unique HitID counts are not necessarily additive because the same biological sequence can generate more than one product class.

---

## PCR target coverage

The biological denominator is the **final retained validation database after sequence-length filtering**.

PCR coverage is calculated from unique HitIDs:

```text
PCR coverage
=
unique retained PCR-positive HitIDs
/
number of sequences in the final retained validation database
```

This avoids counting the same biological target multiple times when MFEprimer predicts multiple retained amplicons for one HitID.

---

## QPCR probe validation

For QPCR designs, Workflow 02 can additionally validate the existing VarVAMP probe with local BLAST+.

The prompt is:

```text
Test the existing VarVAMP probe(s) on the retained PCR products? [y/n]:
```

Probe validation is performed only on retained PCR products.

---

## IUPAC-degenerate probes

If a VarVAMP probe contains ambiguous IUPAC positions, Workflow 02 expands the **existing** probe into every compatible concrete A/C/G/T sequence before BLAST.

For example:

```text
Y = C/T
K = G/T
```

A probe containing one `Y` and one `K` generates four concrete probe variants.

No new ambiguity is proposed and the original probe is not redesigned.

---

## Probe BLAST

Concrete probe variants are searched with:

```text
blastn-short
```

against the retained MFEprimer amplicons.

The reported categories are:

```text
0_MISMATCH
1_MISMATCH
2_MISMATCHES
GE3_MISMATCHES
PARTIAL_ONLY
NO_HIT
```

A **full-length ungapped** probe site requires:

- query start = 1;

- query end = full probe length;

- alignment length = full probe length;

- zero gap openings.

`Full-site` therefore means a full-length ungapped BLAST alignment and may still contain three or more mismatches. It is a descriptive sequence-match category, not a prediction of wet-lab probe performance.

---

## Per-probe-variant analysis

Every concrete probe sequence is evaluated independently.

Workflow 02 reports technical counts among retained amplicons and biological counts at unique-PCR-positive-HitID level.

The files include:

```text
probe_variant_statistics.tsv
probe_variant_hitid_matrix.tsv
```

One target may match more than one concrete probe variant. Therefore per-variant rows must not be added together to calculate global biological coverage.

For final qPCR coverage, the best result across all concrete variants is retained once per unique HitID.

---

## Final qPCR coverage

The denominator is the **final retained validation database after any sequence-length filtering**.

For a final validation set containing `N` sequences:

```text
PCR
= unique retained PCR-positive HitIDs / N
PCR + probe exact
= PCR-positive HitIDs whose best concrete probe variant has 0 mismatches / N
PCR + probe <=1 MM
= PCR-positive HitIDs whose best concrete probe variant has 0 or 1 mismatch / N
PCR + probe <=2 MM
= PCR-positive HitIDs whose best concrete probe variant has 0, 1, or 2 mismatches / N
PCR + full probe site
= PCR-positive HitIDs with a full-length ungapped probe alignment / N
```

These are sequence-level computational coverage metrics. They do not establish wet-lab sensitivity or efficiency.

---

## Workflow 02 output organization

For a project called `hdv`, a validation run may create:

```text
work/
└── hdv/
    └── validation/
        └── 20260906_171000/
            ├── assay_inputs/
            │   ├── primer_pairs/
            │   ├── probes/
            │   ├── schemes/
            │   └── manifest.tsv
            │
            └── database/
                ├── <filtered_database>.fasta
                ├── <formatted_database>_fixed.fasta
                └── <MFEprimer index files>
results/
└── hdv/
    └── validation/
        └── 20260906_171000/
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

Probe-specific files are generated only for QPCR runs in which probe validation is performed.

---

## Validation manifest

Workflow 02 writes:

```text
results/<project>/validation/<run_id>/validation_manifest.json
```

The manifest uses schema version 2 and records information including:

- workflow stage;

- project name;

- validation run ID;

- parent design run ID;

- parent assay-design manifest;

- VarVAMP mode and result directory;

- design FASTA and design SHA-256 from Workflow 01;

- source validation FASTA;

- validation FASTA SHA-256;

- raw database statistics;

- sequence-length filter;

- retained validation FASTA and statistics;

- MFEprimer-formatted FASTA;

- MFEprimer index status and files;

- selected assays;

- MFEprimer parameters;

- coverage denominator;

- TARGET-only filtering decision;

- probe-validation decision;

- summary files;

- tool versions;

- validation log.

This provides an explicit parent-child provenance relationship between the assay-design run and its validation run.

As with Workflow 01, paths located inside the repository are written relative to the repository root whenever possible. This keeps validation provenance portable between Docker, native execution, and a repository moved to another filesystem location.


---

## Reproducibility

For reproducible analysis, retain or report:

1. design sequence database and retrieval date;

2. design sequence accession numbers;

3. inclusion and exclusion criteria;

4. genotype / lineage information;

5. sequence completeness;

6. orientation strategy;

7. sequence topology;

8. MARS use, when applicable;

9. redundancy strategy;

10. CD-HIT-EST parameters, when applicable;

11. MAFFT strategy;

12. alignment QC;

13. trimAl strategy and threshold, when applicable;

14. conservation-analysis thresholds;

15. VarVAMP mode;

16. VarVAMP consensus threshold;

17. primer ambiguity limit;

18. probe ambiguity limit for QPCR;

19. mode-specific VarVAMP parameters;

20. software versions;

21. validation sequence database and retrieval date;

22. validation FASTA SHA-256;

23. sequence-length filtering strategy;

24. MFEprimer indexing and search parameters;

25. TARGET-only versus TARGET + SELF_PRIMING decision;

26. PCR coverage on the final retained validation set;

27. qPCR probe BLAST parameters and mismatch statistics;

28. concrete IUPAC probe variants;

29. parent design run ID and validation run ID;

30. experimental validation of the selected assay.

Workflow 01 records many of these settings in:

```text
assay_design_manifest.json
workflow_summary.txt
workflow.log
```

Workflow 02 records its validation provenance in:

```text
validation_manifest.json
validation.log
inputs/
summary/
```

---

## Important interpretation notes

**Cluster representatives:** if clustering is used for design, validation should ideally be assessed against the original or another appropriate validation dataset rather than only against cluster representatives.

**Circular genomes:** arbitrary FASTA start positions can affect linearized analyses. MARS normalization helps standardize cyclic starts before alignment and assay design.

**Design versus validation data:** the design FASTA and validation FASTA have different roles. Independent validation data provide stronger evidence of generalizability than reusing only the design dataset.

**Primer/probe mismatches:** computational mismatch counts are descriptive. Experimental impact depends on mismatch position, sequence context, oligonucleotide chemistry, melting temperature, reaction conditions, and other factors.

**Full-site probe matches:** a full-length ungapped BLAST site can still contain mismatches and should not be interpreted as proof of efficient probe hybridization.

**MFEprimer secondary structures:** the full-QC hairpin/dimer sections apply to the LEFT/RIGHT primers supplied to MFEprimer. The VarVAMP qPCR probe is evaluated separately by BLAST for sequence matching and is not subjected to a dedicated probe hairpin/dimer calculation in this workflow.

**MFEprimer report parser:** amplicon extraction currently depends on the legacy MFEprimer 3.x text-report structure. Do not infer full MFEprimer 4.x compatibility from index-file recognition.

---

## Useful commands

Validate the Compose configuration:

```bash
docker compose config
```

Build or rebuild the Docker image:

```bash
docker compose build
```

Display Workflow 01 options:

```bash
python3 01_varvamp_assay_design.py --help
```

Display Workflow 02 options:

```bash
python3 02_in_silico_validation.py --help
```

Inspect the generated project tree:

```bash
tree data work results
```

Check Python syntax before committing changes:

```bash
python -m py_compile 01_varvamp_assay_design.py
python -m py_compile 02_in_silico_validation.py
```

---

## Third-party software

The source repository does not vendor third-party executable binaries.

`environment.yml` installs Conda/Bioconda dependencies, while the Dockerfile downloads or builds MARS and MFEprimer from their upstream sources during image construction.

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

Before redistributing a built Docker image, review the licenses and redistribution terms of the third-party software included in that image.


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