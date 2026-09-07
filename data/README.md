# Input Data

This directory contains user-provided nucleotide FASTA datasets used by the two workflow stages.

The directory is divided into two biologically distinct input roles:

```text
data/
├── design/
└── validation/
```

---

## `design/`

FASTA datasets stored under `data/design/` are used by:

```text
01_varvamp_assay_design.py
```

for:

- nucleotide FASTA validation;
- sequence-orientation normalization;
- circular start-position normalization with MARS, when applicable;
- redundancy handling with SeqKit or CD-HIT-EST;
- multiple-sequence alignment with MAFFT;
- alignment quality control;
- optional trimming with trimAl;
- conservation analysis;
- assay design with VarVAMP.

Example:

```text
data/design/hdv_design.fasta
```

FASTA record IDs should be unique so sequences can be tracked reliably through preprocessing and alignment steps.

---

## `validation/`

FASTA datasets stored under `data/validation/` are used by:

```text
02_in_silico_validation.py
```

for:

- validation-database preparation;
- optional sequence-length filtering;
- MFEprimer database formatting and indexing;
- primer-pair specificity and target-coverage analysis;
- qPCR probe analysis with BLAST+, when applicable.

Example:

```text
data/validation/hdv_validation.fasta
```

The validation dataset may be larger and more diverse than the design dataset.

When possible, an independent validation dataset is preferable because it provides stronger evidence of assay generalizability than validation against only the sequences used during assay design.

---

## Design and validation datasets

The two directories have different roles:

```text
data/design/       sequences used to design the assay

data/validation/   sequences used to evaluate the designed assay
```

The same FASTA may technically be used for both stages, but independent or broader validation data are recommended when scientifically appropriate.

For reproducible studies, document relevant metadata such as:

- sequence database or source;
- retrieval date;
- accession numbers;
- search or download strategy;
- inclusion and exclusion criteria;
- genotype, lineage, or taxonomic information;
- sequence completeness;
- analysed genomic region.

---

## Docker usage

The default Docker configuration mounts the repository at:

```text
/app
```

Therefore project input datasets are available inside the container as:

```text
/app/data/design/
/app/data/validation/
```

For normal Docker execution, place biological FASTA inputs inside these directories before running:

```bash
docker compose run --rm design
```

or:

```bash
docker compose run --rm validation
```

Arbitrary files outside the repository are not mounted into the container by default.

---

## Git policy

Biological FASTA datasets are ignored by Git by default.

Large or study-specific sequence datasets should not normally be committed to the repository.

The directory structure can be preserved using files such as:

```text
data/design/.gitkeep
data/validation/.gitkeep
```

Only small, public, non-confidential, and appropriately licensed example datasets should be committed intentionally.

Do not commit:

- confidential sequence data;
- restricted datasets;
- unpublished data without authorization;
- personally identifiable information;
- credentials or private access information.