# Input data

This directory is intended for local nucleotide input datasets and reproducibility metadata.

The repository can use FASTA files stored here, or files located elsewhere on the computer.

Example:

```text
data/my_sequences.fasta
```

Raw sequence files are ignored by the default `.gitignore` to reduce the risk of accidentally publishing:

- large datasets;
- restricted datasets;
- confidential data;
- unpublished sequence data;
- insufficiently documented input files.

---

## Workflow 01 input

The primary input for:

```text
01_varvamp_assay_design.py
```

is a nucleotide FASTA file.

Example:

```text
data/hdv_complete_genomes.fasta
```

The input should contain nucleotide sequences only.

Accepted nucleotide symbols include standard DNA/RNA bases and IUPAC ambiguity codes.

Protein FASTA files are not supported.

FASTA record IDs should be unique so that sequences can be traced through:

- orientation correction;
- circular start-position normalization;
- redundancy reduction;
- multiple sequence alignment;
- assay design.

---

## Sequence topology

Workflow 01 supports both:

```text
linear
```

and:

```text
circular
```

sequence datasets.

For complete circular genomes, arbitrary FASTA start positions may affect linearized alignment and downstream assay coordinates.

When appropriate, Workflow 01 can use MARS to normalize cyclic start positions before the final multiple sequence alignment.

The topology selected during analysis should be documented for reproducibility.

---

## Workflow 02 validation data

The validation workflow:

```text
02_in_silico_validation.py
```

uses a nucleotide FASTA database against which designed primer pairs and qPCR probes are evaluated.

The validation FASTA may be:

- the original unclustered design dataset;
- a broader sequence collection;
- an independently retrieved validation dataset;
- another biologically appropriate target database.

Example:

```text
data/hdv_validation_sequences.fasta
```

For stronger evidence of assay generalizability, an independently assembled validation dataset is preferable to relying only on the exact dataset used during assay design.

---

## Clustered versus original datasets

If CD-HIT-EST or SeqKit is used during Workflow 01, the reduced sequence set is useful for assay design but should not automatically be used as the only validation denominator.

Whenever appropriate, evaluate assay coverage against:

```text
the original complete dataset
```

rather than only:

```text
cluster representatives
```

This helps avoid overestimating biological target coverage.

---

## Reproducibility metadata

For a reproducible public repository, provide an:

```text
accession_numbers.tsv
```

file containing accession and dataset metadata.

A minimal structure may be:

```text
accession	genotype_or_lineage	database	retrieval_date
```

Example:

```text
OR613551.1	HDV	genbank	2026-08-01
OQ200691.1	HDV	genbank	2026-08-01
```

Additional columns may be included when useful, for example:

```text
accession
genotype_or_lineage
subtype
host
country
collection_date
database
retrieval_date
sequence_length
complete_genome
included
exclusion_reason
```

---

## Dataset documentation

For every dataset used for assay design or validation, document as much as possible:

- source database;
- exact retrieval date;
- database query or search terms;
- accession numbers;
- inclusion criteria;
- exclusion criteria;
- genotype, lineage, strain, or subtype information;
- sequence completeness criteria;
- genomic region analysed;
- treatment of partial sequences;
- treatment of exact duplicates;
- clustering or dereplication strategy;
- sequence-orientation handling;
- sequence topology;
- circular-genome normalization, when applicable;
- metadata used for stratified validation;
- any manual filtering performed before analysis.

---

## Design and validation datasets

When the same sequence collection is used for both assay design and validation, results should be interpreted as:

```text
internal target coverage
```

rather than fully independent validation.

When possible, maintain separate records for:

```text
design dataset
validation dataset
```

and document how each dataset was obtained.

---

## Example directory

A local project may use:

```text
data/
├── README.md
├── accession_numbers.tsv
├── design_sequences.fasta
└── validation_sequences.fasta
```

With the default `.gitignore`, the FASTA files remain local while the documentation and accession table can be tracked by Git.

---

## Data privacy and redistribution

Only redistribute sequence data when:

- the source database permits redistribution;
- institutional policies permit redistribution;
- collaborator agreements permit redistribution;
- the data are not confidential or restricted.

Do not commit publicly:

- personally identifiable information;
- clinical identifiers;
- confidential sequence data;
- restricted-access genomic data;
- unpublished collaborator datasets without permission;
- credentials;
- private institutional paths.

When raw sequence redistribution is not appropriate, provide reproducibility metadata such as:

```text
accession_numbers.tsv
```

and describe how users can retrieve the source sequences independently.

---

## Git tracking

The default `.gitignore` is intended to keep local sequence files out of Git, for example:

```gitignore
data/*.fa
data/*.fasta
data/*.fna
data/*.fas
data/*.fastq
data/*.fq
data/*.fastq.gz
data/*.fq.gz
```

while preserving:

```gitignore
!data/README.md
!data/accession_numbers.tsv
```

This keeps the repository lightweight while retaining the information needed to reconstruct the datasets used in an analysis.