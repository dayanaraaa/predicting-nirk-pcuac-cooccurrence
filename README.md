# Predicting NirK-PCuAC Co-occurrence from Protein Sequence Features

## Overview

This repo contains the bioinformatics and machine learning pipeline built for my Master of Data Science dissertation. The project investigates whether sequence-derived protein features can predict the co-occurrence of NirK, a copper-containing nitrite reductase, and PCuAC, a periplasmic copper chaperone. Both proteins are involved in copper-dependent denitrification during anaerobic bacterial respiration.

## Research Question

**Can sequence-derived features from one protein predict whether its partner protein is also present in the same bacterial species?**

The question is framed as two conditional binary classification problems:

- **Model A**:  among species containing NirK, do NirK sequence features predict whether PCuAC is also present? (`nirK_only` vs `both`)
- **Model B**:  among species containing PCuAC, do PCuAC sequence features predict whether NirK is also present? (`pcuac_only` vs `both`)

Framing the problem conditionally means each model asks a biologically meaningful question about a single protein family, rather than mixing unrelated comparisons into one classifier.

**Unit of analysis:** one row = one bacterial species. Features are calculated once per species from a single representative sequence per protein family.

---

## Pipeline

### Stage 1:  Biological context and research question

Define the co-occurrence question and the conditional Model A / Model B framing.

### Stage 2:  Data collection

*Scripts: 2.1a-d, 2.2*

- Define seed proteins for NirK and PCuAC
- Prepare protein database
- Run BLASTp searches against RefSeq to retrieve homologous sequences
- Collect raw candidate sequences

### Stage 3:  Sequence validation and label construction

*Scripts: 3.1 - 3.5b*

- `3.1`:  parse BLAST output into clean accession lists and metadata tables
- `3.2`:  retrieve protein FASTA sequences and metadata
- `3.3`:  align sequences and tag conserved copper-binding residues
- `3.4`:  confirm NirK identity (Type 1 copper site, Type 2 copper site, catalytic residues, length)
- `3.5`:  confirm functional PCuAC identity (Cu(I) motif, length, SignalP-predicted signal peptide)

Validation is motif-based rather than relying on BLAST E-value alone, so that only sequences with the structural features required for function are retained.

### Stage 4:  Dataset construction and feature engineering

*Scripts: 4.0 - 4.2*

- `4.0`:  build sequence-level and species-level tables; select one representative sequence per species per protein family (longest validated sequence, ties resolved by accession order)
- `4.1`:  validate dataset integrity; drop taxonomically unresolved entries from the ML dataset while retaining them for QC
- `4.2`:  calculate sequence features using Biopython `ProteinAnalysis`

Features:

- Protein length
- Molecular weight
- Isoelectric point
- Hydrophobicity (GRAVY)
- Aromaticity
- Amino acid composition

Output: `out/final_feature_table.csv`:  one row per species.

### Stage 5:  Exploratory analysis

*Scripts: 5.1, 5.2, 5.2b*

- `5.1_eda.R`:  dataset structure, class composition, missingness, duplication, feature distributions, and Model A / Model B population sizes. Performed before any modelling or train/test splitting.
- `5.2_genus_structure.py` / `5.2b`:  check whether the signal in the focus feature is driven by taxonomy rather than by co-occurrence itself, using within-genus comparisons and a permutation-based variance null.

The genus-structure check matters because closely related species share both sequence features and gene content by descent. A feature that separates the classes across the whole dataset may only be tracking which genera happen to be in which class.

### Stage 6:  Modelling

*Notebook: `model_A_final.ipynb`*

- Non-parametric group comparisons (Mann-Whitney U) with effect sizes
- Correlation filtering to remove redundant features
- Train/test split
- Logistic regression with L2 regularisation, `C` selected by cross-validation (standardised features)
- Random forest (unscaled features)

### Stage 7:  Evaluation and interpretation

- Accuracy, precision, recall, F1, confusion matrix
- Logistic regression coefficients (standardised units)
- Random forest feature importance

Statistical testing indicates whether a feature *differs* between groups. It does not establish that the feature can *predict* class membership:  the modelling stage addresses that separately.

---

## Tools

**Bioinformatics**
BLASTp · NCBI RefSeq protein database · MAFFT · SignalP · Biopython

**Data science**
Python (pandas, NumPy, scikit-learn, SciPy, matplotlib) · R (tidyverse)

**Computing**
University HPC cluster with SLURM job scheduling

---

## Data availability

Raw BLAST output and downloaded RefSeq sequences are not included in this repo due to size. Scripts in `scripts/2_data_collection/` will regenerate them, though results may shift slightly as source databases are updated. Processed and intermediate data will be added as later stages are finalised.

---
## Future Research Directions
This project identified a measurable difference between the groups and highlighted sequence features that may contribute to this distinction. Several follow-up analyses could build on these findings to determine whether these features represent biologically meaningful signals and to investigate their potential functional consequences.

### 1. Validate sequence-level signals with ESM-2

A natural next step would be to use **ESM-2** to investigate whether a protein language model independently identifies the same regions or sequence features highlighted by the current analysis.

**Key question:**

> Do the features that distinguish the groups in the current analysis also emerge from an independent protein language model?

Agreement between the approaches would provide additional support that these regions contain meaningful sequence information. Conversely, differences between the approaches could identify additional regions or suggest that the original signal may require further investigation.

### 2. Account for phylogenetic relationships

The observed differences may partly reflect **shared evolutionary history** rather than the biological distinction being investigated.

A phylogenetic analysis could therefore assess whether the identified sequence patterns remain associated with the groups after accounting for evolutionary relatedness.

**Key question:**

> Are the observed differences genuinely associated with the biological variable of interest, or can they be explained by phylogenetic structure?

This would help distinguish biologically relevant signals from patterns arising through common ancestry.

### 3. Investigate functional consequences through interaction modelling

If particular sequence regions consistently emerge as distinguishing features, a further question is whether these differences have consequences at the **protein-protein interaction** level.

Interaction modelling could investigate whether the identified sequence differences alter predicted interaction interfaces, binding characteristics, or protein-protein interaction networks.

**Key question:**

> Do the sequence differences identified in this study have potential functional consequences for protein interactions?

This would provide a route from identifying **sequence-level differences** towards understanding their possible **molecular mechanisms**.

### Overall research trajectory

The proposed progression is:

**Identify distinguishing features**
↓
**Test whether ESM-2 identifies the same signals**
↓
**Account for phylogenetic history**
↓
**Investigate potential functional consequences through interaction modelling**

Together, these directions could extend the current work from identifying **what differs** towards understanding **why it differs and whether those differences matter functionally**.

## Limitations

- **Representative sequence choice.** One sequence per species per protein family is used. Where a species carries multiple validated homologues, the longest is chosen, which discards within-species variation.
- **Absence is inferred, not observed.** A species labelled `nirK_only` is one where no PCuAC homologue passed validation. Absence from a search result is not the same as absence from the genome, and incomplete or unsequenced genomes will inflate the "only" classes.
- **Taxonomic non-independence.** Species are not statistically independent observations; shared ancestry links both features and gene content. Stage 5.2 assesses this but does not fully remove it.
- **Class imbalance** in both Model A and Model B populations.
- **No structural or experimental validation.** This project is entirely computational.
