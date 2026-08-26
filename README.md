# Predicting NirK-PCuAC Co-occurrence from Protein Sequence Features

## Overview 
This repo includes the bioinformatics and machine learning pipeline I'm building for my Master of Data Science dissertation. The project investigates whether sequence-derived protein features can be used to predict the co-occurrence of NirK, a copper-containing nitrite reductase, and PCuAC, a periplasmic copper chaperone. These proteins are both involved in copper-dependent denitrification for anaerobic bacterial respiration.

## Research Question
**Can sequence-derived features from NirK and PCuAC proteinsb predict whether a bacterial species contains both proteins?**  
More specifically:
- Can NirK features predict whether PCuAC is present?
- Can PCuAC features predict whether NirK is present?
- Are there any sequence-level signals that could be tied to their co-occurrence?

## Pipeline Overview
### Stage 1: Biological Context and Define Research Question
Framed as a binary classification problem: whether a bacterial species contains both NirK and PCuAC proteins.

### Stage 2: Data Collection 
_Scripts: 2.1a-d, 2.2_  
- Define seed proteins for NirK and PCuAC
- Prepare protein database
- Run large-scale BLASTp searches against RefSeq to retrieve homologous sequences
- Collect raw candidate homologous sequences

### Stage 3: Data Cleaning and Label Construction
_Scripts: 3.1-3.6_
- Standardise BLAST metadata and accession identifiers
- Retrieve protein sequence metadata and FASTA sequences
- Validate candidate homologues using conserved residue analysis
- Remove duplicate, incomplete, or low-quality sequences
- Construct species-level NirK and PCuAC presence/absence labels
- Generate representative protein datasets for downstream feature analysis

### Stage 4: Feature Engineering
_Scripts: 4.1, 4.2_  
- Integrate validated protein-sequences with sequence-derived features
- Physicochemical protein descriptors of interest:
  - Protein length
  - Molecular Weight
  - Isoelectric Point
  - Hydrophobicity (GRAVY score)
  - Amino acid composition
  - Conserved residue features
- Generate final analysis table and machine-learning feature matrices

### Stage 5: Exploratory Data Analysis
_Scripts: _  
- Assess dataset composition and class balance
- Compare NirK-only, PCuAC-only, and co-occurring groups
- Identify feature redundancy through correlation analysis
- Evaluate sequence feature distributions and relationships
- Perform dimensionality reduction using Principal Component Analysis (PCA)
- Generate ESM2 protein language model embeddings to compare sequence-levels
- Use to develop machine learning models for predicting NirK-PCuAC co-occurrence

### Stage 6: Prepare Data for Modelling
- Encode categorical variables
- Scale numerical features where appropriate
- Create training, validation, and test datasets
- Address potential class imbalance

### Stage 7: Model Selection
_Currently considering..._
- Logistic Regression (baseline)
- Decision Tree
- Random Forest
- XGBoost / Gradient Boosting
- Neural Network

### Stage 8: Model Training
- Train selected models using cross-validation
- Optimise hyperparameters
- Evaluate model performance during training

### Stage 9: Model Evaluation
- Accuracy
- Precision
- F1-score
- Confusion matrix

### Stage 10: Model Interpretation
- Feature importance analysis
- Identify sequence features associated with NirK-PCuAC co-occurrence
---

## Tools and Technologies 
**Bioinformatics**
- BLASTp
- NCBI RefSeq Protein DB
- MAFFT
  
**Programming & Data Science**
- Python
- R
- Pandas
- NumPy
- scikit-learn
- matplotlib
  
**Computing Environment**
- University HPC Cluster
- SLURM job scheduling

---
**Data availability**   

Raw data pulled from public databases (RefSeq, BLASTp results) isn't included in this repo due to size. Scripts in scripts/2_data_collection/ will regenerate it, though results may shift slightly as source databases get updated. Processed and intermediate data will be added as later stages are finalised. 
 
---
**Known Limitations & Future Work**

Yes. For a GitHub README, I’d make it concise and frame it as a **research roadmap**, clearly separating what the dissertation established from what could be investigated next.

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

If particular sequence regions consistently emerge as distinguishing features, a further question is whether these differences have consequences at the **protein–protein interaction** level.

Interaction modelling could investigate whether the identified sequence differences alter predicted interaction interfaces, binding characteristics, or protein–protein interaction networks.

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

