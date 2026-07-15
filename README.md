# Predicting NirK-PCuAC Co-occurrence from Protein Sequence Features
**In progress:**
Currently on Stage 3

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
> Currently here!  
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

- The feature set may get trimmed after EDA, some engineered features might turn out not to add predictive value.
- Species-level label construction involves some judgment calls around ambiguous or partial BLAST matches; these are documented in the relevant scripts but may need revisiting.
- No structural or experimental validation is planned, this project is entirely computational-based.
- Model selection, training, and evaluation are still to come, so no results or conclusions exist yet.
- Class imbalance (species with NirK but not PCuAC, or vice versa) may need specific handling once the full dataset is built, this hasn't been assessed yet.
