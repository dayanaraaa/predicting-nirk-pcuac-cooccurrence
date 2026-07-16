# Stage 3: Cleaning, Validation, and Species‑Level Label Construction

## Overview
Stage 3 transforms the raw BLASTp hit tables from Stage 2 into a clean, biologically meaningful dataset of NirK and PCuAC homologues at the species level.

It includes:
- Cleaning BLASTp results
- Validating accessions
- Retrieving FASTA sequences
- Checking conserved catalytic/motif residues
- Collapsing strain‑level entries into species
- Constructing the final species of NirK and PCuAC presence/absence table

## Inputs
- Raw BLASTp hit tables from Stage 2
- RefSeq protein accession metadata
- FASTA sequences retrieved from NCBI
- Residue QC rules for NirK and PCuAC
- Taxonomy information (species, strain, genus, phylum) * step has not been added into workflow yet

## Outputs
- species_master_table
- representative NirK accession
- representative PCuAC accession
- NirK_present (0/1)
- PCuAC_present (0/1)
- class (nirK_only/pcuac_only/both/neither)

- nirK_clean.fasta
- pcuAC_clean.fasta
- Residue QC tables
- catalytic site completeness
- motif presence
- histidine tail counts

---

## Objectives
- Clean raw BLASTp hit tables from Stage 2
- Validate protein accessions and remove invalid or duplicate entries
- Retrieve full protein sequences for NirK and PCuAC candidates
- Perform residue‑level quality checks (catalytic sites, motifs, histidine tail counts)
- Collapse strain‑level entries into species‑level representatives
- Construct the final species‑level presence/absence labels
- Produce clean FASTA files for Stage 4 feature extraction

---

## Workflow

1. Raw BLAST hit cleaning
Remove duplicates, low‑quality alignments, and non‑bacterial hits.
2. Accession validation
Ensure each accession resolves to a valid RefSeq protein.
3. Sequence retrieval
Download FASTA sequences for all candidate NirK and PCuAC hits.
4. Residue-level QC- NirK: T1/T2 copper sites, catalytic histidines
- PCuAC: Cys‑His motif, histidine tail count

5. Species collapsing
Merge strain‑level entries into species‑level representatives.
6. Label construction
Build the final table of species presence: NirK_present or PCuAC_present, or both.
7. Output FASTA files
Cleaned NirK and PCuAC FASTA files for Stage 4.
8. Pass to Stage 4
Feature extraction and integration.

---

## Scripts

Script	Purpose	
3.1_*	Clean raw BLASTp hit tables (remove duplicates, filter low‑quality hits)	
3.2_*	Validate accessions and retrieve RefSeq metadata	
3.3_*	Download FASTA sequences for NirK and PCuAC candidates	
3.4_*	Perform residue‑level QC (catalytic sites, motifs, histidine counts)	
3.5_*	Collapse strain‑level entries into species representatives	
3.6_*	Construct species‑level presence/absence labels and write master table	
3.7_*	Produce cleaned FASTA files for Stage 4	


---

## Inputs

- Raw BLASTp hit tables from Stage 2
- RefSeq protein accession metadata
- FASTA sequences retrieved from NCBI
- Residue QC rules for NirK and PCuAC
- Taxonomy information (species, strain, genus, phylum)


---

## Outputs

Stage 3 produces the cleaned, validated, species‑level dataset:

- species_master_table.csv- species name
- representative NirK accession
- representative PCuAC accession
- NirK_present (0/1)
- PCuAC_present (0/1)
- class (nirK_only/pcuac_only/both/neither)

- nirK_clean.fasta
- pcuAC_clean.fasta
- Residue QC tables- catalytic site completeness
- motif presence
- histidine tail counts

---

## Computational Environment

Stage 3 runs on the same HPC environment as Stage 2:
- SLURM job submission
- Python scripts for accession validation, QC, and FASTA retrieval
- R scripts for taxonomy collapsing and label construction

---

## Notes

A few important considerations for Stage 3:  
- Species collapsing is critical.
RefSeq contains many strain‑level entries; collapsing them ensures each species contributes only one representative sequence.
- Residue QC prevents false positives.
BLASTp alone cannot guarantee functional NirK or PCuAC homologues.
Checking catalytic/motif residues removes spurious hits.
- Presence/absence labels define the ML problem.
The accuracy of Stage 3 directly determines the quality of Stage 4 and Stage 5 outputs.
- NirK is much rarer than PCuAC.
Expect ~2,377 NirK species vs ~7,563 PCuAC species - this imbalance is real biology, not a pipeline error.

