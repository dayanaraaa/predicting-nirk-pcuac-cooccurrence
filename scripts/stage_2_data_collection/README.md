# Stage 2: Data Collection

## Overview
This stage is about finding candidate NirK and PCuAC homologues by running large-scale BLASTp searches against the NCBI RefSeq Protein database.

I start with curated seed sequences for NirK and PCuAC, then use those as queries to search across bacterial protein sequences for potential matches. The raw BLASTp hits that come out of this stage get passed along to Stage 3, where they're cleaned up and validated properly.

---
## Objectives
- Define seed protein sequences for NirK and PCuAC
- Set up the RefSeq Protein database for local BLASTp searches
- Run the BLASTp searches at scale
- Produce raw candidate hit tables for cleaning and validation in the next stage

---
## Workflow
Seed protein sequences >RefSeq Protein database preparation > BLASTp searches > Raw BLASTp hit tables > Stage 3: Cleaning and Label Construction

---
## Scripts

| Script | Purpose |
|---|---|
| `2.1a_*` | Define/retrieve the NirK seed protein sequence |
| `2.1b_*` | Define/retrieve the PCuAC seed protein sequence |
| `2.1c_*` | Run BLASTp searches against the RefSeq Protein database |
| `2.1d_*` | Process and extract the raw BLASTp hits |
| `2.2_*` | Supporting data prep and database handling |

---
## Inputs
- NirK seed protein sequence
- PCuAC seed protein sequence
- NCBI RefSeq Protein database
- BLAST+ tools

---
## Outputs
Raw BLASTp results with candidate homologous hits, including:
- Protein accession identifiers
- Alignment statistics
- Sequence similarity metrics
- Candidate NirK and PCuAC matches

These get passed to Stage 3 for cleaning, metadata retrieval, sequence validation, and species-level label construction.

---
## Computational Environment
This stage runs on a university HPC cluster using SLURM job submission, the BLASTp searches are too computationally heavy to run locally at this scale.

Example submission:
```bash
sbatch 2.1c_blastp_search.sh
```

---
## Notes 
A couple of things worth flagging for this stage specifically:
- BLASTp searches against the full RefSeq Protein database can take a while depending on queue load, so I'd recommend checking walltime limits before submitting rather than after a job gets killed partway through.
- Seed sequences matter a lot here, if the NirK/PCuAC seeds aren't representative enough, the search may miss legitimate homologues or pull in too many false positives. Worth double-checking these against known reference sequences before running at scale.