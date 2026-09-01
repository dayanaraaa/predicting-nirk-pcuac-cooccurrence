Datasets for an MDS 
investigating whether sequence-derived protein features predict the species-level
co-occurrence of copper nitrite reductase (NirK) and its periplasmic copper
chaperone (PCuAC) across bacteria.

## Data files

| File | Rows | Description |
|---|---|---|
| `sequence_master_table.csv` | 23,663 | Validated NirK (7,485) and PCuAC (16,178) sequences |
| `species_master_table.csv` | 11,125 | All species-level entries, incl. unresolved (QC) |
| `species_dataset.csv` | 5,815 | Resolved species only; input to feature extraction |
| `final_feature_table.csv` | 5,815 | 25 features per protein family; input to stage 5 |

Counts correspond to Table 3.1 of the report.

## Notes
- Random seed 11 throughout. Software versions in Appendix A1.
- RefSeq protein database, release 27 June 2026; BLAST+ 2.17.0.
- `species_ml_dataset` in archives is superseded table from an earlier run, kept for
  transparency. It is not read by any script and does not correspond to any
  reported result. All analyses use `species_dataset.csv`.
