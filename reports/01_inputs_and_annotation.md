# 01 — Inputs and annotation verification

**Step scripts:** `scripts/step1_inputs/`
**Outputs:** `data/processed/genbank_features.tsv`, `data/processed/genbank_coordinates_combined.tsv`



## Material

| File | Role |
|---|---|
| `data/raw/PZ518666.gb` … `PZ518671.gb` | all six GenBank records |

(The protein accession map, alignment CSVs and the primer table from other species have been moved to `archive/`; they are not needed for the pipeline.)

Gene and accession map:

| Accession | Gene | Transcript length |
|---|---|---|
| PZ518666 | CserDQD/SDH1.1 | 1885 nt |
| PZ518667 | CserDQD/SDH1.2 | 2779 nt |
| PZ518668 | CserDQD/SDH2 | 2241 nt |
| PZ518669 | CserDQD/SDH3 | 2133 nt |
| PZ518670 | CserDQD/SDH4 | 2083 nt |
| PZ518671 | CserDQD/SDH5 | 1767 nt |

## GenBank accessions and coordinate extraction

The six *Cornus sericea* DQD/SDH gene variants are deposited in NCBI under accessions
PZ518666–PZ518671. These records contain the same sequences as `DQD_SDH_possible_nt.fasta`
and provide annotated CDS and UTR coordinates. `download_genbank.py` retrieves the six
records if they are missing. `parse_genbank.py` extracts those coordinates,
and `extract_genbank_coords.py` writes them to `data/processed/genbank_coordinates.tsv`,
which becomes the coordinate input for the rest of the pipeline.


## Reproduce

```sh
python3 scripts/step1_inputs/download_genbank.py
python3 scripts/step1_inputs/parse_genbank.py
python3 scripts/step1_inputs/extract_genbank_coords.py
```
