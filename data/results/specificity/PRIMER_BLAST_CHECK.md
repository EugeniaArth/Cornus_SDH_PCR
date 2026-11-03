# Primer-BLAST-style local check (blastn-short vs the six Cornus DQD/SDH transcripts)

Parameters: `-task blastn-short -word_size 7 -evalue 1000 -perc_identity 75`.  
Product window: 80-600 bp, total mismatches ≤ 3, perfect 5-nt 3'-terminal match required for each primer.
Database: `DQD_SDH_possible_nt.fasta` (6 transcripts, 12,888 bp).  
Note: this is a LOCAL approximation of Primer-BLAST specificity. It does not scan the whole transcriptome yet (see note at bottom).

## Verdicts

| Pair | Verdict | Target product(s) | Off-target product(s) |
|---|---|---|---|
| CserDQD_SDH1.1 | SPECIFIC | 120 bp (0 mm) | none |
| CserDQD_SDH1.2 | SPECIFIC | 180 bp (0 mm) | none |
| CserDQD_SDH2 | SPECIFIC | 240 bp (0 mm) | none |
| CserDQD_SDH3 | SPECIFIC | 300 bp (0 mm) | none |
| CserDQD_SDH4 | SPECIFIC | 360 bp (0 mm) | none |
| CserDQD_SDH5 | SPECIFIC | 420 bp (0 mm) | none |

## Priming sites per primer (single-primer hits, any identity ≥75 %)

| Primer | Transcript | Hits |
|---|---|---|
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.1 | 4 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.2 | 4 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH2 | 4 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH3 | 6 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH4 | 5 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH5 | 3 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.1 | 5 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.2 | 9 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH2 | 4 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH3 | 2 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH4 | 5 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH5 | 4 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.1 | 5 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.2 | 7 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH2 | 7 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH3 | 3 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH4 | 2 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH5 | 2 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.1 | 2 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.2 | 7 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH2 | 5 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH3 | 7 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH4 | 4 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH5 | 3 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.1 | 4 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.2 | 5 |
| CserDQD_SDH2_Fwd | CserDQD/SDH2 | 6 |
| CserDQD_SDH2_Fwd | CserDQD/SDH3 | 2 |
| CserDQD_SDH2_Fwd | CserDQD/SDH4 | 2 |
| CserDQD_SDH2_Fwd | CserDQD/SDH5 | 3 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.1 | 7 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.2 | 10 |
| CserDQD_SDH2_Rev | CserDQD/SDH2 | 2 |
| CserDQD_SDH2_Rev | CserDQD/SDH3 | 2 |
| CserDQD_SDH2_Rev | CserDQD/SDH4 | 4 |
| CserDQD_SDH2_Rev | CserDQD/SDH5 | 5 |
| CserDQD_SDH3_Fwd | CserDQD/SDH1.2 | 1 |
| CserDQD_SDH3_Fwd | CserDQD/SDH2 | 6 |
| CserDQD_SDH3_Fwd | CserDQD/SDH3 | 5 |
| CserDQD_SDH3_Fwd | CserDQD/SDH4 | 2 |
| CserDQD_SDH3_Fwd | CserDQD/SDH5 | 5 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.1 | 3 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.2 | 4 |
| CserDQD_SDH3_Rev | CserDQD/SDH2 | 5 |
| CserDQD_SDH3_Rev | CserDQD/SDH3 | 5 |
| CserDQD_SDH3_Rev | CserDQD/SDH4 | 5 |
| CserDQD_SDH3_Rev | CserDQD/SDH5 | 2 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.1 | 2 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.2 | 3 |
| CserDQD_SDH4_Fwd | CserDQD/SDH2 | 5 |
| CserDQD_SDH4_Fwd | CserDQD/SDH3 | 2 |
| CserDQD_SDH4_Fwd | CserDQD/SDH4 | 4 |
| CserDQD_SDH4_Fwd | CserDQD/SDH5 | 3 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.1 | 3 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.2 | 4 |
| CserDQD_SDH4_Rev | CserDQD/SDH2 | 5 |
| CserDQD_SDH4_Rev | CserDQD/SDH3 | 7 |
| CserDQD_SDH4_Rev | CserDQD/SDH4 | 2 |
| CserDQD_SDH4_Rev | CserDQD/SDH5 | 2 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.1 | 3 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.2 | 3 |
| CserDQD_SDH5_Fwd | CserDQD/SDH2 | 4 |
| CserDQD_SDH5_Fwd | CserDQD/SDH3 | 6 |
| CserDQD_SDH5_Fwd | CserDQD/SDH4 | 3 |
| CserDQD_SDH5_Fwd | CserDQD/SDH5 | 6 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.1 | 3 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.2 | 3 |
| CserDQD_SDH5_Rev | CserDQD/SDH3 | 2 |
| CserDQD_SDH5_Rev | CserDQD/SDH4 | 5 |
| CserDQD_SDH5_Rev | CserDQD/SDH5 | 6 |

## Full candidate products (allowing mismatches)

| Pair | Transcript | Size | mm(F) | mm(R) | mm total | Class |
|---|---|---|---|---|---|---|
| CserDQD_SDH1.1 | CserDQD/SDH1.1 | 120 | 0 | 0 | 0 | TARGET |
| CserDQD_SDH1.2 | CserDQD/SDH1.2 | 180 | 0 | 0 | 0 | TARGET |
| CserDQD_SDH2 | CserDQD/SDH2 | 240 | 0 | 0 | 0 | TARGET |
| CserDQD_SDH3 | CserDQD/SDH3 | 300 | 0 | 0 | 0 | TARGET |
| CserDQD_SDH4 | CserDQD/SDH4 | 360 | 0 | 0 | 0 | TARGET |
| CserDQD_SDH5 | CserDQD/SDH5 | 420 | 0 | 0 | 0 | TARGET |

## Still missing for a real Primer-BLAST

- Whole-transcriptome database (the Trinity `transcripts.fasta` in the sibling `Cornus_transcriptome` folder) is outside this project directory — it needs explicit read permission.
- The *Cornus* genome/transcriptome is not in NCBI, so the web Primer-BLAST would require pasting the six transcripts as `primer-orientation` template or uploading a custom database.
