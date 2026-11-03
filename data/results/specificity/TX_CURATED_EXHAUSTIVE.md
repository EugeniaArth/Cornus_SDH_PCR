# Exhaustive mismatch scan on the six curated DQD/SDH transcripts

All 12 primers x 6 transcripts x 2 strands; every 20-nt window was compared (no BLAST seeding, no gaps). Every site with ≤6 mismatches is in `curated_exhaustive_sites.tsv`; the table below lists only the best window (fewest mismatches) per primer x transcript.

## Best full-length window per primer x transcript

| Primer | Transcript | mm (of 20) | 3'-perfect run | mismatch positions (primer 5'->3') |
|---|---|---|---|---|
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.1 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.1 | 7 | 0 | 20,12,10,6,5,2,1 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.2 | 2 | 0 | 17,20 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH1.2 | 7 | 0 | 20,12,10,6,5,2,1 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH2 | 8 | 2 | 4,5,6,8,11,14,15,18 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH2 | 8 | 2 | 18,16,13,10,9,8,6,4 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH3 | 7 | 2 | 4,6,8,11,14,17,18 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH3 | 7 | 7 | 13,11,9,7,6,3,1 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH4 | 8 | 4 | 1,2,4,5,6,10,14,16 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH4 | 7 | 5 | 15,12,11,9,8,6,1 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH5 | 8 | 2 | 2,5,10,11,12,16,17,18 |
| CserDQD_SDH1.1_Fwd | CserDQD/SDH5 | 8 | 2 | 18,16,15,11,9,4,3,2 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.1 | 8 | 1 | 2,3,8,9,10,13,17,19 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.1 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.2 | 8 | 1 | 2,3,8,9,10,13,17,19 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH1.2 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH1.1_Rev | CserDQD/SDH2 | 8 | 1 | 8,10,11,13,14,17,18,19 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH2 | 7 | 2 | 18,17,13,11,10,5,4 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH3 | 8 | 0 | 2,6,9,12,13,17,18,20 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH3 | 6 | 2 | 18,17,13,10,5,4 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH4 | 7 | 1 | 4,5,11,15,17,18,19 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH4 | 8 | 0 | 20,18,16,15,13,11,10,3 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH5 | 8 | 1 | 1,2,4,6,9,11,13,19 |
| CserDQD_SDH1.1_Rev | CserDQD/SDH5 | 8 | 0 | 20,19,16,13,11,10,8,4 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.1 | 8 | 3 | 2,5,7,9,12,15,16,17 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.1 | 7 | 0 | 20,18,15,14,12,9,5 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.2 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH1.2 | 7 | 0 | 20,18,15,14,12,9,5 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH2 | 8 | 0 | 4,8,9,15,17,18,19,20 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH2 | 8 | 0 | 20,11,7,6,5,4,3,2 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH3 | 7 | 2 | 2,3,7,12,15,16,18 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH3 | 8 | 0 | 20,19,17,10,7,5,2,1 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH4 | 7 | 0 | 5,6,7,12,15,16,20 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH4 | 8 | 1 | 19,18,17,16,11,10,6,5 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH5 | 7 | 0 | 2,6,11,16,17,18,20 |
| CserDQD_SDH1.2_Fwd | CserDQD/SDH5 | 7 | 0 | 20,18,14,12,11,5,1 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.1 | 8 | 0 | 1,2,3,12,13,14,19,20 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.1 | 8 | 3 | 17,16,15,14,12,9,3,1 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.2 | 8 | 0 | 1,2,3,12,13,14,19,20 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH1.2 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH1.2_Rev | CserDQD/SDH2 | 9 | 0 | 2,3,4,6,10,12,13,19,20 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH2 | 6 | 0 | 20,17,15,12,6,1 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH3 | 8 | 3 | 1,3,4,10,13,15,16,17 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH3 | 7 | 5 | 15,14,12,8,6,2,1 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH4 | 9 | 2 | 1,2,4,5,7,10,11,15,18 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH4 | 8 | 1 | 19,18,16,13,11,5,4,2 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH5 | 8 | 3 | 1,2,8,12,13,14,15,17 |
| CserDQD_SDH1.2_Rev | CserDQD/SDH5 | 7 | 1 | 19,17,15,10,9,6,1 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.1 | 7 | 2 | 3,6,7,8,14,16,18 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.1 | 9 | 7 | 13,12,9,8,7,6,4,3,1 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.2 | 7 | 2 | 3,6,7,8,14,16,18 |
| CserDQD_SDH2_Fwd | CserDQD/SDH1.2 | 7 | 0 | 20,19,18,15,14,11,4 |
| CserDQD_SDH2_Fwd | CserDQD/SDH2 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH2_Fwd | CserDQD/SDH2 | 8 | 7 | 13,12,11,7,6,4,2,1 |
| CserDQD_SDH2_Fwd | CserDQD/SDH3 | 5 | 5 | 7,8,9,10,15 |
| CserDQD_SDH2_Fwd | CserDQD/SDH3 | 9 | 2 | 18,17,15,10,8,7,6,5,4 |
| CserDQD_SDH2_Fwd | CserDQD/SDH4 | 5 | 9 | 6,8,9,10,11 |
| CserDQD_SDH2_Fwd | CserDQD/SDH4 | 8 | 3 | 17,11,7,6,5,4,2,1 |
| CserDQD_SDH2_Fwd | CserDQD/SDH5 | 8 | 1 | 2,3,8,9,15,17,18,19 |
| CserDQD_SDH2_Fwd | CserDQD/SDH5 | 9 | 2 | 18,17,14,13,8,7,6,5,4 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.1 | 9 | 0 | 4,8,14,15,16,17,18,19,20 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.1 | 8 | 0 | 20,19,18,13,10,7,6,4 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.2 | 8 | 0 | 5,10,13,14,17,18,19,20 |
| CserDQD_SDH2_Rev | CserDQD/SDH1.2 | 8 | 0 | 20,19,18,13,10,7,6,4 |
| CserDQD_SDH2_Rev | CserDQD/SDH2 | 8 | 0 | 6,7,9,12,14,16,18,20 |
| CserDQD_SDH2_Rev | CserDQD/SDH2 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH2_Rev | CserDQD/SDH3 | 9 | 1 | 6,7,10,12,13,16,17,18,19 |
| CserDQD_SDH2_Rev | CserDQD/SDH3 | 4 | 0 | 20,19,10,1 |
| CserDQD_SDH2_Rev | CserDQD/SDH4 | 6 | 6 | 3,4,5,8,9,14 |
| CserDQD_SDH2_Rev | CserDQD/SDH4 | 7 | 1 | 19,15,14,13,11,2,1 |
| CserDQD_SDH2_Rev | CserDQD/SDH5 | 9 | 0 | 1,6,7,8,10,11,12,15,20 |
| CserDQD_SDH2_Rev | CserDQD/SDH5 | 7 | 1 | 19,13,10,7,6,5,4 |
| CserDQD_SDH3_Fwd | CserDQD/SDH1.1 | 7 | 0 | 7,10,11,14,16,19,20 |
| CserDQD_SDH3_Fwd | CserDQD/SDH1.1 | 7 | 6 | 14,10,9,8,6,3,2 |
| CserDQD_SDH3_Fwd | CserDQD/SDH1.2 | 7 | 0 | 7,10,11,14,16,19,20 |
| CserDQD_SDH3_Fwd | CserDQD/SDH1.2 | 7 | 6 | 14,10,9,8,6,3,2 |
| CserDQD_SDH3_Fwd | CserDQD/SDH2 | 3 | 7 | 7,10,13 |
| CserDQD_SDH3_Fwd | CserDQD/SDH2 | 8 | 0 | 20,18,17,15,14,10,8,7 |
| CserDQD_SDH3_Fwd | CserDQD/SDH3 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH3_Fwd | CserDQD/SDH3 | 8 | 2 | 18,14,13,11,9,8,6,1 |
| CserDQD_SDH3_Fwd | CserDQD/SDH4 | 3 | 7 | 1,7,13 |
| CserDQD_SDH3_Fwd | CserDQD/SDH4 | 8 | 2 | 18,16,14,13,12,8,5,3 |
| CserDQD_SDH3_Fwd | CserDQD/SDH5 | 5 | 0 | 8,9,18,19,20 |
| CserDQD_SDH3_Fwd | CserDQD/SDH5 | 8 | 0 | 20,19,14,13,12,7,3,1 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.1 | 8 | 3 | 2,3,4,5,6,9,14,17 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.1 | 7 | 2 | 18,17,16,15,12,9,3 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.2 | 8 | 3 | 2,3,4,5,6,9,14,17 |
| CserDQD_SDH3_Rev | CserDQD/SDH1.2 | 7 | 2 | 18,17,16,15,12,9,3 |
| CserDQD_SDH3_Rev | CserDQD/SDH2 | 7 | 0 | 2,6,10,12,17,18,20 |
| CserDQD_SDH3_Rev | CserDQD/SDH2 | 4 | 5 | 15,6,5,3 |
| CserDQD_SDH3_Rev | CserDQD/SDH3 | 8 | 0 | 1,2,3,12,13,16,17,20 |
| CserDQD_SDH3_Rev | CserDQD/SDH3 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH3_Rev | CserDQD/SDH4 | 7 | 0 | 4,5,6,8,12,19,20 |
| CserDQD_SDH3_Rev | CserDQD/SDH4 | 3 | 8 | 12,8,6 |
| CserDQD_SDH3_Rev | CserDQD/SDH5 | 8 | 0 | 3,5,8,11,12,13,15,20 |
| CserDQD_SDH3_Rev | CserDQD/SDH5 | 7 | 0 | 20,19,15,10,9,6,2 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.1 | 7 | 1 | 1,2,11,12,13,15,19 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.1 | 7 | 0 | 20,19,17,15,13,11,4 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.2 | 7 | 1 | 1,2,11,12,13,15,19 |
| CserDQD_SDH4_Fwd | CserDQD/SDH1.2 | 7 | 0 | 20,19,17,15,13,11,4 |
| CserDQD_SDH4_Fwd | CserDQD/SDH2 | 4 | 6 | 1,7,9,14 |
| CserDQD_SDH4_Fwd | CserDQD/SDH2 | 8 | 0 | 20,15,14,13,12,8,6,5 |
| CserDQD_SDH4_Fwd | CserDQD/SDH3 | 6 | 1 | 1,4,7,9,14,19 |
| CserDQD_SDH4_Fwd | CserDQD/SDH3 | 7 | 0 | 20,19,15,12,4,3,1 |
| CserDQD_SDH4_Fwd | CserDQD/SDH4 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH4_Fwd | CserDQD/SDH4 | 8 | 0 | 20,14,12,10,9,4,3,2 |
| CserDQD_SDH4_Fwd | CserDQD/SDH5 | 7 | 1 | 1,7,10,14,15,18,19 |
| CserDQD_SDH4_Fwd | CserDQD/SDH5 | 7 | 0 | 20,16,14,13,11,7,5 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.1 | 9 | 0 | 3,4,7,9,10,16,17,18,20 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.1 | 8 | 0 | 20,19,15,14,12,10,7,6 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.2 | 7 | 1 | 7,10,12,13,16,18,19 |
| CserDQD_SDH4_Rev | CserDQD/SDH1.2 | 8 | 0 | 20,19,15,14,12,10,7,6 |
| CserDQD_SDH4_Rev | CserDQD/SDH2 | 7 | 1 | 2,6,7,10,15,18,19 |
| CserDQD_SDH4_Rev | CserDQD/SDH2 | 4 | 2 | 18,14,13,12 |
| CserDQD_SDH4_Rev | CserDQD/SDH3 | 8 | 0 | 5,6,8,11,13,15,19,20 |
| CserDQD_SDH4_Rev | CserDQD/SDH3 | 5 | 2 | 18,14,13,12,3 |
| CserDQD_SDH4_Rev | CserDQD/SDH4 | 8 | 0 | 2,6,7,10,15,18,19,20 |
| CserDQD_SDH4_Rev | CserDQD/SDH4 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH4_Rev | CserDQD/SDH5 | 7 | 5 | 2,6,7,8,11,12,15 |
| CserDQD_SDH4_Rev | CserDQD/SDH5 | 8 | 2 | 18,16,14,13,12,9,6,3 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.1 | 3 | 1 | 1,8,19 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.1 | 9 | 2 | 18,16,15,13,10,7,5,2,1 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.2 | 3 | 1 | 1,8,19 |
| CserDQD_SDH5_Fwd | CserDQD/SDH1.2 | 9 | 2 | 18,16,15,13,10,7,5,2,1 |
| CserDQD_SDH5_Fwd | CserDQD/SDH2 | 5 | 1 | 2,7,8,13,19 |
| CserDQD_SDH5_Fwd | CserDQD/SDH2 | 9 | 3 | 17,15,13,12,11,8,4,2,1 |
| CserDQD_SDH5_Fwd | CserDQD/SDH3 | 8 | 2 | 1,2,3,6,7,10,16,18 |
| CserDQD_SDH5_Fwd | CserDQD/SDH3 | 7 | 2 | 18,15,14,13,6,5,1 |
| CserDQD_SDH5_Fwd | CserDQD/SDH4 | 7 | 1 | 2,7,8,11,12,13,19 |
| CserDQD_SDH5_Fwd | CserDQD/SDH4 | 9 | 1 | 19,16,14,13,12,11,6,5,2 |
| CserDQD_SDH5_Fwd | CserDQD/SDH5 | 0 | 20 | - **(designed site)** |
| CserDQD_SDH5_Fwd | CserDQD/SDH5 | 9 | 0 | 20,17,15,13,12,11,9,7,2 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.1 | 9 | 0 | 4,7,8,13,14,15,18,19,20 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.1 | 7 | 1 | 19,16,9,8,6,5,1 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.2 | 9 | 0 | 4,7,8,13,14,15,18,19,20 |
| CserDQD_SDH5_Rev | CserDQD/SDH1.2 | 7 | 1 | 19,16,9,8,6,5,1 |
| CserDQD_SDH5_Rev | CserDQD/SDH2 | 8 | 1 | 1,5,7,8,11,13,15,19 |
| CserDQD_SDH5_Rev | CserDQD/SDH2 | 8 | 0 | 20,19,18,14,10,9,7,1 |
| CserDQD_SDH5_Rev | CserDQD/SDH3 | 9 | 0 | 1,3,7,9,13,17,18,19,20 |
| CserDQD_SDH5_Rev | CserDQD/SDH3 | 8 | 0 | 20,13,12,10,9,4,3,1 |
| CserDQD_SDH5_Rev | CserDQD/SDH4 | 7 | 9 | 3,4,5,6,8,10,11 |
| CserDQD_SDH5_Rev | CserDQD/SDH4 | 8 | 0 | 20,14,13,11,8,5,4,2 |
| CserDQD_SDH5_Rev | CserDQD/SDH5 | 8 | 1 | 1,5,8,13,14,15,18,19 |
| CserDQD_SDH5_Rev | CserDQD/SDH5 | 0 | 20 | - **(designed site)** |

## Candidate amplicons under four rule sets

### strict  (<=3 mm, perfect 5-nt 3' end, 80-600 bp)
- products: 6  |  cross-products (pair on a non-target transcript): 0
  - CserDQD_SDH1.1 → **CserDQD/SDH1.1** 120 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH1.2 → **CserDQD/SDH1.2** 180 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH2 → **CserDQD/SDH2** 240 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH3 → **CserDQD/SDH3** 300 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH4 → **CserDQD/SDH4** 360 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH5 → **CserDQD/SDH5** 420 bp (mm 0/0, 3'-run 20/20) [TARGET]

### medium  (<=4 mm, perfect 8-nt 3' end, 60-1000 bp)
- products: 6  |  cross-products (pair on a non-target transcript): 0
  - CserDQD_SDH1.1 → **CserDQD/SDH1.1** 120 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH1.2 → **CserDQD/SDH1.2** 180 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH2 → **CserDQD/SDH2** 240 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH3 → **CserDQD/SDH3** 300 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH4 → **CserDQD/SDH4** 360 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH5 → **CserDQD/SDH5** 420 bp (mm 0/0, 3'-run 20/20) [TARGET]

### loose   (<=6 mm, perfect 8-nt 3' end, 60-1000 bp)
- products: 6  |  cross-products (pair on a non-target transcript): 0
  - CserDQD_SDH1.1 → **CserDQD/SDH1.1** 120 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH1.2 → **CserDQD/SDH1.2** 180 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH2 → **CserDQD/SDH2** 240 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH3 → **CserDQD/SDH3** 300 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH4 → **CserDQD/SDH4** 360 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH5 → **CserDQD/SDH5** 420 bp (mm 0/0, 3'-run 20/20) [TARGET]

### permissive (<=6 mm each, <=12 mm total, perfect 5-nt 3' end, 60-1000 bp)
- products: 10  |  cross-products (pair on a non-target transcript): 4
  - CserDQD_SDH1.1 → **CserDQD/SDH1.1** 120 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH1.2 → **CserDQD/SDH1.2** 180 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH2 → **CserDQD/SDH2** 240 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH3 → **CserDQD/SDH2** 636 bp (mm 6/4, 3'-run 6/5) [CROSS-PRODUCT]
  - CserDQD_SDH3 → **CserDQD/SDH2** 309 bp (mm 3/4, 3'-run 7/5) [CROSS-PRODUCT]
  - CserDQD_SDH3 → **CserDQD/SDH3** 300 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH3 → **CserDQD/SDH4** 630 bp (mm 6/3, 3'-run 6/8) [CROSS-PRODUCT]
  - CserDQD_SDH3 → **CserDQD/SDH4** 300 bp (mm 3/3, 3'-run 7/8) [CROSS-PRODUCT]
  - CserDQD_SDH4 → **CserDQD/SDH4** 360 bp (mm 0/0, 3'-run 20/20) [TARGET]
  - CserDQD_SDH5 → **CserDQD/SDH5** 420 bp (mm 0/0, 3'-run 20/20) [TARGET]

