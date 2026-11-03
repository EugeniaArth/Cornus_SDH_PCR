# Exhaustive transcriptome scan (no BLAST seeding)

Scanned every one of 48.1 Mb of the assembly at 47.6 M valid 20-nt windows, for all 12 primers on both strands. Reported every site with ≤6 mismatches.

## Sites per primer (mismatch distribution)

| Primer | strand | 0 mm | ≤1 | ≤2 | ≤3 | ≤4 | ≤5 | ≤6 | best 3'-tail at ≤3 mm |
|---|---|---|---|---|---|---|---|---|---|
| CserDQD_SDH1.1_Fwd | plus | 0 | 0 | 1 | 2 | 34 | 325 | 2239 | 0 |
| CserDQD_SDH1.1_Fwd | minus | 0 | 0 | 0 | 1 | 14 | 172 | 1377 | 0 |
| CserDQD_SDH1.1_Rev | plus | 0 | 0 | 0 | 0 | 17 | 222 | 1530 | 0 |
| CserDQD_SDH1.1_Rev | minus | 1 | 1 | 1 | 5 | 31 | 273 | 2110 | 0 |
| CserDQD_SDH1.2_Fwd | plus | 1 | 1 | 1 | 6 | 45 | 368 | 2412 | 0 |
| CserDQD_SDH1.2_Fwd | minus | 0 | 0 | 0 | 1 | 37 | 296 | 2029 | 0 |
| CserDQD_SDH1.2_Rev | plus | 0 | 0 | 0 | 0 | 17 | 170 | 1297 | 0 |
| CserDQD_SDH1.2_Rev | minus | 1 | 1 | 1 | 4 | 34 | 284 | 1873 | 0 |
| CserDQD_SDH2_Fwd | plus | 1 | 1 | 1 | 3 | 21 | 179 | 1426 | 0 |
| CserDQD_SDH2_Fwd | minus | 0 | 0 | 0 | 1 | 21 | 174 | 1402 | 0 |
| CserDQD_SDH2_Rev | plus | 0 | 0 | 0 | 2 | 17 | 145 | 1106 | 0 |
| CserDQD_SDH2_Rev | minus | 1 | 1 | 1 | 3 | 22 | 255 | 1760 | 0 |
| CserDQD_SDH3_Fwd | plus | 1 | 1 | 1 | 10 | 50 | 392 | 2505 | 0 |
| CserDQD_SDH3_Fwd | minus | 0 | 0 | 0 | 2 | 24 | 229 | 1754 | 0 |
| CserDQD_SDH3_Rev | plus | 0 | 0 | 0 | 1 | 8 | 144 | 1145 | 0 |
| CserDQD_SDH3_Rev | minus | 1 | 1 | 1 | 2 | 27 | 225 | 1573 | 0 |
| CserDQD_SDH4_Fwd | plus | 1 | 1 | 1 | 8 | 54 | 479 | 3232 | 0 |
| CserDQD_SDH4_Fwd | minus | 0 | 0 | 0 | 4 | 34 | 332 | 2449 | 0 |
| CserDQD_SDH4_Rev | plus | 0 | 0 | 0 | 2 | 19 | 196 | 1538 | 0 |
| CserDQD_SDH4_Rev | minus | 1 | 2 | 2 | 2 | 20 | 210 | 1457 | 0 |
| CserDQD_SDH5_Fwd | plus | 1 | 1 | 1 | 3 | 16 | 179 | 1468 | 0 |
| CserDQD_SDH5_Fwd | minus | 0 | 0 | 0 | 0 | 19 | 146 | 1201 | 0 |
| CserDQD_SDH5_Rev | plus | 0 | 0 | 0 | 1 | 19 | 129 | 1189 | 0 |
| CserDQD_SDH5_Rev | minus | 1 | 1 | 1 | 1 | 18 | 153 | 1365 | 0 |

## Candidate products (exhaustive, 60-1500 bp)

- all pairs (≤8 total mismatches): **10**
- Primer-BLAST-strict (≤3 mm, perfect 5-nt 3' ends): **5**
- extension-competent relaxed set (≤5 mm, perfect 8-nt 3' ends): **5**

| Pair | Contig | Size | mm F/R | 3'-tail F/R | Class | Annotation |
|---|---|---|---|---|---|---|
| CserDQD_SDH1.1 | NODE_10980_length_2779_cov_121.715385_g1243_i5 | 1012 | 2/0 | 0/20 | TARGET_CONTIG | - |
| CserDQD_SDH1.2 | NODE_10980_length_2779_cov_121.715385_g1243_i5 | 180 | 0/0 | 20/20 | TARGET_CONTIG | - |
| CserDQD_SDH2 | NODE_18517_length_2241_cov_255.196168_g8755_i0 | 240 | 0/0 | 20/20 | TARGET_CONTIG | - |
| CserDQD_SDH3 | NODE_20612_length_2133_cov_1535.833013_g8101_i1 | 300 | 0/0 | 20/20 | TARGET_CONTIG | - |
| CserDQD_SDH3 | NODE_20612_length_2133_cov_1535.833013_g8101_i1 | 651 | 3/0 | 3/20 | TARGET_CONTIG | - |
| CserDQD_SDH4 | NODE_21648_length_2083_cov_1398.299410_g8244_i1 | 360 | 0/0 | 20/20 | TARGET_CONTIG | - |
| CserDQD_SDH5 | NODE_29040_length_1767_cov_22.706636_g12409_i1 | 420 | 0/0 | 20/20 | TARGET_CONTIG | - |
| CserDQD_SDH3 | NODE_21648_length_2083_cov_1398.299410_g8244_i1 | 300 | 3/3 | 7/8 | UNRELATED_CONTIG | - |
| CserDQD_SDH3 | NODE_18517_length_2241_cov_255.196168_g8755_i0 | 309 | 3/4 | 7/5 | UNRELATED_CONTIG | - |
| CserDQD_SDH4 | NODE_18517_length_2241_cov_255.196168_g8755_i0 | 375 | 4/4 | 6/2 | UNRELATED_CONTIG | - |

## Reading

- `mm` counts mismatches over the full 20 nt, so a low-identity partial match is not a priming site; the `3'-tail` column is what actually permits extension.

- Every site is reported, including those BLAST could not seed.

