# Backup CserDQD/SDH3 pairs — zero competent sites on SDH2/SDH4

Mandatory filters: 20-nt, exact-unique to SDH3, GC 40–60 %, no homopolymer >4 nt, acceptable 3′ tail, no ≥7-nt self-complementary hairpin, product 300 bp, ΔTm ≤ 2.0 °C, no ≥8-nt 3′-terminal complementarity with the other ten primers, and zero competent foreign sites (≤6 mismatches with ≥5-nt perfect 3′ run) at primer level (Tier A) or pair level (Tier B).

## Recommended backup (best-ranked)

- **tier**: A
- **fwd_seq**: TCATTAGCAATGGCCTCAGG
- **rev_seq**: ATCACCTTCATACTGACCGC
- **fwd_pos**: 232-251
- **rev_pos**: 512-531
- **fwd_region**: 5UTR
- **rev_region**: CDS
- **product**: 300
- **fwd_gc**: 50.0
- **rev_gc**: 50.0
- **fwd_tm**: 53.4
- **rev_tm**: 53.6
- **delta_tm**: 0.2
- **worst_foreign_fwd**: CserDQD/SDH2 7mm tail2
- **worst_foreign_rev**: CserDQD/SDH2 2mm tail1
- **competent_foreign_fwd**: 0
- **competent_foreign_rev**: 0
- **very_permissive_foreign**: 2
- **alt_products_on_SDH3**: none
- **max_3prime_dimer_vs_other_pairs**: 2
- **score**: 7.66

## In-service SDH3 pair, for comparison

- forward `TGGAGCCGACATTGTGAAGA` 55.2 °C, reverse `TGACCGACAGGCTTTCCAAT` 55.1 °C, product 300 bp
- competent foreign sites: 4 forward + 2 reverse (they create the SDH4 300 bp and SDH2 309 bp cross-pairings)

## Genome-wide verification of the recommended backup pair

Exhaustive seed-free scan of all 48.1 Mb in 29124 transcripts; every site with ≤4 mismatches is listed.

### `TCATTAGCAATGGCCTCAGG` — 64 site(s) with ≤4 mismatches
- NODE_20612_length_2133_cov_1535.833013_g8101_i1  strand +  start 23274132  0 mm  3′-run 20 nt
- NODE_45672_length_1214_cov_496.389700_g21067_i0  strand -  start 39494611  3 mm  3′-run 10 nt
- NODE_1279_length_5075_cov_28.709511_g374_i5  strand +  start 2443457  3 mm  3′-run 7 nt
- NODE_13263_length_2592_cov_362.226111_g4385_i1  strand -  start 16582851  3 mm  3′-run 4 nt
- NODE_19596_length_2184_cov_264.091803_g8206_i1  strand -  start 22431166  3 mm  3′-run 4 nt
- NODE_14085_length_2528_cov_124.572408_g2739_i1  strand -  start 17412489  3 mm  3′-run 3 nt
- NODE_26211_length_1876_cov_109.020252_g12111_i0  strand +  start 27823787  3 mm  3′-run 0 nt
- NODE_54324_length_970_cov_162.310532_g25930_i0  strand -  start 42748660  4 mm  3′-run 10 nt
- NODE_156_length_7927_cov_111.965220_g67_i0  strand -  start 441926  4 mm  3′-run 9 nt
- NODE_43579_length_1277_cov_339.790717_g20005_i0  strand -  start 38535880  4 mm  3′-run 8 nt
- NODE_225_length_7354_cov_109.220945_g103_i0  strand -  start 613115  4 mm  3′-run 7 nt
- NODE_5930_length_3420_cov_57.202611_g2954_i0  strand -  start 8694381  4 mm  3′-run 7 nt

### `ATCACCTTCATACTGACCGC` — 49 site(s) with ≤4 mismatches
- NODE_20612_length_2133_cov_1535.833013_g8101_i1  strand -  start 23274412  0 mm  3′-run 20 nt
- NODE_18517_length_2241_cov_255.196168_g8755_i0  strand -  start 21470179  2 mm  3′-run 1 nt
- NODE_43298_length_1286_cov_389.512530_g19870_i0  strand +  start 38408306  3 mm  3′-run 1 nt
- NODE_21648_length_2083_cov_1398.299410_g8244_i1  strand -  start 24141025  3 mm  3′-run 1 nt
- NODE_29040_length_1767_cov_22.706636_g12409_i1  strand -  start 29909132  3 mm  3′-run 1 nt
- NODE_6525_length_3323_cov_145.985339_g3252_i0  strand +  start 9420262  4 mm  3′-run 11 nt
- NODE_16381_length_2369_cov_152.101724_g4750_i4  strand -  start 19543062  4 mm  3′-run 9 nt
- NODE_48423_length_1133_cov_173.251845_g22500_i0  strand +  start 40607449  4 mm  3′-run 7 nt
- NODE_52198_length_1026_cov_9.942682_g24624_i0  strand +  start 42015108  4 mm  3′-run 7 nt
- NODE_15759_length_2411_cov_40.693480_g7528_i0  strand -  start 18967877  4 mm  3′-run 7 nt
- NODE_11266_length_2755_cov_65.266814_g4193_i1  strand +  start 14477831  4 mm  3′-run 6 nt
- NODE_14385_length_2505_cov_54.660423_g6058_i1  strand +  start 17699931  4 mm  3′-run 6 nt

## All passing backups, ranked

| tier | fwd_seq | rev_seq | fwd_pos | rev_pos | fwd_region | rev_region | fwd_tm | rev_tm | delta_tm | fwd_gc | rev_gc | worst_foreign_fwd | worst_foreign_rev | very_permissive_foreign | alt_products_on_SDH3 | max_3prime_dimer_vs_other_pairs | score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | TCATTAGCAATGGCCTCAGG | ATCACCTTCATACTGACCGC | 232-251 | 512-531 | 5UTR | CDS | 53.4 | 53.6 | 0.2 | 50.0 | 50.0 | CserDQD/SDH2 7mm tail2 | CserDQD/SDH2 2mm tail1 | 2 | none | 2 | 7.66 |
| A | CACCCAGAAGAAGGGATGGT | GCAGACTCTCTCTTATGCGT | 1552-1571 | 1832-1851 | CDS | CDS | 54.7 | 53.4 | 1.3 | 55.0 | 50.0 | CserDQD/SDH4 3mm tail1 | CserDQD/SDH1.1 7mm tail1 | 0 | none | 4 | 10.49 |
| A | AAGTCCCAACCTACCTACCC | AACTGTGGTGCCCATCAATG | 49-68 | 329-348 | 5UTR | CDS | 53.7 | 54.5 | 0.8 | 55.0 | 50.0 | CserDQD/SDH4 7mm tail0 | CserDQD/SDH1.1 7mm tail2 | 3 | none | 3 | 17.53 |
| A | TTGGTGAACCCTCCTCTGAG | CAATCTCCTTCGCCTTTCCC | 89-108 | 369-388 | 5UTR | CDS | 54.4 | 54.0 | 0.4 | 55.0 | 55.0 | CserDQD/SDH4 8mm tail2 | CserDQD/SDH4 5mm tail1 | 2 | none | 2 | 18.1 |
| A | TGCCAATCGTACATTCGAGC | TCAAAGCCCTGTCGGATCAA | 1467-1486 | 1747-1766 | CDS | CDS | 54.8 | 54.9 | 0.1 | 50.0 | 50.0 | CserDQD/SDH2 5mm tail1 | CserDQD/SDH1.1 5mm tail2 | 6 | none | 2 | 19.96 |
| A | TGTTGATCCGACAGGGCTTT | TTCATGTCCTGGCTCTTCCT | 1745-1764 | 2025-2044 | CDS | 3UTR | 55.1 | 54.1 | 1.1 | 50.0 | 50.0 | CserDQD/SDH1.1 4mm tail2 | CserDQD/SDH2 6mm tail2 | 5 | none | 2 | 20.6 |
| A | GCCAATCGTACATTCGAGCG | TTCAAAGCCCTGTCGGATCA | 1468-1487 | 1748-1767 | CDS | CDS | 56.0 | 54.9 | 1.1 | 55.0 | 50.0 | CserDQD/SDH2 5mm tail2 | CserDQD/SDH1.1 6mm tail1 | 3 | none | 2 | 21.36 |
| A | CCAATCGTACATTCGAGCGA | GTTCAAAGCCCTGTCGGATC | 1469-1488 | 1749-1768 | CDS | CDS | 54.3 | 54.6 | 0.3 | 50.0 | 55.0 | CserDQD/SDH2 5mm tail3 | CserDQD/SDH1.1 6mm tail0 | 5 | none | 2 | 22.0 |
| A | TTTCCACCCAGAAGAAGGGA | ACTCTCTCTTATGCGTGTCT | 1548-1567 | 1828-1847 | CDS | CDS | 53.2 | 51.9 | 1.3 | 50.0 | 45.0 | CserDQD/SDH4 2mm tail4 | CserDQD/SDH4 6mm tail1 | 4 | none | 3 | 25.28 |
| A | AGTCCCAACCTACCTACCCT | CAACTGTGGTGCCCATCAAT | 50-69 | 330-349 | 5UTR | CDS | 54.8 | 54.5 | 0.3 | 55.0 | 50.0 | CserDQD/SDH1.1 8mm tail3 | CserDQD/SDH1.1 6mm tail1 | 6 | none | 2 | 25.47 |
| A | CATTAGCAATGGCCTCAGGA | CATCACCTTCATACTGACCG | 233-252 | 513-532 | 5UTR | CDS | 53.4 | 51.9 | 1.5 | 50.0 | 50.0 | CserDQD/SDH2 6mm tail3 | CserDQD/SDH2 2mm tail0 | 6 | none | 2 | 26.71 |
| A | ATGTTGATCCGACAGGGCTT | TCATGTCCTGGCTCTTCCTT | 1744-1763 | 2024-2043 | CDS | 3UTR | 55.2 | 54.1 | 1.2 | 50.0 | 50.0 | CserDQD/SDH1.1 4mm tail1 | CserDQD/SDH4 6mm tail1 | 7 | none | 3 | 26.89 |
| A | CGTTTGCTACTTCGGATGCT | GCTTAAGTGCATCTTGGCGC | 260-279 | 540-559 | CDS | CDS | 54.4 | 56.3 | 1.9 | 50.0 | 55.0 | CserDQD/SDH1.1 6mm tail2 | CserDQD/SDH4 7mm tail0 | 4 | none | 4 | 27.15 |
| A | CTTATGGTGCACTTGAGGCG | TAAGTGCAGCCTCTTTGTGG | 902-921 | 1182-1201 | CDS | CDS | 55.3 | 53.6 | 1.6 | 55.0 | 50.0 | CserDQD/SDH4 3mm tail0 | CserDQD/SDH2 2mm tail0 | 5 | none | 1 | 27.38 |
| A | TCGTTTGCTACTTCGGATGC | CTTAAGTGCATCTTGGCGCT | 259-278 | 539-558 | CDS | CDS | 54.4 | 54.6 | 0.2 | 50.0 | 50.0 | CserDQD/SDH4 6mm tail0 | CserDQD/SDH5 7mm tail2 | 9 | none | 3 | 28.63 |
| A | GTTGGTGAACCCTCCTCTGA | AATCTCCTTCGCCTTTCCCA | 88-107 | 368-387 | 5UTR | CDS | 54.7 | 54.2 | 0.5 | 55.0 | 50.0 | CserDQD/SDH1.2 7mm tail0 | CserDQD/SDH4 5mm tail2 | 7 | none | 2 | 28.88 |
| A | ACCCAGAAGAAGGGATGGTC | CGCAGACTCTCTCTTATGCG | 1553-1572 | 1833-1852 | CDS | CDS | 54.4 | 54.6 | 0.2 | 55.0 | 55.0 | CserDQD/SDH4 4mm tail0 | CserDQD/SDH1.1 7mm tail8 | 8 | none | 2 | 35.93 |
| A | ACTGGCTGACAAAGTTGACG | GCTGGCAAACCAGTGAAGTT | 1497-1516 | 1777-1796 | CDS | CDS | 54.4 | 54.5 | 0.1 | 50.0 | 50.0 | CserDQD/SDH1.1 5mm tail4 | CserDQD/SDH2 4mm tail1 | 12 | none | 3 | 37.48 |
| A | TTGTCTCAGGCAGTTCAGTC | AGATTGCCAATAGCCTCAGC | 417-436 | 697-716 | CDS | CDS | 53.3 | 53.8 | 0.5 | 50.0 | 50.0 | CserDQD/SDH4 7mm tail7 | CserDQD/SDH2 4mm tail1 | 12 | none | 2 | 38.96 |
| A | GTTGATCCGACAGGGCTTTG | TTTCATGTCCTGGCTCTTCC | 1746-1765 | 2026-2045 | CDS | 3UTR | 54.9 | 53.0 | 1.9 | 55.0 | 50.0 | CserDQD/SDH1.1 4mm tail3 | CserDQD/SDH4 6mm tail2 | 9 | none | 1 | 39.64 |
| A | TTGAAGATGGACTACGAGGC | GGGTTGCATCAATGTTGGGT | 1322-1341 | 1602-1621 | CDS | CDS | 53.0 | 54.5 | 1.6 | 50.0 | 50.0 | CserDQD/SDH2 2mm tail0 | CserDQD/SDH4 7mm tail3 | 11 | none | 3 | 39.79 |
| A | CTCAAACGAGTAGTAGTGGA | CTCCCAATTCCATAGCTAGC | 278-297 | 558-577 | CDS | CDS | 49.9 | 51.2 | 1.4 | 45.0 | 50.0 | CserDQD/SDH5 6mm tail0 | CserDQD/SDH2 4mm tail0 | 8 | none | 2 | 41.49 |
| A | TCCCAAGTCCCAACCTACCT | GTGGTGCCCATCAATGGAGT | 45-64 | 325-344 | 5UTR | CDS | 55.4 | 56.4 | 0.9 | 55.0 | 55.0 | CserDQD/SDH2 7mm tail4 | CserDQD/SDH2 7mm tail1 | 8 | none | 4 | 41.61 |
| A | AGATCAGGCTGGATTGTCTC | CCTCAGCAGATGGAGTATTG | 404-423 | 684-703 | CDS | CDS | 52.8 | 51.5 | 1.3 | 50.0 | 50.0 | CserDQD/SDH5 3mm tail0 | CserDQD/SDH2 2mm tail0 | 11 | none | 2 | 42.04 |
| A | ATGATTTCCACCCAGAAGAA | TCTCTTATGCGTGTCTTGAC | 1544-1563 | 1824-1843 | CDS | CDS | 49.3 | 51.2 | 1.9 | 40.0 | 45.0 | CserDQD/SDH4 3mm tail0 | CserDQD/SDH1.1 7mm tail3 | 5 | none | 2 | 45.16 |
| A | AAGCCAGATACCAAAGTCTA | GTAAGCAGTTAGCGTTAAGT | 991-1010 | 1271-1290 | CDS | CDS | 48.8 | 48.8 | 0.0 | 40.0 | 40.0 | CserDQD/SDH4 5mm tail1 | CserDQD/SDH5 5mm tail1 | 5 | none | 4 | 45.64 |
| A | ATTTCCACCCAGAAGAAGGG | CTCTCTCTTATGCGTGTCTT | 1547-1566 | 1827-1846 | CDS | CDS | 52.2 | 50.6 | 1.6 | 50.0 | 45.0 | CserDQD/SDH4 2mm tail3 | CserDQD/SDH4 5mm tail2 | 10 | none | 4 | 46.75 |
| A | CATCTGTAGGAATGAAACCC | GAATTCAAAGCAACTTGGGC | 1586-1605 | 1866-1885 | CDS | 3UTR | 48.9 | 50.9 | 1.9 | 45.0 | 45.0 | CserDQD/SDH4 1mm tail0 | CserDQD/SDH5 7mm tail3 | 7 | none | 2 | 47.06 |
| A | TGGTGAACCCTCCTCTGAGA | CCAATCTCCTTCGCCTTTCC | 90-109 | 370-389 | 5UTR | CDS | 55.5 | 54.0 | 1.5 | 55.0 | 55.0 | CserDQD/SDH2 7mm tail1 | CserDQD/SDH4 5mm tail0 | 10 | none | 2 | 47.55 |
| A | ATGGTCCTTGCAAACACTAC | CTGTGATGGATAAACGCAGA | 1567-1586 | 1847-1866 | CDS | 3UTR | 51.5 | 51.5 | 0.0 | 45.0 | 45.0 | CserDQD/SDH2 3mm tail2 | CserDQD/SDH4 5mm tail1 | 11 | none | 2 | 48.18 |
