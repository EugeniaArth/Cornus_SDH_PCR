# Backup CserDQD/SDH3 pair — verification against the in-service set

The same exhaustive machinery is run on four configurations: the six pairs in service and three SDH3 replacements. All 36 forward×reverse combinations of each 6-plex are examined, on the six curated transcripts and on all 48.1 Mb of the assembly (seed-free, every 20-nt window, both strands).

Note: only **five** of the six genes have their own contig in the clustered assembly (SDH1.1 is collapsed into NODE_10980 with SDH1.2), so the maximum number of detectable designed amplicons transcriptome-wide is 5, not 6.


## in-service (baseline)

- SDH3 primers in service: F `TGGAGCCGACATTGTGAAGA` (738-757) / R `TGACCGACAGGCTTTCCAAT` (1018-1037), 300 bp

### Stage 1 — six curated transcripts, all 36 forward×reverse combinations

| rule | designed amplicons | other products | foreign-paralog / cross-pair detail |
|---|---|---|---|
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 6 | 2 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH3+R:SDH4 on SDH4 948 bp (F 3 mm/3′7, R 0 mm/3′20) *[cross]* |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 6 | 6 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH2+R:SDH3 on SDH2 494 bp (F 0 mm/3′20, R 4 mm/3′5) *[cross]*; F:SDH2+R:SDH3 on SDH3 488 bp (F 5 mm/3′5, R 0 mm/3′20) *[cross]*; F:SDH3+R:SDH2 on SDH2 382 bp (F 6 mm/3′6, R 0 mm/3′20) *[cross]* (+2 more) |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 6 | 3 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8); F:SDH2+R:SDH3 on SDH4 488 bp (F 5 mm/3′9, R 3 mm/3′8) *[cross]* |
| permissive      (each ≤6 mm, 3′ runs ≥5 nt) | 6 | 16 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8); F:SDH2+R:SDH1.2 on SDH3 571 bp (F 5 mm/3′5, R 7 mm/3′5) *[cross]*; F:SDH2+R:SDH3 on SDH2 308 bp (F 0 mm/3′20, R 8 mm/3′6) *[cross]* (+12 more) |
| very permissive (each ≤8 mm, 3′ runs ≥3 nt) | 6 | 64 | F:SDH1.1+R:SDH1.1 on SDH1.1 504 bp (F 8 mm/3′3, R 0 mm/3′20); F:SDH1.1+R:SDH1.1 on SDH1.1 756 bp (F 8 mm/3′4, R 0 mm/3′20); F:SDH1.1+R:SDH1.2 on SDH1.1 428 bp (F 8 mm/3′3, R 8 mm/3′8) *[cross]*; F:SDH1.1+R:SDH1.2 on SDH1.1 60 bp (F 8 mm/3′4, R 8 mm/3′3) *[cross]* (+60 more) |

### Stage 2 — whole transcriptome, every window, same rules as TX_EXHAUSTIVE_SCAN.md

| primer | sites ≤2 mm | ≤4 mm | ≤6 mm | longest 3′ run at ≤4 mm |
|---|---|---|---|---|
| `AAGAGTCTGGTCCATTTCGG` | 2 | 39 | 2995 | 20 |
| `ACCTCAGTCACAGGAATTCG` | 1 | 37 | 2554 | 20 |
| `ACTCAGCAACCCAATAGTGC` | 1 | 51 | 3170 | 20 |
| `ACTTGCTCAGTCCATAGGTG` | 1 | 35 | 2669 | 20 |
| `AGAGGAAGCACTAGGAGGTT` | 1 | 88 | 5681 | 20 |
| `ATAAGCACCTTACTTGAAGC` | 1 | 48 | 3640 | 20 |
| `GCAAGCAGATTTCAGAGCTG` | 1 | 82 | 4441 | 20 |
| `GCACTTCGATTAGCCATGGA` | 1 | 42 | 2828 | 20 |
| `GTTCACTGGTGGATTAAGTG` | 1 | 48 | 3616 | 16 |
| `TGACCGACAGGCTTTCCAAT` ← **backup** | 1 | 35 | 2718 | 20 |
| `TGCCACATCAGTGATATCCG` | 1 | 39 | 2866 | 20 |
| `TGGAGCCGACATTGTGAAGA` ← **backup** | 1 | 74 | 4259 | 20 |

| rule | designed | cross-pair exact | on-family (paralog) | off-family (unrelated contigs) | worst example |
|---|---|---|---|---|---|
| designed        (both primers 0 mm) | 5 | 1 | 0 | 0 | — |
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 5 | 1 | 1 | 0 | F:SDH3+R:SDH4→NODE_21648_length_20… 948 bp, 3 mm, 3′ 7/20 |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 5 | 1 | 5 | 0 | F:SDH3+R:SDH4→NODE_21648_length_20… 948 bp, 3 mm, 3′ 7/20 |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 5 | 1 | 1 | 0 | F:SDH2+R:SDH3→NODE_21648_length_20… 488 bp, 8 mm, 3′ 9/8 |
| any             (total ≤8 mm, no 3′ requirement) | 5 | 1 | 14 | 0 | F:SDH1.1+R:SDH1.2→NODE_10980_length_27… 878 bp, 2 mm, 3′ 0/20 |


## #1 top-ranked backup

- SDH3 primers under test: F `TCATTAGCAATGGCCTCAGG` (232-251, 53.4 °C) / R `ATCACCTTCATACTGACCGC` (512-531, 53.6 °C), 300 bp

### Stage 1 — six curated transcripts, all 36 forward×reverse combinations

| rule | designed amplicons | other products | foreign-paralog / cross-pair detail |
|---|---|---|---|
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 6 | 2 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8) |
| permissive      (each ≤6 mm, 3′ runs ≥5 nt) | 6 | 5 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8); F:SDH2+R:SDH1.2 on SDH3 571 bp (F 5 mm/3′5, R 7 mm/3′5) *[cross]*; F:SDH3+R:SDH1.2 on SDH3 889 bp (F 0 mm/3′20, R 7 mm/3′5) *[cross]* (+1 more) |
| very permissive (each ≤8 mm, 3′ runs ≥3 nt) | 6 | 38 | F:SDH1.1+R:SDH1.1 on SDH1.1 504 bp (F 8 mm/3′3, R 0 mm/3′20); F:SDH1.1+R:SDH1.1 on SDH1.1 756 bp (F 8 mm/3′4, R 0 mm/3′20); F:SDH1.1+R:SDH1.2 on SDH1.1 428 bp (F 8 mm/3′3, R 8 mm/3′8) *[cross]*; F:SDH1.1+R:SDH1.2 on SDH1.1 60 bp (F 8 mm/3′4, R 8 mm/3′3) *[cross]* (+34 more) |

### Stage 2 — whole transcriptome, every window, same rules as TX_EXHAUSTIVE_SCAN.md

| primer | sites ≤2 mm | ≤4 mm | ≤6 mm | longest 3′ run at ≤4 mm |
|---|---|---|---|---|
| `AAGAGTCTGGTCCATTTCGG` | 2 | 39 | 2995 | 20 |
| `ACCTCAGTCACAGGAATTCG` | 1 | 37 | 2554 | 20 |
| `ACTCAGCAACCCAATAGTGC` | 1 | 51 | 3170 | 20 |
| `ACTTGCTCAGTCCATAGGTG` | 1 | 35 | 2669 | 20 |
| `AGAGGAAGCACTAGGAGGTT` | 1 | 88 | 5681 | 20 |
| `ATAAGCACCTTACTTGAAGC` | 1 | 48 | 3640 | 20 |
| `ATCACCTTCATACTGACCGC` ← **backup** | 2 | 49 | 2934 | 20 |
| `GCAAGCAGATTTCAGAGCTG` | 1 | 82 | 4441 | 20 |
| `GCACTTCGATTAGCCATGGA` | 1 | 42 | 2828 | 20 |
| `GTTCACTGGTGGATTAAGTG` | 1 | 48 | 3616 | 16 |
| `TCATTAGCAATGGCCTCAGG` ← **backup** | 1 | 64 | 4082 | 20 |
| `TGCCACATCAGTGATATCCG` | 1 | 39 | 2866 | 20 |

| rule | designed | cross-pair exact | on-family (paralog) | off-family (unrelated contigs) | worst example |
|---|---|---|---|---|---|
| designed        (both primers 0 mm) | 5 | 1 | 0 | 0 | — |
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 5 | 1 | 0 | 0 | — |
| any             (total ≤8 mm, no 3′ requirement) | 5 | 1 | 4 | 0 | F:SDH1.1+R:SDH1.2→NODE_10980_length_27… 878 bp, 2 mm, 3′ 0/20 |


## #2 zero foreign sites even at ≤8 mm/3-nt run

- SDH3 primers under test: F `CACCCAGAAGAAGGGATGGT` (1552-1571, 54.7 °C) / R `GCAGACTCTCTCTTATGCGT` (1832-1851, 53.4 °C), 300 bp

### Stage 1 — six curated transcripts, all 36 forward×reverse combinations

| rule | designed amplicons | other products | foreign-paralog / cross-pair detail |
|---|---|---|---|
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 6 | 2 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8) |
| permissive      (each ≤6 mm, 3′ runs ≥5 nt) | 6 | 3 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8); F:SDH2+R:SDH1.2 on SDH3 571 bp (F 5 mm/3′5, R 7 mm/3′5) *[cross]* |
| very permissive (each ≤8 mm, 3′ runs ≥3 nt) | 6 | 33 | F:SDH1.1+R:SDH1.1 on SDH1.1 504 bp (F 8 mm/3′3, R 0 mm/3′20); F:SDH1.1+R:SDH1.1 on SDH1.1 756 bp (F 8 mm/3′4, R 0 mm/3′20); F:SDH1.1+R:SDH1.2 on SDH1.1 428 bp (F 8 mm/3′3, R 8 mm/3′8) *[cross]*; F:SDH1.1+R:SDH1.2 on SDH1.1 60 bp (F 8 mm/3′4, R 8 mm/3′3) *[cross]* (+29 more) |

### Stage 2 — whole transcriptome, every window, same rules as TX_EXHAUSTIVE_SCAN.md

| primer | sites ≤2 mm | ≤4 mm | ≤6 mm | longest 3′ run at ≤4 mm |
|---|---|---|---|---|
| `AAGAGTCTGGTCCATTTCGG` | 2 | 39 | 2995 | 20 |
| `ACCTCAGTCACAGGAATTCG` | 1 | 37 | 2554 | 20 |
| `ACTCAGCAACCCAATAGTGC` | 1 | 51 | 3170 | 20 |
| `ACTTGCTCAGTCCATAGGTG` | 1 | 35 | 2669 | 20 |
| `AGAGGAAGCACTAGGAGGTT` | 1 | 88 | 5681 | 20 |
| `ATAAGCACCTTACTTGAAGC` | 1 | 48 | 3640 | 20 |
| `CACCCAGAAGAAGGGATGGT` ← **backup** | 1 | 73 | 4039 | 20 |
| `GCAAGCAGATTTCAGAGCTG` | 1 | 82 | 4441 | 20 |
| `GCACTTCGATTAGCCATGGA` | 1 | 42 | 2828 | 20 |
| `GCAGACTCTCTCTTATGCGT` ← **backup** | 1 | 36 | 2643 | 20 |
| `GTTCACTGGTGGATTAAGTG` | 1 | 48 | 3616 | 16 |
| `TGCCACATCAGTGATATCCG` | 1 | 39 | 2866 | 20 |

| rule | designed | cross-pair exact | on-family (paralog) | off-family (unrelated contigs) | worst example |
|---|---|---|---|---|---|
| designed        (both primers 0 mm) | 5 | 1 | 0 | 0 | — |
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 5 | 1 | 0 | 0 | — |
| any             (total ≤8 mm, no 3′ requirement) | 5 | 1 | 5 | 0 | F:SDH1.1+R:SDH1.2→NODE_10980_length_27… 878 bp, 2 mm, 3′ 0/20 |


## #3 highest-Tm backup (hot annealing)

- SDH3 primers under test: F `TCCCAAGTCCCAACCTACCT` (45-64, 55.4 °C) / R `GTGGTGCCCATCAATGGAGT` (325-344, 56.4 °C), 300 bp

### Stage 1 — six curated transcripts, all 36 forward×reverse combinations

| rule | designed amplicons | other products | foreign-paralog / cross-pair detail |
|---|---|---|---|
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 6 | 1 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]* |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 6 | 2 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8) |
| permissive      (each ≤6 mm, 3′ runs ≥5 nt) | 6 | 3 | F:SDH1.2+R:SDH1.1 on SDH1.2 314 bp (F 0 mm/3′20, R 0 mm/3′20) *[cross]*; F:SDH1.2+R:SDH1.2 on SDH1.2 238 bp (F 0 mm/3′20, R 8 mm/3′8); F:SDH2+R:SDH1.2 on SDH3 571 bp (F 5 mm/3′5, R 7 mm/3′5) *[cross]* |
| very permissive (each ≤8 mm, 3′ runs ≥3 nt) | 6 | 37 | F:SDH1.1+R:SDH1.1 on SDH1.1 504 bp (F 8 mm/3′3, R 0 mm/3′20); F:SDH1.1+R:SDH1.1 on SDH1.1 756 bp (F 8 mm/3′4, R 0 mm/3′20); F:SDH1.1+R:SDH1.2 on SDH1.1 428 bp (F 8 mm/3′3, R 8 mm/3′8) *[cross]*; F:SDH1.1+R:SDH1.2 on SDH1.1 60 bp (F 8 mm/3′4, R 8 mm/3′3) *[cross]* (+33 more) |

### Stage 2 — whole transcriptome, every window, same rules as TX_EXHAUSTIVE_SCAN.md

| primer | sites ≤2 mm | ≤4 mm | ≤6 mm | longest 3′ run at ≤4 mm |
|---|---|---|---|---|
| `AAGAGTCTGGTCCATTTCGG` | 2 | 39 | 2995 | 20 |
| `ACCTCAGTCACAGGAATTCG` | 1 | 37 | 2554 | 20 |
| `ACTCAGCAACCCAATAGTGC` | 1 | 51 | 3170 | 20 |
| `ACTTGCTCAGTCCATAGGTG` | 1 | 35 | 2669 | 20 |
| `AGAGGAAGCACTAGGAGGTT` | 1 | 88 | 5681 | 20 |
| `ATAAGCACCTTACTTGAAGC` | 1 | 48 | 3640 | 20 |
| `GCAAGCAGATTTCAGAGCTG` | 1 | 82 | 4441 | 20 |
| `GCACTTCGATTAGCCATGGA` | 1 | 42 | 2828 | 20 |
| `GTGGTGCCCATCAATGGAGT` ← **backup** | 1 | 47 | 3363 | 20 |
| `GTTCACTGGTGGATTAAGTG` | 1 | 48 | 3616 | 16 |
| `TCCCAAGTCCCAACCTACCT` ← **backup** | 2 | 73 | 3641 | 20 |
| `TGCCACATCAGTGATATCCG` | 1 | 39 | 2866 | 20 |

| rule | designed | cross-pair exact | on-family (paralog) | off-family (unrelated contigs) | worst example |
|---|---|---|---|---|---|
| designed        (both primers 0 mm) | 5 | 1 | 0 | 0 | — |
| strict          (total ≤3 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| competent       (total ≤6 mm, 3′ runs ≥5 nt) | 5 | 1 | 0 | 0 | — |
| long anchor     (total ≤8 mm, 3′ runs ≥8 nt) | 5 | 1 | 0 | 0 | — |
| any             (total ≤8 mm, no 3′ requirement) | 5 | 1 | 4 | 0 | F:SDH1.1+R:SDH1.2→NODE_10980_length_27… 878 bp, 2 mm, 3′ 0/20 |


Full product tables: `backup_curated_products.tsv` (six transcripts) and `backup_tx_products.tsv` (whole transcriptome).
