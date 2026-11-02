# Cornus sericea DQD/SDH Primer Design

Six gene‑specific PCR primer pairs designed to distinguish the six *Cornus sericea* DQD/SDH genes (`CserDQD/SDH1.1`, `1.2`, `2`, `3`, `4`, `5`). All primers are 20 nt, produce amplicons spaced 60 bp apart (120 bp → 420 bp) for multiplex gel resolution, and have been validated for specificity across the whole *Cornus* transcriptome (29 k transcripts, 48 Mb).

## Quick links

| File | Description |
|------|-------------|
| [`FINAL_primer_pairs.tsv`](FINAL_primer_pairs.tsv) | Machine‑readable primer table (gene, sequence, position, GC %, Tm, amplicon size). |
| [`FINAL_primers.fasta`](FINAL_primers.fasta) | Twelve primers in FASTA format, ready for ordering‑portal bulk upload. |
| [`FULL_REPORT.md`](FULL_REPORT.md) | Comprehensive design report (goal, input data, pipeline, specificity checks, backup SDH3 pair, experimental recommendations). |
| [`data/results/`](data/results/) | Copies of the final primer files. |
| [`data/input/`](data/input/) | Input sequences, coordinates, transcriptome assembly. |
| [`scripts/`](scripts/) | Reproducible Python scripts for primer generation, validation, and specificity checking. |
| [`reports/`](reports/) | Detailed validation reports (including four‑layer specificity checks). |
| [`backup/`](backup/) | Intermediate files and older versions (archived). |

## Primer sequences at a glance

| Gene | Forward (5′→3′) | Reverse (5′→3′) | Amplicon (bp) |
|------|-----------------|-----------------|---------------|
| SDH1.1 | `GTTCACTGGTGGATTAAGTG` | `ATAAGCACCTTACTTGAAGC` | 120 |
| SDH1.2 | `GCAAGCAGATTTCAGAGCTG` | `ACTCAGCAACCCAATAGTGC` | 180 |
| SDH2 | `GCACTTCGATTAGCCATGGA` | `TGCCACATCAGTGATATCCG` | 240 |
| SDH3 | `TGGAGCCGACATTGTGAAGA` | `TGACCGACAGGCTTTCCAAT` | 300 |
| SDH4 | `AGAGGAAGCACTAGGAGGTT` | `AAGAGTCTGGTCCATTTCGG` | 360 |
| SDH5 | `ACTTGCTCAGTCCATAGGTG` | `ACCTCAGTCACAGGAATTCG` | 420 |

**Tm:** SantaLucia 1998 nearest‑neighbour values at 250 nM primer, 50 mM Na⁺.

## Key findings

1. **Specificity:** Zero off‑target priming in the whole transcriptome under competent rules (≤6 mismatches, ≥5‑nt perfect 3′‑run).
2. **SDH3 cross‑amplification risk:** The original SDH3 pair can co‑amplify SDH4 (300 bp) and SDH2 (309 bp) under permissive conditions. A **backup SDH3 pair** with zero competent foreign sites is provided (`CACCCAGAAGAAGGGATGGT` / `GCAGACTCTCTCTTATGCGT`).
3. **SDH1.1 uniqueness:** SDH1.1 shares most of its sequence with SDH1.2; the forward primer is unique, the reverse primer is shared but cannot amplify SDH1.2 because the forward primer is absent there.
4. **Multiplex‑ready:** Amplicons spaced ~60 bp apart; all six bands resolve on 2 % agarose.

## How to use

1. **Order primers** using `FINAL_primers.fasta` (or replace SDH3 with the backup pair if desired).
2. **PCR setup:**
   - Multiplex gel: 55–58 °C annealing, 2 % agarose, 120 V for 60 min.
   - Single‑plex RT‑qPCR: 54–56 °C annealing, SYBR/EvaGreen chemistry.
3. **Validation:**
   - Confirm each single‑plex yields only the expected band on mixed cDNA.
   - Sequence the SDH1.1 120 bp band (distinguish from potential 1 kb SDH1.2 product).
   - If using the original SDH3 pair, sequence the 300 bp band to check for SDH4 co‑amplification.

## Reproducing the analysis

See [`FULL_REPORT.md`](FULL_REPORT.md) section “Appendix: How to reproduce” for command‑line instructions.

## License and citation

This work is provided as‑is for academic use. If you use these primers in a publication, please cite the GenBank accessions (PZ518666–PZ518671) and mention that primer design was performed with a custom Python pipeline (available in `scripts/`).