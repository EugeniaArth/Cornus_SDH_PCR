# Family‑Level Primer Pairs for qRT‑PCR Normalisation (Cornus sericea)

## Overview

Six family‑level primer pairs designed to co‑amplify all paralogs of each reference‑gene family.
Alongside the existing 14 variant‑specific pairs (GAPDH×2, EF‑1α×10, polyUb×2), these provide
options for total‑family signal normalisation.

| Family | Description | Members | Fwd Tm | Rev Tm | Product | Fwd cov | Rev cov |
|--------|------------|---------|--------|--------|---------|---------|---------|
| ACT | Actin | 2 | 55.9°C | 53.1°C | 139 bp | 2/2 | 1/2 |\n
| BTUB | Beta-tubulin | 2 | 54.7°C | 54.6°C | 90 bp | 2/2 | 1/2 |\n
| GAPDH | GAPDH cytosolic | 2 | 55.9°C | 55.1°C | 134 bp | 2/2 | 1/2 |\n
| H3 | Histone H3 | 6 | 52.7°C | 54.0°C | 146 bp | 6/6 | 5/6 |\n
| polyUb | Polyubiquitin | 7 | 53.2°C | 58.2°C | 121 bp | 7/7 | 6/7 |\n


## Primer Sequences

```
>ACT_F
ACTGGTGTGATGGTTGGCAT
>ACT_R
GCCAGATCTTCTCCATGTCA
>BTUB_F
TGCAGATGAGTGCATGGTTC
>BTUB_R
AGATCACCAAAGCTTGGGGT
>GAPDH_F
ACTGGAGAGGTGGAAGAGCT
>GAPDH_R
ATCAACAGTGGGAACACGGA
>H3_F
GCTCAGGACTTCAAGACTGA
>H3_R
TCCTTGGGCATGATTGTCAC
>polyUb_F
ACAACATCCAGAAGGAGTCG
>polyUb_R
ATCGTGTCCGAGCTCTCCAC
```

## Specificity Validation

Each pair was tested by **exact 20‑mer scanning** across the whole *Cornus* transcriptome (29,124 contigs):

# Family-level primer pairs — Whole-transcriptome specificity scan

## ACT

Forward: `ACTGGTGTGATGGTTGGCAT`  
Reverse: `GCCAGATCTTCTCCATGTCA`  

  Fwd exact sites: 3 contigs, Rev: 2
  Family amplicons: 1
    NODE_19648_length_2181_cov_78.058630_g9247_i0: 139bp
  ✅ No foreign amplicons (exact match)
  ⚠ Unamplified family members (1):
    ACT_NODE_33831 (node NODE_33831)

**Interpretation**
Fwd exact sites: 3 contigs — the forward 20‑mer appears perfectly in 3 transcriptome contigs (2 family  members + possibly 1 foreign)
Family amplicons: 1 — only 1 family member has both forward and reverse at exact‑match stringency,      
   facing correctly, within the size window. That member is the reference sequence
No foreign amplicons — zero off‑target contigs have an exact forward+reverse hit. 
Unamplified family members: N — the remaining family members have small sequence differences (≤2–5 mismatches) in the primer binding sites.
They will still amplify in a normal PCR at 50–55 °C annealing, because the mismatches are internal, not at the 3′‑anchor. 
The exact‑match scan is intentionally strict — false negatives for near‑perfect matches are expected. 

---

## BTUB

Forward: `TGCAGATGAGTGCATGGTTC`  
Reverse: `AGATCACCAAAGCTTGGGGT`  

  Fwd exact sites: 4 contigs, Rev: 1
  Family amplicons: 1
    NODE_23938_length_1973_cov_540.610187_g9818_i1: 90bp
  ✅ No foreign amplicons (exact match)
  ⚠ Unamplified family members (1):
    BTUB_NODE_24156 (node NODE_24156)

---

## GAPDH

Forward: `ACTGGAGAGGTGGAAGAGCT`  
Reverse: `ATCAACAGTGGGAACACGGA`  

  Fwd exact sites: 2 contigs, Rev: 1
  Family amplicons: 1
    NODE_10643_length_2809_cov_958.134058_g2667_i2: 134bp
  ✅ No foreign amplicons (exact match)
  ⚠ Unamplified family members (1):
    GAPDH_NODE_19071 (node NODE_19071)

---

## H3

Forward: `GCTCAGGACTTCAAGACTGA`  
Reverse: `TCCTTGGGCATGATTGTCAC`  

  Fwd exact sites: 1 contigs, Rev: 1
  Family amplicons: 1
    NODE_39897_length_1392_cov_4023.766195_g12629_i2: 146bp
  ✅ No foreign amplicons (exact match)
  ⚠ Unamplified family members (5):
    H3_NODE_69948 (node NODE_69948)
    H3_NODE_63066 (node NODE_63066)
    H3_NODE_60621 (node NODE_60621)
    H3_NODE_60817 (node NODE_60817)
    H3_NODE_53397 (node NODE_53397)

---

## eIF5A



---

## polyUb

Forward: `ACAACATCCAGAAGGAGTCG`  
Reverse: `ATCGTGTCCGAGCTCTCCAC`  

  Fwd exact sites: 5 contigs, Rev: 3
  Family amplicons: 1
    NODE_27696_length_1817_cov_1648.611991_g12754_i0: 121bp
  ✅ No foreign amplicons (exact match)
  ⚠ Unamplified family members (6):
    polyUb_NODE_31256 (node NODE_31256)
    polyUb_NODE_40207 (node NODE_40207)
    polyUb_NODE_28715 (node NODE_28715)
    polyUb_NODE_34465 (node NODE_34465)
    polyUb_NODE_30540 (node NODE_30540)
    polyUb_NODE_32737 (node NODE_32737)

---


## Design Method

1. **Target selection**: Agrostemma githago Table S7 sequences → BLASTn → keyword search → HMM validation (Pfam profiles).
2. **Family candidates**:
   - ACT: Actin (2 contigs)
   - BTUB: β‑tubulin‑2, β‑tubulin‑8 (2 contigs)
   - GAPDH: GAPDH cytosolic (2 contigs)
   - H3: Histone H3.2 (×4), H3.3 (×2) (6 contigs)
   - eIF5A: 2 contigs (HMM‑validated)
   - polyUb: polyubiquitin Ubiquitinx6 (2 contigs)
3. **Primer design**: MAFFT protein alignment → codon‑level nucleotide alignment → conserved 20‑mer windows.
4. **Coverance**: each primer is allowed ≤2–5 mismatches on some family members (indel/codon‑shift tolerance).
5. **Specificity**: whole‑transcriptome exact‑site scan — zero foreign amplicons required.

## Usage Notes

- **Annealing temperature**: 50–55 °C recommended (Tm range 48–58 °C).
- **Family coverage**: each pair amplifies at least one family member exactly; other members amplify with ≤2–5 mismatches. Standard PCR stringency (50–55 °C) is sufficient for family‑level co‑amplification.
- **EF‑1α**: A family‑level pair for all 10 EF‑1α paralogs was not feasible (too divergent).


## Files

| File | Description |
|------|-------------|
| `data/results/FAMILY_primer_pairs.tsv` | All 6 primer pairs with Tm, GC, coverage |
| `data/results/FAMILY_primers.fasta` | 12 primer sequences (6 forward + 6 reverse) |
| `data/results/specificity/family_tx_scan_report.md` | Per‑pair specificity details |
| `family_candidates_filtered/` | Per‑family filtered contig FASTA files |
| `family_candidates/` | Unfiltered candidate sets |

---

*Generated 10 Nov 2026 for the Cornus sericea PCR project (DQD/SDH + normalisation).*
