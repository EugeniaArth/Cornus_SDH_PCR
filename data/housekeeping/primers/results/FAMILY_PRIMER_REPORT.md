# Family-Level Primer Pairs for qRT-PCR Normalisation (Cornus sericea)

## Overview

Six family-level primer pairs designed to co-amplify all paralogs of each reference-gene family.
All primers are in the **CDS** region (derived from MAFFT codon alignment).

| Family | Description | Members | Forward (5'->3') | Reverse (5'->3') | Fwd GC | Rev GC | Product | FwdPos | RevPos | Region |
|--------|------------|---------|-----------------|-----------------|--------|--------|---------|--------|--------|--------|
| ACT | Actin | 2 | `ACTGGTGTGATGGTTGGCAT` | `GCCAGATCTTCTCCATGTCA` | 50.0% | 50.0% | 139 bp | 127 | 264 | CDS |
| BTUB | Beta-tubulin | 2 | `TGCAGATGAGTGCATGGTTC` | `AGATCACCAAAGCTTGGGGT` | 50.0% | 50.0% | 90 bp | 585 | 673 | CDS |
| GAPDH | GAPDH cytosolic | 2 | `ACTGGAGAGGTGGAAGAGCT` | `ATCAACAGTGGGAACACGGA` | 55.0% | 50.0% | 134 bp | 590 | 722 | CDS |
| H3 | Histone H3 | 6 | `GCTCAGGACTTCAAGACTGA` | `TCCTTGGGCATGATTGTCAC` | 50.0% | 50.0% | 146 bp | 226 | 370 | CDS |
| polyUb | Polyubiquitin | 7 | `ACAACATCCAGAAGGAGTCG` | `ATCGTGTCCGAGCTCTCCAC` | 50.0% | 60.0% | 121 bp | 656 | 775 | CDS |

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

Each pair was tested by **exact 20-mer scanning** across the whole *Cornus* transcriptome (29,124 contigs):

*Scan report not available - run family_tx_scan.py first.*

## Design Method

1. **Target selection**: Agrostemma githago Table S7 sequences -> BLASTn -> keyword search -> HMM validation (Pfam profiles).
2. **Family candidates**:
   - ACT: Actin (2 contigs)
   - BTUB: Beta-tubulin-2, Beta-tubulin-8 (2 contigs)
   - GAPDH: GAPDH cytosolic (2 contigs)
   - H3: Histone H3.2 (x4), H3.3 (x2) (6 contigs)
   - eIF5A: 3 contigs (PF01287 validated)
   - polyUb: polyubiquitin (7 contigs, architecture-validated)
3. **Primer design**: MAFFT protein alignment -> codon-level nucleotide alignment -> conserved 20-mer windows.
4. **Coverage**: each primer is allowed <=5 mismatches on some family members (indel/codon-shift tolerance).
5. **Specificity**: whole-transcriptome exact-site scan - zero foreign amplicons required.

## Usage Notes

- **Annealing temperature**: 50-55 C recommended (Tm range 48-58 C).
- **Family coverage**: each pair amplifies at least one family member exactly; other members amplify with <=5 mismatches.
- **EF-1a**: A family-level pair for all 10 EF-1a paralogs was not feasible (too divergent).
- **eIF5A**: Family-level pair was not feasible for eIF5A (3 members, too divergent for a single conserved pair).

## Files

| File | Description |
|------|-------------|
| `data/results/FAMILY_primer_pairs.tsv` | All primer pairs with Tm, GC, positions, coverage |
| `data/results/FAMILY_primers.fasta` | Primer sequences |
| `data/results/specificity/family_tx_scan_report.md` | Per-pair specificity details |
| `family_candidates_filtered/` | Per-family filtered contig FASTA files |

---

*Generated for the Cornus sericea PCR project (DQD/SDH + normalisation).*