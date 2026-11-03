# Final Cornus DQD/SDH Primer Set (de novo)

All primers are 20 nt unless otherwise noted. Product sizes are spread
120 → 420 bp in ~60 bp steps so the six genes can be resolved together
on a 2% agarose gel or by fragment analyzer.

## Primer pairs

| Gene | Accession | Direction | Sequence (5'→3') | Region | Position | GC% | Tm(°C) | Amplicon bp |
|---|---|---|---|---|---|---|---|---|
| CserDQD/SDH1.1 | PZ518666 | F | GTTCACTGGTGGATTAAGTG | CDS | 1570-1589 | 45.0 | 49.5 | 120 |
|  |  | R | ATAAGCACCTTACTTGAAGC | 3UTR | 1670-1689 | 40.0 | 48.6 |  |
| CserDQD/SDH1.2 | PZ518667 | F | GCAAGCAGATTTCAGAGCTG | 3UTR | 2268-2287 | 50.0 | 53.1 | 180 |
|  |  | R | ACTCAGCAACCCAATAGTGC | 3UTR | 2428-2447 | 50.0 | 54.0 |  |
| CserDQD/SDH2 | PZ518668 | F | GCACTTCGATTAGCCATGGA | CDS | 492-511 | 50.0 | 54.0 | 240 |
|  |  | R | TGCCACATCAGTGATATCCG | CDS | 712-731 | 50.0 | 54.0 |  |
| CserDQD/SDH3 | PZ518669 | F | TGGAGCCGACATTGTGAAGA | CDS | 738-757 | 50.0 | 55.2 | 300 |
|  |  | R | TGACCGACAGGCTTTCCAAT | CDS | 1018-1037 | 50.0 | 55.1 |  |
| CserDQD/SDH4 | PZ518670 | F | AGAGGAAGCACTAGGAGGTT | CDS | 1216-1235 | 50.0 | 53.1 | 360 |
|  |  | R | AAGAGTCTGGTCCATTTCGG | CDS | 1556-1575 | 50.0 | 52.8 |  |
| CserDQD/SDH5 | PZ518671 | F | ACTTGCTCAGTCCATAGGTG | CDS | 1048-1067 | 50.0 | 53.0 | 420 |
|  |  | R | ACCTCAGTCACAGGAATTCG | CDS | 1448-1467 | 50.0 | 53.2 |  |

## Specificity

In-silico exact-match PCR against all six Cornus DQD/SDH transcripts:

- ✓ No pair amplifies a non-target Cornus DQD/SDH transcript.

## Primer-dimer check

- ✓ No 8-nt 3'-complementarity among the 12 primers.

## Notes on design logic

- **Both-unique pairs** (genes SDH1.2, 2, 3, 4, 5): both forward and reverse 20-mers exist only in that transcript and in no other Cornus DQD/SDH transcript. This gives true gene-specific amplification.
- **SDH1.1 special case**: SDH1.1 owns only 38 unique nucleotides (1567-1604, 19 unique 20-mer starts; see `data/processed/unique_20mer_spans.tsv`) near the CDS 3' end. Because that window is too narrow to hold two 20-mers 120+ bp apart, the SDH1.1 pair uses a *unique forward anchor* (TGTTCACTGGTGGATTAAGT, inside the 38-nt unique window) plus a *shared reverse* that also exists in SDH1.2. Specificity still holds for SDH1.1 because extension cannot be initiated from SDH1.2 without the unique forward primer, but expect SDH1.1 to be the only clean amplicon in this pair.

## Why only ONE primer pair per gene (2, 3, 4, 5)?

- The paralogs are >90% identical across most of their CDS, so most 20-mers are shared. Only short patches (~20-40 nt each) in the 5'UTR, near variable codons in the CDS, and in the 3'UTR are truly unique per gene.
- Once a unique 20-mer is chosen as the gene-specific anchor, a second unique 20-mer must lie 120-500 bp away on the same transcript to complete a pair. That distance filter leaves a small set of options, usually only a handful per gene, so the pipeline returns one high-quality pair per gene rather than many alternatives.
- Additional pairs *could* be generated per gene by relaxing either the both-unique rule (allowing one shared partner) or the size spacing rule; if you want backup pairs, we can rerun with those relaxed settings.
