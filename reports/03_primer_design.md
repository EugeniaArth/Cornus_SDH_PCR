# 03 — Primer design and the final six-pair set

**Step scripts:** `scripts/step3_design/design_final_pairs.py`, `scripts/step3_design/validate_final_primers.py`
**Outputs:** `data/results/FINAL_primer_pairs.tsv`, `data/results/FINAL_primers.fasta`,
`data/results/VALIDATION_auto_report.md`

## Design rules

| Property | Requirement | Achieved |
|---|---|---|
| Length | 20 nt for every primer | 20 nt, all 12 |
| Uniqueness | primer must occur in one transcript only | 5 pairs both-unique; SDH1.1 by unique anchor (see below) |
| GC | 40-50 % | 40-50 % |
| Tm | nearest-neighbour, matched within a pair | 48.6-55.2 °C, ΔTm ≤ 1.6 °C |
| Amplicon | one distinct size per gene, resolvable on a gel | 120/180/240/300/360/420 bp, 60 bp apart |
| Dimers | no 8-nt 3'-complementarity between any two primers | none among the 12 primers |
| Cross-amplification | no pair product on a non-target transcript | none at 0 mismatches (see report 04 for relaxed settings) |

Melting temperature is calculated with the SantaLucia 1998 unified nearest-neighbour table
(enthalpy/entropy sums, salt correction for ~50 mM Na⁺, 250 nM primer). The Marmur-Doty
approximation is also reported in the tables for comparison, but selection used the
nearest-neighbour value.

## How candidates were generated

`design_final_pairs.py` enumerates every 20-mer of all six transcripts, keeps only those absent
from the other five (reverse complement included), then searches for forward/reverse pairs that
face each other inside the required amplicon window. Candidates are scored on Tm match, GC,
absence of runs and dimers, and on keeping the six products on the 60 bp ladder. The size ladder
is what makes this a usable multiplex: one tube, six products, each identifiable by length.

Five genes (SDH1.2, SDH2, SDH3, SDH4, SDH5) yielded both-unique pairs. SDH1.1 did not, and could
not: its entire unique territory is 38 nt at positions 1567-1604 (report 02), which cannot hold
two 20-mers 120 bp apart. The fallback the script allows is one unique anchor plus one shared
partner.

## The set

| Gene | Accession | Primer | Sequence (5'→3') | Pos | Region | GC% | Tm °C | Product |
|---|---|---|---|---|---|---|---|---|
| CserDQD/SDH1.1 | PZ518666 | F | GTTCACTGGTGGATTAAGTG | 1570-1589 | CDS | 45.0 | 49.5 | 120 bp |
| | | R | ATAAGCACCTTACTTGAAGC | 1670-1689 | 3'UTR | 40.0 | 48.6 | |
| CserDQD/SDH1.2 | PZ518667 | F | GCAAGCAGATTTCAGAGCTG | 2268-2287 | 3'UTR | 50.0 | 53.1 | 180 bp |
| | | R | ACTCAGCAACCCAATAGTGC | 2428-2447 | 3'UTR | 50.0 | 54.0 | |
| CserDQD/SDH2 | PZ518668 | F | GCACTTCGATTAGCCATGGA | 492-511 | CDS | 50.0 | 54.0 | 240 bp |
| | | R | TGCCACATCAGTGATATCCG | 712-731 | CDS | 50.0 | 54.0 | |
| CserDQD/SDH3 | PZ518669 | F | TGGAGCCGACATTGTGAAGA | 738-757 | CDS | 50.0 | 55.2 | 300 bp |
| | | R | TGACCGACAGGCTTTCCAAT | 1018-1037 | CDS | 50.0 | 55.1 | |
| CserDQD/SDH4 | PZ518670 | F | AGAGGAAGCACTAGGAGGTT | 1216-1235 | CDS | 50.0 | 53.1 | 360 bp |
| | | R | AAGAGTCTGGTCCATTTCGG | 1556-1575 | CDS | 50.0 | 52.8 | |
| CserDQD/SDH5 | PZ518671 | F | ACTTGCTCAGTCCATAGGTG | 1048-1067 | CDS | 50.0 | 53.0 | 420 bp |
| | | R | ACCTCAGTCACAGGAATTCG | 1448-1467 | CDS | 50.0 | 53.2 | |

Ordered DNA (desal) is sufficient for routine PCR; no 5' tails or blocked bases are used, so the
multiplex stays cheap.

## SDH1.1 read carefully

The SDH1.1 reverse site (1670-1689) also exists in SDH1.2 — that is why the two rows can look like
they belong to different genes when you line them up against the sequences. Specificity comes from
the forward primer: extension can only start where the unique SDH1.1 forward anchor anneals, so
SDH1.2 contributes nothing on its own. In practice SDH1.1 is the band to trust least without
confirmation, and the one to sequence first.

## Validation

`validate_final_primers.py` recomputes length, GC and both Tm estimates, checks every primer
against all six transcripts for non-target products, and checks all 12×12 pairwise 3'
complementarity. Result (`data/results/VALIDATION_auto_report.md`):

* no pair amplifies a non-target Cornus DQD/SDH transcript at 0 mismatches;
* no 8-nt 3'-complementarity among the 12 primers;
* all 12 primers inside the GC and Tm targets.

## Known weak point

Under relaxed matching (≤6 mismatches per primer with a short perfect 3' run) the SDH3 pair can
also prime on SDH4 and SDH2. That finding is documented in report 04, and report 05 gives a
replacement SDH3 pair that removes it. The other five pairs show no such behaviour.

## Reproduce

```sh
python3 scripts/step3_design/design_final_pairs.py      # regenerates the set byte-for-byte
python3 scripts/step3_design/validate_final_primers.py
```
