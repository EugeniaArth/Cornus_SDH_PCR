# 02 — Sequence survey: where the six transcripts differ

**Step scripts:** `scripts/step2_survey/extract_cornus_utrs.py`, `scripts/step2_survey/unique_20mer_spans.py`
**Outputs:** `data/processed/cornus_5utr.fasta`, `cornus_3utr.fasta`, `cornus_cds.fasta`,
`cornus_utr_summary.tsv`, `cornus_unique_5utr_kmers.tsv`, `cornus_unique_3utr_kmers.tsv`,
`unique_20mer_spans.tsv`

The earlier variable‑positions analysis (compare to *Arabidopsis* SDH) was based on alignment CSV files that have been moved to `archive/`. It is not required for primer selection and is omitted from the current pipeline; the uniqueness of 20‑mers in UTRs provides sufficient discrimination.

## Why survey before designing

A gene-specific primer has to sit in a region that is unique to its target. With six paralogs
this similar (the bifunctional DQD/SDH domain is close to identical across them, which is why
they were picked out of the assembly as a family), most of the CDS is shared, so the informative
territory is mostly untranslated. The survey cut the transcripts into UTR and CDS compartments
and counted how many 20-mers are unique to one gene, which is the raw material the design step
draws from.

## Transcript layout

| Gene | Transcript | CDS | 5' UTR | 3' UTR |
|---|---|---|---|---|
| CserDQD/SDH1.1 | 1885 | 1497 | 130 | 258 |
| CserDQD/SDH1.2 | 2779 | 1503 | 130 | 1146 |
| CserDQD/SDH2 | 2241 | 1566 | 224 | 451 |
| CserDQD/SDH3 | 2133 | 1599 | 240 | 294 |
| CserDQD/SDH4 | 2083 | 1584 | 133 | 366 |
| CserDQD/SDH5 | 1767 | 1560 | 112 | 95 |

## Unique 20-mers per compartment

Count of 20-mers that occur in exactly one of the six transcripts:

| Gene | 5' UTR | 3' UTR |
|---|---|---|
| CserDQD/SDH1.1 | 0 | 0 |
| CserDQD/SDH1.2 | 0 | 888 |
| CserDQD/SDH2 | 205 | 432 |
| CserDQD/SDH3 | 221 | 275 |
| CserDQD/SDH4 | 114 | 347 |
| CserDQD/SDH5 | 93 | 76 |

Two findings shaped the whole design:

* **SDH1.1 and SDH1.2 have identical 5' UTRs** (both 130 nt, byte-identical), so neither has any
  unique 5'-UTR 20-mer. Any primer placed there amplifies both.
* **SDH1.1 has no unique 20-mer in either UTR.** Its 258-nt 3' UTR is contained inside the
  1146-nt 3' UTR of SDH1.2. The only sequence that distinguishes SDH1.1 is a short block near the
  3' end of the CDS.

The amino-acid survey (`analyze_variable_positions.py`) confirmed this from the protein side: the
six genes differ at scattered positions, and the windows that are both variable and unique in
nucleotide terms number 72 (`variable_unique_windows.tsv`). Most are too narrow, or too close
together, to host a primer pair.

## Unique territory, per transcript

`unique_20mer_spans.py` counts 20-mers absent from the other five transcripts (reverse complements
included, because a shared site can sit on either strand), and reports where they lie:

| Gene | Transcript | Unique 20-mers | Span | Longest contiguous stretch |
|---|---|---|---|---|
| CserDQD/SDH1.1 | 1885 | 19 | 1567-1604 (38 nt) | 19 starts |
| CserDQD/SDH1.2 | 2779 | 913 | 1567-2779 (1213 nt) | 911 starts |
| CserDQD/SDH2 | 2241 | 2208 | whole transcript | 1535 starts |
| CserDQD/SDH3 | 2133 | 2110 | whole transcript | 740 starts |
| CserDQD/SDH4 | 2083 | 2050 | whole transcript | 728 starts |
| CserDQD/SDH5 | 1767 | 1747 | whole transcript | 982 starts |

Table: `data/processed/unique_20mer_spans.tsv`.

## Consequence for SDH1.1

SDH1.1 owns exactly 38 nucleotides: positions 1567-1604, just upstream of the CDS end at 1627.
Inside that window 19 different 20-mers start. Two primers 120 bp apart cannot fit into 38 nt, so
a fully unique SDH1.1 pair is impossible on these transcripts. The design settled on a unique
forward anchor inside that window (1570-1589) plus a reverse primer whose site also exists in
SDH1.2 (see report 03). This is a property of the data, not an oversight, and it is the reason
SDH1.1 is the pair to check first on a gel.

## Reproduce

```sh
python3 scripts/step2_survey/extract_cornus_utrs.py
python3 scripts/step2_survey/analyze_variable_positions.py
python3 scripts/step2_survey/unique_20mer_spans.py
```
