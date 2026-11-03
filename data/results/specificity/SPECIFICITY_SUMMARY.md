# Specificity of the final six-pair DQD/SDH primer set — full check summary

Date: 2026-10-29. All analyses ran **locally** (blastn 2.17.0 + numpy); nothing was uploaded.
This document is an *additional* validation layer. The design deliverables
(`FINAL_primer_pairs.tsv`, `FINAL_primers.fasta`, `FINAL_primers_report.md`)
were **not modified**.

Templates used

| Source | Content |
|---|---|
| `DQD_SDH_possible_nt.fasta` | the six curated transcripts PZ518666–PZ518671 (12,888 bp) — the wet-lab templates |
| `Transcriptome/final.clust_transcripts_longest_iso.fasta` | whole *C. sericea* transcriptome, 29,124 transcripts, 48.1 Mb |
| `annotated_final_nr.gff3` / `Transcriptome/final.clust_annotation_longest_iso.gff3` | contig annotation (nr + UniProt best hits) |

---

## Check 1 — BLAST-based, six curated transcripts
`run_primer_blast_check.py` → `PRIMER_BLAST_CHECK.md`

`blastn -task blastn-short -word_size 7 -perc_identity 75`, product 80–600 bp,
full-length primer alignment, ≤3 total mismatches, perfect 5-nt 3′ terminal match.

| Pair | Verdict | Product | Off-target products |
|---|---|---|---|
| CserDQD_SDH1.1 | SPECIFIC | 120 bp | none |
| CserDQD_SDH1.2 | SPECIFIC | 180 bp | none |
| CserDQD_SDH2 | SPECIFIC | 240 bp | none |
| CserDQD_SDH3 | SPECIFIC | 300 bp | none |
| CserDQD_SDH4 | SPECIFIC | 360 bp | none |
| CserDQD_SDH5 | SPECIFIC | 420 bp | none |

## Check 2 — BLAST-based, whole transcriptome
`run_transcriptome_check.py` → `TRANSCRIPTOME_CHECK.md`, `tx_blast_hits.tsv`

52,371 raw primer hits over 29,124 transcripts; only **11** were full-length 20/20
alignments and **all 11 were 100 % identity at the designed positions on the five
target contigs**. Candidate products outside the DQD/SDH loci: **0**.

Megablast mapping of the curated transcripts onto the assembly:

| Curated transcript | Contig | Identity | Aligned |
|---|---|---|---|
| CserDQD/SDH1.1 | NODE_10980_length_2779 (g1243_i5) | 100 % | 1,585 bp (84 % of the transcript, in 2 blocks) |
| CserDQD/SDH1.2 | NODE_10980_length_2779 | 100 % | 2,779 bp (full) |
| CserDQD/SDH2 | NODE_18517_length_2241 | 100 % | full |
| CserDQD/SDH3 | NODE_20612_length_2133 | 100 % | full |
| CserDQD/SDH4 | NODE_21648_length_2083 | 100 % | full |
| CserDQD/SDH5 | NODE_29040_length_1767 | 100 % | full |

Consequence: SDH1.1 has **no separate contig** — CD-HIT/Trinity collapsed it with
SDH1.2 (they are identical over the first 1,585 nt), so the SDH1.1 pair produces no
product inside the assembly and had to be judged on the curated sequence (Checks 1 and 4).

The assembly contains exactly **five** bifunctional DQD/SDH-annotated contigs
(= the six genes) plus one 3-dehydroquinate **synthase** contig (NODE_24349, DHQS,
a pathway neighbour). No primer has any site with ≤6 mismatches on DHQS.

## Check 3 — Exhaustive seed-free scan, whole transcriptome
`tx_exhaustive_scan.py` → `TX_EXHAUSTIVE_SCAN.md`, `exhaustive_sites.tsv`

BLAST needs an exact 7-mer seed, so a mismatch-rich site can escape it. This scan
compared every 20-nt window: 47.6 M valid positions × 12 primers × 2 strands,
no gaps, non-ACGT treated as never-matching, boundary-spanning windows discarded.

Sites per primer (cumulative, whole transcriptome):

| mismatches over 20 nt | sites |
|---|---|
| 0 | 11 (the designed sites, all on the 5 target contigs) |
| ≤1 | 12 |
| ≤2 | 13 |
| ≤3 | 64 |
| ≤4 | 618 |
| ≤5 | 5,677 |
| ≤6 | 41,437 |

Candidate amplicons (same contig, 60–1500 bp, ≤8 total mismatches): **10**, of which

* 7 on the target/family contigs = the 5 designed amplicons + 2 marginal variants,
* 3 on other DQD/SDH contigs (SDH3 pair on the SDH4 and SDH2 contigs, SDH4 pair on the SDH2 contig),
* **0 on any unrelated transcript**.

Relaxing the product filter further (BLAST variants A/B/C in `TX_SENSITIVITY.md`:
≤4 / ≤6 / ≤8 mismatches, 3′-terminal match 5 / 4 / 3 nt, 60–1500 bp) still returns
only the five expected amplicons — never an unrelated-contig product.

## Check 4 — Exhaustive seed-free scan, six curated transcripts
`curated_exhaustive_scan.py` → `TX_CURATED_EXHAUSTIVE.md`, `curated_exhaustive_sites.tsv`

Every window of every primer against all six transcripts, best window tabulated with
mismatch positions and the length of the perfectly matched run counting back from the
primer 3′ end.

| Rule set | Products | Cross-products |
|---|---|---|
| strict — ≤3 mm total, perfect 5-nt 3′ end, 80–600 bp | 6 (the designed ones) | 0 |
| medium — ≤4 mm total, perfect 8-nt 3′ end, 60–1000 bp | 6 | 0 |
| loose — ≤6 mm total, perfect 8-nt 3′ end, 60–1000 bp | 6 | 0 |
| permissive — any mm ≤6, 3′ end ≥5 nt, 60–1000 bp | — | 4 (see below) |

Each pair yields exactly **one** product on its own transcript (120/180/240/300/360/420 bp),
with no alternative internal sites — so no spuriously sized bands and no primer-dimer-driven
side products on the targets.

The four permissive cross-pairings, all originating from the **SDH3 pair**:

| Pair | Non-target transcript | Size | mm (F/R) | 3′ run (F/R) |
|---|---|---|---|---|
| SDH3 | CserDQD/SDH4 | **300 bp** | 3 / 3 | 7 / 8 |
| SDH3 | CserDQD/SDH2 | 309 bp | 3 / 4 | 7 / 5 |
| SDH3 | CserDQD/SDH4 | 630 bp | 6 / 3 | 6 / 8 |
| SDH3 | CserDQD/SDH2 | 636 bp | 6 / 4 | 6 / 5 |

---

## Two things to keep in mind

**1. SDH1.1 vs SDH1.2 is decided by the 3′-terminal nucleotide.**
`CserDQD_SDH1.1_Fwd` (GTTCACTGGTGGATTAAGTG) aligns to SDH1.2 with only 2 mismatches,
at primer positions 17 and **20**:

```
SDH1.1 1570-1589 : GTTCACTGGTGGATTAAGTG   (perfect, 3' run 20)
SDH1.2 1570-1589 : GTTCACTGGTGGATTAGGTA   (mm at 17 A/G and at 20 G/A, 3' run 0)
```

A mismatch on the very last base is the most damaging kind (no 3′→5′ cleaning in Taq),
and there is a second safety margin: the shared reverse primer sits 1,012 bp away on
SDH1.2, so even if that terminal mismatch were extended, the product would be ~1.0 kb,
not 120 bp — it separates on the gel and is badly amplified in qPCR. Still, SDH1.1 is
the least robust member of the set: confirm its band by sequencing, and treat the
120 bp band as SDH1.1 only if the ~1 kb smear/band is absent.

**2. SDH3 primers can, at low stringency, co-amplify SDH4 (same 300 bp) and SDH2 (309 bp).**
These sites carry 3–4 mismatches each but keep 7–8 nt perfect 3′ runs, so they are rejected
by every standard Primer-BLAST rule yet are not physically impossible at a soft annealing
temperature. The SDH3 pair is also the highest-Tm pair of the set (55.2 / 55.1 °C,
SantaLucia NN), so it tolerates heat better than the others:

* multiplex gel: anneal ≥58 °C rather than 55 °C (the SDH1.1 pair at 49.5/48.6 °C is the
  weak link — if its band fades at 58 °C, re-order that pair with a slightly longer anchor);
* SDH3 singleplex qPCR: 58–60 °C annealing, plus a melt-curve check (the SDH4 co-product
  has a different sequence and should show a distinct Tm);
* definitive test: excise the 300 bp band and sequence it.

## Scope limits (stated plainly)

* The search space is the **transcriptome**, which is correct for RT-PCR/RT-qPCR on cDNA.
  Genomic DNA, introns and intergenic regions are not covered, so a no-RT control may still
  show a product if the amplicon region does not span an intron.
* The assembly is clustered to longest isoforms; allelic variants and low-coverage
  transcripts collapsed out of it are under-represented (this is exactly why SDH1.1 shows
  no contig and had to be assessed on the curated sequence).
* Non-ACGT positions (sequencing errors, `N`) were counted as mismatches, i.e. treated
  conservatively.

## Files

| File | Content |
|---|---|
| `run_primer_blast_check.py`, `PRIMER_BLAST_CHECK.md`, `blast_hits.tsv`, `product_check.tsv` | Check 1 |
| `run_transcriptome_check.py`, `TRANSCRIPTOME_CHECK.md`, `tx_blast_hits.tsv`, `tx_product_check.tsv`, `tx_priming_sites.tsv`, `tx_selfmap.tsv` | Check 2 |
| `tx_sensitivity_check.py`, `TX_SENSITIVITY.md` | threshold stress test of Check 2 |
| `tx_exhaustive_scan.py`, `TX_EXHAUSTIVE_SCAN.md`, `exhaustive_sites.tsv`, `exhaustive_products.tsv` | Check 3 |
| `curated_exhaustive_scan.py`, `TX_CURATED_EXHAUSTIVE.md`, `curated_exhaustive_sites.tsv`, `curated_exhaustive_products.tsv` | Check 4 |
| `cornus_dqdsdh_6tx.n*`, `cornus_tx.n*` | BLAST databases (regenerable) |
