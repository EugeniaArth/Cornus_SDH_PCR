# 04 — Specificity checks (four layers) and the SDH3 finding

**Step scripts:** `scripts/step4_specificity/run_primer_blast_check.py` (Check 1),
`run_transcriptome_check.py` (Check 2), `tx_exhaustive_scan.py` (Check 3),
`tx_sensitivity_check.py` (Check 3b), `curated_exhaustive_scan.py` (Check 4)
**Outputs:** everything in `data/results/specificity/`

All checks ran locally with `blastn` (Bioinformatics Toolkit suite) plus numpy-based positional
scans. No web Primer-BLAST service was used, so nothing was uploaded.

## Why four layers

A BLAST search answers "does this primer have a near-perfect match elsewhere". It does not answer
"can this primer pair, with a few mismatches and a short good 3' end, still make a product".
Those are different questions, and the second one needs an exhaustive scan of every possible
forward×reverse combination. The layers go from fast-and-strict to slow-and-exhaustive:

| Check | Scope | Method | What it can see |
|---|---|---|---|
| 1 | 6 curated transcripts | blastn + Primer-BLAST-style product filter (full-length primer, 80-600 bp, ≤3 total mismatches, perfect 5-nt 3' end) | near-perfect off-target pairs |
| 2 | whole transcriptome, 29 124 contigs / 48.1 Mb | `blastn -task blastn-short -word_size 7 -perc_identity 70` + same product filter | priming anywhere in the assembly |
| 3 | whole transcriptome | positional scan of every start (47.6 M positions), no BLAST threshold | mismatch-tolerant products BLAST would drop |
| 3b | whole transcriptome | threshold probe: how many products at 3 progressively looser rule sets | whether the filters, not the primers, explain the clean result |
| 4 | 6 curated transcripts | exhaustive forward×reverse pairing at four rule sets | mismatch-tolerant cross-amplification inside the family |

## Results

**Check 1** — all six pairs verdict `SPECIFIC`: one product on the intended transcript, zero
off-target. (`PRIMER_BLAST_CHECK.md`, `product_check.tsv`, `blast_hits.tsv`.)

**Check 2** — 52 371 raw primer hits, 11 strong single-primer sites across the assembly. Under the
product filter, one product per pair, all on the expected contig, none on unrelated contigs
(`TRANSCRIPTOME_CHECK.md`, `tx_product_check.tsv`, `tx_blast_hits.tsv`, `tx_priming_sites.tsv`).
SDH1.1 reports zero products, and that is expected rather than alarming:  NODE_10980 is SDH1.2 longest variant 
with intron retention, SDH1.1 - after manual intron removal

**Check 3** — 41 437 priming sites at ≤6 mismatches; 10 candidate products over four rule sets.
Two of them belong to the SDH3 pair and land on family contigs: 309 bp on the SDH2 contig
(NODE_18517) and 300 bp on the SDH4 contig (NODE_21648). This is where the SDH3 risk was first
detected. (`TX_EXHAUSTIVE_SCAN.md`, `exhaustive_products.tsv`, `exhaustive_sites.tsv`.)

**Check 3b** — 5 products under each of three progressively looser variants, consistent with
Check 2 once the thresholds relax; documents how sensitive the conclusion is to the filter
(`TX_SENSITIVITY.md`).

**Check 4** — exhaustive pairing on the six curated transcripts. Its own three rule sets
(strict ≤3 total mismatches with a perfect 8-nt 3' end; medium ≤4 mm; loose ≤6 mm) all report
`cross=0`. Adding a permissive set (≤6 mismatches per primer, ≤12 total, perfect 5-nt 3' run)
yields 10 products of which 4 are cross-products, all from SDH3:

| Pair | Product on | Size | Mismatches F/R | Perfect 3' run F/R |
|---|---|---|---|---|
| CserDQD_SDH3 | CserDQD/SDH4 | 300 bp | 3 / 3 | 7 / 8 |
| CserDQD_SDH3 | CserDQD/SDH2 | 309 bp | 3 / 4 | 7 / 5 |
| CserDQD_SDH3 | CserDQD/SDH4 | 630 bp | 6 / 3 | 6 / 8 |
| CserDQD_SDH3 | CserDQD/SDH2 | 636 bp | 6 / 4 | 6 / 5 |

(`curated_exhaustive_products.tsv`, `TX_CURATED_EXHAUSTIVE.md`, `curated_exhaustive_sites.tsv`.)

## Provenance note, corrected

Earlier versions of the write-up credited this risk to Check 4. That was wrong in two ways, and it
matters for anyone re-reading the evidence:

* **Check 3 found it first.** The 300 bp and 309 bp products are in
  `exhaustive_products.tsv` from the transcriptome scan.
* **Check 4, as originally written, could not have found it.** Its strict/medium/loose rule sets all
  demand a perfect 8-nt 3' anchor. The SDH3 sites carry 7-nt and 6-nt anchors at 3-6 mismatches, so
  they fall in a dead zone between the rules: too many mismatches for the 8-nt-anchor sets, too few
  for the loosest set to matter. Reporting `cross=0` there was correct but not the whole picture.

The permissive rule came from `step5_sdh3_backup/verify_sdh3_backup.py`, which is what produced the
four cross-pairings quoted in the write-up. `curated_exhaustive_scan.py` now carries that rule set
itself, so Check 4 reproduces the finding on its own instead of pointing at another script. Running
the strict/medium/loose rows still gives identical numbers to the archived run; the permissive rows
are the addition.

## Reading the risk

At ordinary stringency the SDH3 pair is specific (Check 1, Check 2, Check 4 strict/medium/loose).
The cross-products need 3-6 mismatches, and both primers anneal with a decent 3' end, so they can
appear if the annealing temperature is dropped. Two of them (300 and 309 bp) sit close to the
intended 300 bp band, so a gel alone may not separate them.

Practical consequences:

1. anneal at 55-58 °C, not lower, for the 6-plex;
2. sequence-purify the 300 bp band rather than trusting the size;
3. if the SDH3 band must be clean by itself, use the backup pair in report 05 — it has zero
   competent foreign sites even under the permissive and very-permissive rules.

## Reproduce

```sh
python3 scripts/step4_specificity/run_primer_blast_check.py    # Check 1 (builds its own BLAST DB)
python3 scripts/step4_specificity/run_transcriptome_check.py   # Check 2   (needs data/external/)
python3 scripts/step4_specificity/tx_exhaustive_scan.py        # Check 3   (needs data/external/)
python3 scripts/step4_specificity/tx_sensitivity_check.py      # Check 3b  (needs data/external/)
python3 scripts/step4_specificity/curated_exhaustive_scan.py   # Check 4
```
