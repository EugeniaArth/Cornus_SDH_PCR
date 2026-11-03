# 05 — Backup SDH3 pair with zero competent foreign sites

**Step scripts:** `scripts/step5_sdh3_backup/design_sdh3_backup.py`,
`scripts/step5_sdh3_backup/verify_sdh3_backup.py`
**Outputs:** `data/results/BACKUP_SDH3_PAIR.txt`,
`data/results/specificity/SDH3_BACKUP_PAIRS.tsv`, `SDH3_backup_candidates.fasta`,
`SDH3_BACKUP_REPORT.md`, `SDH3_BACKUP_VERIFY.md`, `backup_curated_products.tsv`,
`backup_tx_products.tsv`

## What had to be fixed

The in-service SDH3 pair (`TGGAGCCGACATTGTGAAGA` / `TGACCGACAGGCTTTCCAAT`) makes four
mismatch-tolerant cross-products inside the family (report 04). The goal was a replacement SDH3
pair that keeps the multiplex intact and leaves no competent priming site on SDH2 or SDH4.

## Search criteria

`design_sdh3_backup.py` required, for every candidate:

* 20 nt, exactly unique to SDH3;
* GC 40-60 %, no homopolymer run >4 nt, acceptable 3' tail;
* no self-complementary hairpin of ≥7 nt;
* product exactly 300 bp, so the ladder position does not move;
* ΔTm ≤ 2.0 °C between the two primers;
* no ≥8-nt 3'-terminal complementarity with the other ten primers in the multiplex;
* **zero competent foreign sites** — no site anywhere in the six transcripts with ≤6 mismatches
  and a perfect 3' run ≥5 nt. The last condition is the point of the exercise: a single primer
  that can initiate on a paralog is enough to spoil the assay, so it is screened at primer level
  (tier A) and at pair level (tier B).

**Result: 174 candidates, 71 in tier A (primer-level clean) and 103 in tier B (pair-level clean).**
Full ranking in `SDH3_BACKUP_PAIRS.tsv`, sequences in `SDH3_backup_candidates.fasta`.

## Two candidates, and why the second one was picked

| | top-ranked by the script | **recommended** |
|---|---|---|
| Forward | TCATTAGCAATGGCCTCAGG | **CACCCAGAAGAAGGGATGGT** |
| Reverse | ATCACCTTCATACTGACCGC | **GCAGACTCTCTCTTATGCGT** |
| Position on SDH3 | 232-251 (5'UTR) / 512-531 (CDS) | 1552-1571 (CDS) / 1832-1851 (CDS) |
| GC | 50 % / 50 % | 55 % / 50 % |
| Tm (NN) | 53.4 °C / 53.6 °C, ΔTm 0.2 | 54.7 °C / 53.4 °C, ΔTm 1.3 |
| Product | 300 bp | 300 bp |
| Competent foreign sites (≤6 mm, ≥5-nt run) | 0 | 0 |
| Very-permissive foreign sites (≤8 mm, ≥3-nt run) | **2** | **0** |

The score function favours ΔTm and Tm uniformity, which puts the 5'UTR/CDS pair first. The
recommendation goes one step further: `verify_sdh3_backup.py` re-audits candidates at
`≤8 mismatches with a ≥3-nt 3' run`, and only the second candidate holds at zero. It also sits
entirely in the CDS, so there is no UTR-junction ambiguity about where it anneals. That is
`#2 zero foreign sites even at ≤8 mm/3-nt run` in the verification report.

## Verification of the swap

`verify_sdh3_backup.py` substitutes the candidate into the six-plex and rescans all 36 forward×reverse
combinations, then repeats the scan on the whole transcriptome.

Six curated transcripts — number of non-designed products in the tube:

| Rule | in-service set | with backup SDH3 |
|---|---|---|
| strict (≤3 mm total, ≥5-nt run) | 2 | 1 |
| competent (≤6 mm total, ≥5-nt run) | 6 | 1 |
| long anchor (≤8 mm total, ≥8-nt run) | 3 | 2 |
| permissive (≤6 mm each, ≥5-nt run) | 16 | 3 |
| very permissive (≤8 mm each, ≥3-nt run) | 64 | 33 |

Whole transcriptome with the backup pair: 5 designed amplicons detectable (SDH1.1 has no contig of
its own, see report 04), zero on-family and zero off-family products at the competent rule; the
only remaining signal at the loosest rule is a 878 bp SDH1.1/SDH1.2 combination that already exists
in the current set.

**One residual, unchanged by this swap:** F(SDH1.2) + R(SDH1.1) gives a 314 bp product on SDH1.2
with both primers perfectly matched. It comes from the SDH1.1 shared-reverse design, not from
SDH3, and it is present in the baseline as well. It is the reason SDH1.1 and SDH1.2 should not be
read as a quantitative pair without checking.

## Using it

Replace `CserDQD_SDH3_Fwd` and `CserDQD_SDH3_Rev` in `FINAL_primers.fasta` with the two sequences
above. Same 300 bp band, same gel pattern, compatible with the 55-58 °C annealing window. Ordering
information and the full rationale are in `data/results/BACKUP_SDH3_PAIR.txt`.

## Reproduce

```sh
python3 scripts/step5_sdh3_backup/design_sdh3_backup.py    # needs data/external/ for the txome scan
python3 scripts/step5_sdh3_backup/verify_sdh3_backup.py
```
