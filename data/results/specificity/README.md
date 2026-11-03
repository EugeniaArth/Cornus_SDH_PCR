# Specificity evidence — index and provenance

All five specificity checks of the six *CserDQD/SDH* primer pairs, plus the backup SDH3 analysis.
Promoted from `backup/primer_blast_check/` on 2026-11-02 so the evidence sits with the deliverables,
then re-run here to prove the relocated scripts reproduce it (see "Reproducibility status" below).

## Reports

| File | Check | Produced by |
|---|---|---|
| `SPECIFICITY_SUMMARY.md` | Consolidated verdict over all four checks | hand-written consolidation |
| `PRIMER_BLAST_CHECK.md` | 1 — BLAST-like, six curated transcripts | `run_primer_blast_check.py` |
| `TRANSCRIPTOME_CHECK.md` | 2 — BLAST-like, whole transcriptome (29 124 tx) | `run_transcriptome_check.py` |
| `TX_EXHAUSTIVE_SCAN.md` | 3 — exhaustive seed-free scan, whole transcriptome | `tx_exhaustive_scan.py` |
| `TX_SENSITIVITY.md` | 3b — threshold sensitivity of Check 3 | `tx_sensitivity_check.py` |
| `TX_CURATED_EXHAUSTIVE.md` | 4 — exhaustive scan, six curated transcripts | `curated_exhaustive_scan.py` |
| `SDH3_BACKUP_REPORT.md` | backup SDH3 pair search (174 pairs) | `design_sdh3_backup.py` |
| `SDH3_BACKUP_VERIFY.md` | backup pair in full 6-plex context | `verify_sdh3_backup.py` |

## Data tables

| File | Contents |
|---|---|
| `curated_exhaustive_sites.tsv` | every 20-nt window with ≤6 mm per primer × transcript × strand, with 3′-anchor length and mismatch positions |
| `curated_exhaustive_products.tsv` | candidate amplicons under the four rule sets (strict / medium / loose / **permissive**) |
| `exhaustive_sites.tsv`, `exhaustive_products.tsv` | same, transcriptome-wide (47.6 M positions) |
| `tx_priming_sites.tsv`, `tx_selfmap.tsv`, `tx_blast_hits.tsv`, `tx_product_check.tsv` | Check 2 raw + filtered output |
| `blast_hits.tsv`, `product_check.tsv` | Check 1 raw + filtered output |
| `backup_curated_products.tsv` | 6-plex cross-pairing matrix for in-service and backup SDH3 pairs |
| `backup_tx_products.tsv` | backup pair against the whole transcriptome |
| `SDH3_BACKUP_PAIRS.tsv`, `SDH3_backup_candidates.fasta` | the 174 ranked backup pairs |

## Provenance correction (important)

The SDH3 cross-amplification risk (SDH3 pair co-amplifying SDH4 and SDH2) was originally credited to
**Check 4**. That was wrong. As shipped, `curated_exhaustive_scan.py` had only three rule sets:

| rule | mm cap | 3′ anchor required |
|---|---|---|
| strict | ≤3 | ≥5 nt |
| medium | ≤4 | **≥8 nt** |
| loose | ≤6 | **≥8 nt** |

The anchor demand jumps from 5 nt to 8 nt as soon as mismatches exceed 3, so the biologically relevant band
"≤6 mm with a ~5–7 nt anchor" was never tested. The SDH3 forward sites on SDH4 and SDH2 have a **7-nt** anchor —
one base under the 8-nt gate — and were silently dropped. Check 4 therefore printed `cross=0` and its archived
`TX_CURATED_EXHAUSTIVE.md` / `curated_exhaustive_products.tsv` (18 rows, 0 `CROSS-PRODUCT`) contradict the claim.

The risk is real. Its actual sources:

* **Check 3**, `tx_exhaustive_scan.py` → `exhaustive_products.tsv` — first detection. Its product filter is
  `mm_total ≤ 8` with **no** 3′-tail gate (tail lengths are recorded as informational columns only), which is why
  SDH3 → `NODE_21648` (SDH4) 300 bp and SDH3 → `NODE_18517` (SDH2) 309 bp appear there.
* **`verify_sdh3_backup.py`** → `backup_curated_products.tsv`, rows `config = in-service (baseline)`,
  `rule = permissive` — defines the missing rule (`each ≤6 mm, 3′ runs ≥5 nt`, total ≤12) and emits all four rows.

On 2026-11-02 the same rule set was added to `curated_exhaustive_scan.py`, and the script was repathed so it
locates the project root and writes here. Re-running it reproduces the documented finding exactly:

```
permissive (<=6 mm each, <=12 mm total, perfect 5-nt 3' end, 60-1000 bp): products=10 cross=4
  CserDQD_SDH3 -> CserDQD/SDH4 300 bp  mm 3/3  tail 7/8
  CserDQD_SDH3 -> CserDQD/SDH2 309 bp  mm 3/4  tail 7/5
  CserDQD_SDH3 -> CserDQD/SDH4 630 bp  mm 6/3  tail 6/8
  CserDQD_SDH3 -> CserDQD/SDH2 636 bp  mm 6/4  tail 6/5
```

`TX_CURATED_EXHAUSTIVE.md`, `curated_exhaustive_sites.tsv` and `curated_exhaustive_products.tsv` in this folder are
that **re-run**; the 2026-10-29 originals remain unmodified in `backup/primer_blast_check/`.

The two rules are not identical in scope: Check 4 pairs Fwd×Rev **of the same primer pair** on the six curated
transcripts, while `verify_sdh3_backup.py` also tests all 36 cross-pair combinations in the tube. That is why the
verify table lists extra in-service rows (SDH2 382 bp, SDH2 123 bp, SDH4 948 bp) that Check 4 correctly does not
show — they involve a reverse primer belonging to a *different* pair (SDH2 382 bp, SDH4 948 bp), or a site carrying
8 mismatches, outside Check 4's ≤6 mm per-primer cap (SDH2 123 bp).

## Reproducibility status

Every file in this folder except the archival originals was regenerated in place by the scripts in
`scripts/step4_specificity/` and `scripts/step5_sdh3_backup/`, after the project was restructured.
All of those scripts now resolve their inputs through `scripts/paths.py`, build their own BLAST
databases in `data/external/blastdb/`, and write their reports here.

Regression against the archived originals in `backup/primer_blast_check/`: 22 files byte-identical;
`curated_exhaustive_products.tsv` and `TX_CURATED_EXHAUSTIVE.md` carry the added `permissive` rows
described above; `backup_curated_products.tsv` and the two k-mer tables had nondeterministic row
order that has since been fixed in the scripts (identical content, verified order-insensitively).

Two large raw dumps are deliberately **not committed**, because the scripts rewrite them in seconds:
`tx_blast_hits.tsv` (5.0 MB, Check 2) and `exhaustive_sites.tsv` (4.9 MB, Check 3). They are listed
in `.gitignore`. Checks 2, 3 and 3b need the assembly in `data/external/`; Checks 1, 4 and the
six-transcript part of step 5 run on committed data alone.

One naming note: the reports generated by Check 2 and Check 3 print the template as
`Transcriptome/final.clust_transcripts_longest_iso.fasta`, which is where the assembly lived in the
original working directory. In this repository the same file is read from
`data/external/transcriptome.fasta` (and the annotation from `data/external/transcriptome_annotation.gff3`).
The text was left as generated so these files stay byte-identical to the archived originals.
