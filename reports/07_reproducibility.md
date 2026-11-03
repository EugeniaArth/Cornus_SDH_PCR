# 07 — Reproducibility, inventory, and what is still open

## Environment

| Tool | Version used |
|---|---|
| Python | 3.14.7 |
| numpy | 2.5.2 |
| BLAST+ (`blastn`, `makeblastdb`) | 2.17.0+ |

Only numpy is required beyond the standard library. Biopython is **not** used; FASTA parsing and
reverse-complement are implemented in `scripts/step3_design/design_final_pairs.py` and
`scripts/step2_survey/unique_20mer_spans.py`. No network access is needed except to fetch the two
missing GenBank records noted below.

## Layout

```
Cornus_PCR_git/
├── README.md                       front page: what this is, how to run it
├── .gitignore
├── reports/                        01-07, one per pipeline stage
├── scripts/
│   ├── paths.py                    single source of truth for every file path
│   ├── run_all.sh                  whole pipeline in order, --no-tx to skip the heavy checks
│   ├── step1_inputs/               annotation parse, coordinate check, prior-primer screen
│   ├── step2_survey/               UTR/CDS split, unique-k-mer and variable-window counts
│   ├── step3_design/               pair generation + validation (the deliverable)
│   ├── step4_specificity/          Checks 1, 2, 3, 3b, 4
│   └── step5_sdh3_backup/          backup SDH3 search and verification
└── data/
    ├── raw/                        curated cDNA, 6 GenBank records, alignment CSVs, primer table
    ├── processed/                  derived tables (UTR/CDS, k-mers, spans, coordinates)
    ├── results/                    deliverables: FINAL_primers.fasta, pairs TSV, backup pair,
    │                               VALIDATION_auto_report.md
    │                               └── specificity/  evidence for Checks 1-4 and the SDH3 work
    └── external/                   NOT committed: transcriptome FASTA/GFF, blastdb/ (16 MB)
```

Scripts resolve paths through `scripts/paths.py`, so they run from any working directory. Nothing
hard-codes an absolute location.

## Full run

```sh
./scripts/run_all.sh --no-tx     # stages 1-4 on committed data only: 14 s measured
./scripts/run_all.sh             # adds stage 5 and Checks 2, 3, 3b: about 2 min total
```

`--no-tx` skips stage 5 as well, because `verify_sdh3_backup.py` rescans the assembly. The candidate
search (`design_sdh3_backup.py`) would run without the assembly, just without its transcriptome
verification step.

Or stage by stage:

```sh
python3 scripts/step1_inputs/parse_genbank.py
python3 scripts/step1_inputs/compare_genbank_coords.py
python3 scripts/step1_inputs/extract_cornus_utrs.py
python3 scripts/step2_survey/analyze_variable_positions.py
python3 scripts/step2_survey/unique_20mer_spans.py
python3 scripts/step3_design/design_final_pairs.py
python3 scripts/step3_design/validate_final_primers.py
python3 scripts/step4_specificity/run_primer_blast_check.py
python3 scripts/step4_specificity/curated_exhaustive_scan.py
python3 scripts/step5_sdh3_backup/design_sdh3_backup.py     # transcriptome part needs data/external/
python3 scripts/step5_sdh3_backup/verify_sdh3_backup.py     # needs data/external/
python3 scripts/step4_specificity/run_transcriptome_check.py
python3 scripts/step4_specificity/tx_exhaustive_scan.py
python3 scripts/step4_specificity/tx_sensitivity_check.py
```

Measured run times on this machine: `run_all.sh --no-tx` 13.6 s; Check 2 ~7 s; Check 3 ~22 s;
Check 3b ~30 s; `design_sdh3_backup.py` ~5 s; `verify_sdh3_backup.py` ~30 s. The one-off heavy cost
is `makeblastdb` on the 48 MB transcriptome, cached afterwards in `data/external/blastdb/`.

## Regression check after the move

The repository was created by relocating a working project, so every artifact was regenerated in
the new layout and compared against the archived original:

| Artifact | Result |
|---|---|
| `FINAL_primers.fasta`, `FINAL_primer_pairs.tsv` | byte-identical |
| `genbank_coordinates_combined.tsv` | byte-identical, once both missing GenBank records were re-fetched |
| `unique_20mer_spans.tsv`-derived processed tables (3 of 5) | byte-identical |
| specificity evidence (Check 1-5, 26 files) | 22 byte-identical, 3 explained below, 3 excluded from git |
| `genbank_features.tsv` | grew from 4 records to 6: the two re-fetched accessions added features |
| `curated_exhaustive_products.tsv` (Check 4) | 10 added rows: the `permissive` rule set. Strict/medium/loose rows unchanged |
| `TX_CURATED_EXHAUSTIVE.md` (Check 4) | differs only by the same added rule set |
| `backup_curated_products.tsv` | identical content; row order of tied rows was nondeterministic, now fixed |
| `cornus_unique_5utr_kmers.tsv`, `cornus_unique_3utr_kmers.tsv` | identical content, row order was set-iteration noise, now sorted |

Two of those differences were genuine defects found by the regression test, and both are fixed in
the scripts rather than patched in the data:

* **Nondeterministic row order.** `extract_cornus_utrs.py` and `verify_sdh3_backup.py` iterated
  Python sets, so output rows shuffled between runs. Sorting was added
  (`sorted(...)` on the k-mer lists, and a full-tuple tiebreak on the product rows). Three
  consecutive runs of `verify_sdh3_backup.py` now produce one hash.
* **Broken path logic.** After the earlier move into `scripts/primer_blast_check/`, the Check 4
  script computed `PARENT = dirname(HERE)` and looked for inputs one level up, which no longer
  existed. It now goes through `paths.py`, and the other 14 scripts were converted the same way.

One script-level wording fix came out of the same exercise: the SDH1.1 unique region was described
as "a ~19 nt block" in `validate_final_primers.py` and "~35 bp" in `design_final_pairs.py`. The
measured value is a 38-nt window (1567-1604) containing 19 unique 20-mer starts. Both docstrings
and the generated report now say that, backed by
`data/processed/unique_20mer_spans.tsv`.

## Files deliberately not in git

| Pattern | Why |
|---|---|
| `data/external/*` | 48 MB transcriptome FASTA, two GFF3 annotation files, 16 MB BLAST databases |
| `data/results/specificity/tx_blast_hits.tsv` | 5.0 MB raw blastn output, rewritten by Check 2 in 7 s |
| `data/results/specificity/exhaustive_sites.tsv` | 4.9 MB site list, rewritten by Check 3 in 22 s |
| `*.tsv.gz` | compressed duplicates of tables that are committed uncompressed |
| `archive/` | local scratch |

## What lives outside git, in `archive/`

`Cornus_PCR_git/archive/` (gitignored) holds redundant copies that were moved rather than deleted:
three `.tsv.gz` duplicates of committed evidence tables, a duplicate of the coordinate table that
belongs in `data/processed/`, and Python bytecode directories. The larger superseded material sits
outside the repository entirely — `PCR_SDH/archive/old_layout_20261102/` (59 MB, with
`MANIFEST.txt` mapping every item to its replacement in this repo), plus `PCR_SDH/backup/` from
earlier rounds.

## Open items

1. **Nothing has been validated in the lab yet.** Every claim in reports 03-06 is in-silico.
   Report 06 lists the wet-lab confirmations, with sequencing of the 120 bp and 300 bp bands as the
   two that matter most.
2. **The transcriptome checks see five of six genes.** SDH1.1 is collapsed into the SDH1.2 contig by
   CD-HIT, so its amplicon cannot be confirmed in that assembly. Confirming SDH1.1 needs genomic
   sequence or a non-clustered assembly.
3. **Primer sites are untested against other genotypes.** A SNP under any of the 12 primers in a
   different *C. sericea* accession would change amplification, and the transcriptome used here is
   one sample.

Closed during reorganisation: the two missing GenBank records (PZ518670, PZ518671) were re-fetched
from NCBI and committed, so the coordinate table is now regenerable and verified byte-identical.

## Correction log

| Date | Change |
|---|---|
| 2026-11-02 | Moved accession map, alignment CSVs and primer table to `archive/`; updated `compare_genbank_coords.py` and `extract_cornus_utrs.py` to read coordinates from GenBank records only; added `download_genbank.py`. Variable‑positions analysis omitted (alignment CSV missing). |
| 2026-11-02 | Whole pipeline rerun after the moves: `run_all.sh` exits 0 in 85 s, no failures, `FINAL_primers.fasta` hashes identically to the deliverable. Fixed a basename-relative path bug in `compare_genbank_coords.py` found by that rerun. |
| 2026-11-02 | PZ518670 and PZ518671 re-fetched from NCBI into `data/raw/`; coordinate table regenerated byte-identically, `genbank_features.tsv` now covers all six genes. |
| 2026-11-02 | Redundant files (3 gz duplicates, duplicate coordinate table, bytecode) moved to `archive/`; superseded pre-repository layout moved to `PCR_SDH/archive/old_layout_20261102/` with a manifest. Nothing deleted. |
| 2026-11-02 | SDH3 cross-amplification risk re-credited to Check 3 (first detection) instead of Check 4, which had no rule able to see it; the permissive rule set was added to Check 4 so the claim is reproducible there. See report 04, "Provenance note". |
| 2026-11-02 | Specificity evidence promoted from `backup/primer_blast_check/` into `data/results/specificity/`, then regenerated in place to prove the relocated scripts reproduce it. |
| 2026-11-02 | Pipeline relocated to this repository, split into step1-step5 folders, path handling centralised in `scripts/paths.py`, two nondeterministic output bugs fixed. |
| 2026-11-02 | Vector-screen work moved out to the sibling `vector_work/` directory; reference material (photos, InterProScan output, literature spreadsheets) moved to the sibling `reference/`. Neither belongs in a primer-design repository. |
