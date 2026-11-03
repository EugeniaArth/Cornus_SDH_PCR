# Gene-specific PCR primers for six *Cornus sericea* DQD/SDH genes

Six primer pairs that discriminate `CserDQD/SDH1.1`, `SDH1.2`, `SDH2`, `SDH3`, `SDH4` and `SDH5`
in one multiplex, designed from the *C. sericea* transcriptome and checked for specificity
entirely locally.

**Deliverables**

| File | Contents |
|---|---|
| `data/results/FINAL_primers.fasta` | the 12 primers, ready for bulk ordering |
| `data/results/FINAL_primer_pairs.tsv` | sequences, positions, regions, GC, Tm, product sizes |
| `data/results/BACKUP_SDH3_PAIR.txt` | drop-in replacement SDH3 pair with no competent foreign sites |
| `data/results/VALIDATION_auto_report.md` | auto-generated validation of the set |
| `data/results/specificity/` | evidence for all specificity checks, machine-readable |

## Start here

| Report | Read it for |
|---|---|
| [01](reports/01_inputs_and_annotation.md) | which transcripts are the six genes, coordinate verification, why no existing lab primer could be reused |
| [02](reports/02_sequence_survey.md) | where the paralogs actually differ; the SDH1.1 constraint |
| [03](reports/03_primer_design.md) | design rules and the final twelve primers |
| [04](reports/04_specificity_checks.md) | the four-layer specificity screen and the SDH3 cross-amplification finding |
| [05](reports/05_sdh3_backup_pair.md) | the backup SDH3 pair and how it was verified |
| [06](reports/06_lab_plan.md) | ordering, cycling, gel expectations, failure modes |
| [07](reports/07_reproducibility.md) | environment, how to rerun, inventory, open items |

## The result in one table

| Gene | Forward | Reverse | Product |
|---|---|---|---|
| CserDQD/SDH1.1 | GTTCACTGGTGGATTAAGTG | ATAAGCACCTTACTTGAAGC | 120 bp |
| CserDQD/SDH1.2 | GCAAGCAGATTTCAGAGCTG | ACTCAGCAACCCAATAGTGC | 180 bp |
| CserDQD/SDH2 | GCACTTCGATTAGCCATGGA | TGCCACATCAGTGATATCCG | 240 bp |
| CserDQD/SDH3 | TGGAGCCGACATTGTGAAGA | TGACCGACAGGCTTTCCAAT | 300 bp |
| CserDQD/SDH4 | AGAGGAAGCACTAGGAGGTT | AAGAGTCTGGTCCATTTCGG | 360 bp |
| CserDQD/SDH5 | ACTTGCTCAGTCCATAGGTG | ACCTCAGTCACAGGAATTCG | 420 bp |

All 20 nt. Sizes are spaced 60 bp apart so the six amplicons separate on 2 % agarose.
GenBank records: PZ518666 (SDH1.1), PZ518667 (SDH1.2), PZ518668 (SDH2), PZ518669 (SDH3),
PZ518670 (SDH4), PZ518671 (SDH5).

**Two things to know before ordering** (details in reports 04-06):

* The SDH3 pair is specific at normal stringency but retains mismatch-tolerant sites on SDH4
  (300 bp) and SDH2 (309 bp). Anneal at 58 °C, sequence the band, or order the backup pair.
* SDH1.1 has only 38 unique nucleotides in these transcripts, so its pair uses a unique forward
  anchor and a reverse site shared with SDH1.2. A 314 bp SDH1.2 band can appear in the tube.

## Quickstart

```sh
python3 -m pip install numpy            # the only non-stdlib requirement
./scripts/run_all.sh --no-tx            # stages 1-4 on committed data (Checks 1 and 4)
./scripts/run_all.sh                    # adds stage 5 and the transcriptome-wide Checks 2, 3, 3b
```

The transcriptome-wide work reads the assembly from `data/external/`, which is not committed, so
a fresh clone needs those files placed there first — see `data/external/README.md` for names and
symlink commands. Those scripts stop with a clear message when a file is absent instead of
producing a partial result.

Runs 1-4 of the pipeline need only what is committed in `data/raw/`. Each script also works on its
own; they are in `scripts/step1_inputs/` … `scripts/step5_sdh3_backup/`.

Scripts locate their own files through `scripts/paths.py`, so any script can be run from any
working directory, and the folder can be moved or cloned anywhere without editing paths.

## Reproducing the deliverable

`design_final_pairs.py` regenerates `FINAL_primers.fasta` and `FINAL_primer_pairs.tsv` byte for
byte from `data/raw/DQD_SDH_possible_nt.fasta` plus the verified coordinate table. That was
confirmed after the project was restructured: the regenerated primers hash identically to the
committed ones, and 22 of the 26 specificity evidence files do as well
(report 07 lists the exceptions and why).

## Not in this repository

* 48 MB transcriptome FASTA and GFF3 annotation, and the BLAST databases built from them
  (`data/external/`, gitignored).
* Vector-screening work for the NCBI submission — that is a separate task with its own tools and
  outputs, kept in the sibling `vector_work/` directory.
* Photos, InterProScan output and literature spreadsheets, kept in the sibling `reference/`.
* Superseded material: `archive/` inside this folder (redundant duplicates, gitignored) and
  `PCR_SDH/archive/old_layout_20261102/` outside it (the pre-repository layout, with a manifest
  mapping every file to its replacement here). Earlier rounds remain in `PCR_SDH/backup/`.
  Nothing was deleted during the reorganisation.
