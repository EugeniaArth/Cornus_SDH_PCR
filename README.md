# Gene-specific PCR primers for six *Cornus sericea* DQD/SDH genes

Six primer pairs that discriminate `CserDQD/SDH1.1`, `SDH1.2`, `SDH2`, `SDH3`, `SDH4` and `SDH5`
in one multiplex, designed from the *C. sericea* transcriptome and checked for specificity
entirely locally.

NB!!! The transcriptome-wide work reads the assembly needs to be downloaded from
https://drive.google.com/drive/folders/1aDozkY-kDjnAsJzh3PdtklgKbpfJFE5D?usp=sharing
to folder `data/external/`, which is not committed

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


Runs 1-4 of the pipeline need only what is committed in `data/raw/`. Each script also works on its
own; they are in `scripts/step1_inputs/` … `scripts/step5_sdh3_backup/`.

Scripts locate their own files through `scripts/paths.py`, so any script can be run from any
working directory, and the folder can be moved or cloned anywhere without editing paths.

## Reproducing the deliverable

`design_final_pairs.py` regenerates `FINAL_primers.fasta` and `FINAL_primer_pairs.tsv` byte for
byte from `data/raw/Cornus_DQD_SDHs.fasta` plus the verified coordinate table. That was
confirmed after the project was restructured: the regenerated primers hash identically to the
committed ones, and 22 of the 26 specificity evidence files do as well
(report 07 lists the exceptions and why).

## Not in this repository

* 48 MB transcriptome FASTA and GFF3 annotation, and the BLAST databases built from them
  (`data/external/`, gitignored).
  
## Housekeeping genes

Six family‑level primer pairs designed to co‑amplify all paralogs of each reference‑gene family.


- **ACT** — Actin (2 contigs)
- **BTUB** — β-tubulin (2 contigs)
- **GAPDH** — Glyceraldehyde-3-phosphate dehydrogenase, cytosolic (2 contigs)
- **H3** — Histone H3 (6 contigs: H3.2 ×4, H3.3 ×2)
- **eIF5A** — Eukaryotic translation initiation factor 5A (2 contigs)
- **polyUb** — Polyubiquitin (2 contigs, Ubiquitinx6)

| Folder | Strategy | Pairs | Scope | Amplicons |
| `housekeeping/` | **Family‑level** | 6 | total‑family signal | 90–151 bp |

Scripts for housekeeping genes search are in scripts/step6_housekeeping folder 

## Design Method

1. **Target curation**:
     BLASTn of *Agrostemma githago* ortholog sequences (Table S7 (57 contigs across 9 families), 
  reference:  Bielecka, Monika, et al. "Reference gene stability in Agrostemma githago using quantitative real-time PCR." 
   International Journal of Molecular Sciences 27.2 (2026): 889.)
   - Keyword search in Cornus GFF3 annotation using product‑name patterns
   - HMM validation with Pfam profiles 
   PF00125 (Histone) -> H3 (filtered to canonical H3.2/H3.3 by NR product)
   PF01287 (eIF-5a) -> eIF5A (corrected; replaces wrongly-used PF01331/mRNA_cap_enzyme)
   PF00009 (GTP_EFTU) -> EF-1a / TEF1 / organellar EF-Tu (classified by NR product)
   PF00240 (ubiquitin) -> ubiquitin fold; polyubiquitin by domain architecture
   - Product-annotation based filtering to canonical isoforms only

2. **Primer design**:
   - MAFFT protein alignment → codon-level nucleotide alignment
   - Sliding 20-mer window: count mismatches across all family members using real transcript sequences
   - Pick forward + reverse with max coverage, Tm 50‑60 °C, GC 40‑60%, amplicon 70‑250 bp
   - Each primer is allowed ≤2–5 mismatches on some family members (indel/codon‑shift tolerance)
   - All 6 pairs validated: zero exact foreign amplicons across the whole transcriptome (29,124 contigs)

## Summary Table

| Family | BLAST hits | Keyword hits | Intersection | HMM validated | Notes |
|--------|-----------|-------------|-------------|---------------|-------|
| ACT (Actin) | 2 | 155 | 2 | 0 |  |
| BTUB (Beta-tubulin) | 3 | 32 | 2 | 0 |  |
| EF1a (Elongation factor 1-alpha) | 4 | 26 | 4 | 0 | EF1a: 0 |
| GAPDH (GAPDH cytosolic) | 2 | 7 | 2 | 0 |  |
| H3 (Histone H3) | 0 | 11 | 0 | 6 | canonical H3: 6 (45 raw PF00125 hits) |
| eIF5A (eIF5A) | 0 | 3 | 0 | 3 | eIF5A (PF01287): 3 |
| TEF1 (Translation elongation factor 1) | 0 | 12 | 0 | 0 | TEF1: 0 |
| UBQ (Ubiquitin / polyubiquitin) | 0 | 472 | 0 | 54 | polyUb: 12/54 |


## Result

| Gene | Forward | Reverse | Product |
|---|---|---|---|
| ACT | ACTGGTGTGATGGTTGGCAT | GCCAGATCTTCTCCATGTCA | 139 bp |
| BTUB | TGCAGATGAGTGCATGGTTC | AGATCACCAAAGCTTGGGGT | 90 bp |
| GAPDH | ACTGGAGAGGTGGAAGAGCT | ATCAACAGTGGGAACACGGA | 134 bp |
| H3 | GCTCAGGACTTCAAGACTGA | TCCTTGGGCATGATTGTCAC | 146 bp |
| polyUb  | ACAACATCCAGAAGGAGTCG | ATCGTGTCCGAGCTCTCCAC | 121 bp |

- **EF‑1α**: A family‑level pair for all 10 EF‑1α paralogs was not feasible (too divergent).
-**TEF1**: confirmed absent** in *C. sericea* (0 of 46 PF00009 GTP_EFTU hits classify as TEF1 by NR product)


## Key Files

| Path | Description |
|------|-------------|
| `data/housekeeping/family_candidates_filtered/ ACT_filtered.fasta` … `polyUb_filtered.fasta` | Filtered contig sequences per family (input, from curation) |
| `data/housekeeping/family_candidates_filtered/filtered_summary.tsv` | Product annotation for each curated contig |
| `data/housekeeping/primers/results/FAMILY_primer_pairs.tsv` | 6 family-level primer pairs |
| `data/housekeeping/primers/results/FAMILY_primers.fasta` | 12 primer sequences |
| `data/housekeeping/primers/results/FAMILY_PRIMER_REPORT.md` | Full report with design, validation, and recommendations |
| `data/housekeeping/primers/specificity/sensitivuty_report.md` | Per-pair specificity details |


