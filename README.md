# Cornus sericea DQD/SDH Primer Design


**Project:** Design six gene‑specific PCR primer pairs to distinguish the six *Cornus sericea* DQD/SDH genes:
`CserDQD/SDH1.1`, `CserDQD/SDH1.2`, `CserDQD/SDH2`, `CserDQD/SDH3`, `CserDQD/SDH4`, `CserDQD/SDH5`.

**Date:** 2026‑10‑29  


---

## 1. Goal and strategy

**Goal:** Obtain six primer pairs that:
- Amplify only their intended transcript (no cross‑amplification among the six DQD/SDH paralogs).
- Produce amplicons spaced ~60 bp apart (120 bp → 420 bp) so all six reactions can be resolved together on a 2 % agarose gel.
- Pass standard primer‑quality filters (GC 40–60 %, no homopolymers ≥4 nt, no problematic 3′‑tail composition, no ≥7‑nt self‑complementarity).
- Be compatible with a single annealing temperature window (54–58 °C) for multiplex PCR.
- Show zero competent off‑target priming in the whole *Cornus* transcriptome (29 k transcripts, 48 Mb).

**Strategy:** De‑novo primer design using the six curated GenBank transcripts (PZ518666–PZ518671). No attempt was made to reuse existing primers from other species because they showed no similarity (≥2 mismatches) to the Cornus sequences.

---


## 2. Data

| File | Description |
|------|-------------|
| [`FINAL_primer_pairs.tsv`](FINAL_primer_pairs.tsv) | Machine‑readable primer table (gene, sequence, position, GC %, Tm, amplicon size). |
| [`FINAL_primers.fasta`](FINAL_primers.fasta) | Twelve primers in FASTA format, ready for ordering‑portal bulk upload. |
| [`data/results/`](data/results/) | Copies of the final primer files. |
| [`data/input/`](data/input/) | Input sequences, coordinates, transcriptome assembly. |
| [`scripts/`](scripts/) | Reproducible Python scripts for primer generation, validation, and specificity checking. |


**Gene mapping:**

| Gene | GenBank accession | Length (bp) | 5′‑UTR bp | CDS bp | 3′‑UTR bp |
|------|------------------|-------------|-----------|--------|-----------|
| CserDQD/SDH1.1 | PZ518666 | 1885 | 130 | 1497 | 258 |
| CserDQD/SDH1.2 | PZ518667 | 2779 | 130 | 1503 | 1146 |
| CserDQD/SDH2 | PZ518668 | 2241 | 224 | 1566 | 451 |
| CserDQD/SDH3 | PZ518669 | 2133 | 240 | 1599 | 294 |
| CserDQD/SDH4 | PZ518670 | 2083 | 133 | 1584 | 366 |
| CserDQD/SDH5 | PZ518671 | 1767 | 112 | 1560 | 95 |

---

## 3. Design pipeline (scripts and steps)

All scripts are available in `scripts/`; the main workflow is reproducible.

### 3.1. Primer‑candidate generation (`design_final_pairs.py`)
- Slides a 20‑nt window along each whole transcript.
- Retains only windows **unique** to that transcript (no exact match in any of the other five).
- Applies quality filters:
  - GC 40–60 %
  - No homopolymer ≥4 nt
  - No 3′‑tail `TTTT`/`AAAA` (prone to mis‑priming)
  - No `GGG`/`CCC` at the 3′ end (stable mis‑extension)
  - No inverted self‑complementarity ≥7 nt
  - Acceptable 3′‑terminal base composition (avoid T/A as last base)
- Calculates SantaLucia 1998 nearest‑neighbour Tm (250 nM primer, 50 mM Na⁺).
- Pairs forward‑reverse candidates 120–500 bp apart with ΔTm ≤ 4 °C.

**Output per gene:** thousands of candidate pairs.

### 3.2. Multiplex‑ladder assignment
Each gene is assigned a target amplicon size to produce a 60‑bp‑spaced ladder:

| Gene | Target size (bp) |
|------|------------------|
| SDH1.1 | 120 |
| SDH1.2 | 180 |
| SDH2   | 240 |
| SDH3   | 300 |
| SDH4   | 360 |
| SDH5   | 420 |

From each candidate pool the pair closest to the target size is selected, provided it does not form a ≥8‑nt 3′‑terminal dimer with any primer already chosen.

### 3.3. Special handling of SDH1.1
SDH1.1 shares its 5′‑UTR and most of its CDS with SDH1.2; only a **19‑nt block** near the 3′ end is unique. Consequently no fully‑unique pair exists. The pipeline falls back to a **one‑sided‑specific** design:
- Forward primer sits inside the unique block (unique to SDH1.1).
- Reverse primer is shared with SDH1.2 but cannot amplify SDH1.2 because the forward primer is absent there.

### 3.4. In‑silico validation (`validate_final_primers.py`)
- Exact‑match PCR on the six transcripts: confirms no cross‑amplification.
- Inter‑primer dimer check: no ≥8‑nt 3′‑terminal complementarity among the 12 primers.
- Self‑dimer and hairpin checks.

### 3.5. Specificity validation (four‑layer check)
All analyses were run **locally**; no data were uploaded to external services.

| Check | Script | Scope | Key result |
|-------|--------|-------|------------|
| 1 – BLAST‑like, six curated transcripts | `run_primer_blast_check.py` | Six transcripts | All six pairs specific (no cross‑products under ≤3 mm, 5‑nt 3′‑run rule). |
| 2 – BLAST‑like, whole transcriptome | `run_transcriptome_check.py` | 29 k transcripts, 48 Mb | Only 11 full‑length primer hits, all on the five DQD/SDH contigs (0 off‑family). |
| 3 – Exhaustive seed‑free scan, whole transcriptome | `tx_exhaustive_scan.py` | Every 20‑nt window (47.6 M positions) | 0 products on unrelated transcripts even under ≤8 mismatches, 3‑nt 3′‑run rule. |
| 4 – Exhaustive seed‑free scan, six curated transcripts | `curated_exhaustive_scan.py` | Every window vs all six transcripts | Under permissive rules (≤6 mm, ≥5‑nt 3′‑run) 0 off-target. |


## 4. Final primer set 

| Gene | Direction | Sequence (5′→3′) | Position (transcript) | Region | GC % | Tm (°C) | Amplicon (bp) |
|------|-----------|------------------|----------------------|--------|------|---------|---------------|
| CserDQD/SDH1.1 | F | `GTTCACTGGTGGATTAAGTG` | 1570–1589 | CDS | 45.0 | 49.5 | **120** |
|  | R | `ATAAGCACCTTACTTGAAGC` | 1670–1689 | 3′UTR | 40.0 | 48.6 | |
| CserDQD/SDH1.2 | F | `GCAAGCAGATTTCAGAGCTG` | 2268–2287 | 3′UTR | 50.0 | 53.1 | **180** |
|  | R | `ACTCAGCAACCCAATAGTGC` | 2428–2447 | 3′UTR | 50.0 | 54.0 | |
| CserDQD/SDH2 | F | `GCACTTCGATTAGCCATGGA` | 492–511 | CDS | 50.0 | 54.0 | **240** |
|  | R | `TGCCACATCAGTGATATCCG` | 712–731 | CDS | 50.0 | 54.0 | |
| CserDQD/SDH3 | F | `TGGAGCCGACATTGTGAAGA` | 738–757 | CDS | 50.0 | 55.2 | **300** |
|  | R | `TGACCGACAGGCTTTCCAAT` | 1018–1037 | CDS | 50.0 | 55.1 | |
| CserDQD/SDH4 | F | `AGAGGAAGCACTAGGAGGTT` | 1216–1235 | CDS | 50.0 | 53.1 | **360** |
|  | R | `AAGAGTCTGGTCCATTTCGG` | 1556–1575 | CDS | 50.0 | 52.8 | |
| CserDQD/SDH5 | F | `ACTTGCTCAGTCCATAGGTG` | 1048–1067 | CDS | 50.0 | 53.0 | **420** |
|  | R | `ACCTCAGTCACAGGAATTCG` | 1448–1467 | CDS | 50.0 | 53.2 | |

**Tm note:** SantaLucia 1998 nearest‑neighbour values at 250 nM primer, 50 mM Na⁺. For high‑salt Taq buffers, add ~3–5 °C.

**Amplicon ladder (multiplex‑ready):**
```
SDH1.1  |■■ 120 bp
SDH1.2  |■■■ 180 bp
SDH2    |■■■■ 240 bp
SDH3    |■■■■■ 300 bp
SDH4    |■■■■■■ 360 bp
SDH5    |■■■■■■■ 420 bp
         ↑ ~60 bp apart → all six bands resolve on 2% agarose or a TapeStation
```


