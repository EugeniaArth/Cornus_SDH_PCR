# Final Cornus DQD/SDH Primer Set (de novo)

Six gene-specific PCR primer pairs designed against the six *Cornus sericea*
DQD/SDH transcripts (PZ518666–PZ518671). All primers are 20 nt (inside the
usual 18–24 nt window). Amplicons are spread across 120 → 420 bp in ~60 bp
steps so all six reactions can be resolved together on a 2 % agarose gel or a
fragment analyser.

Sources:
- Transcripts: `DQD_SDH_possible_nt.fasta`
- Feature coordinates (verified against GenBank): `genbank_coordinates_combined.tsv`
- Design script: `design_final_pairs.py`
- Validation script: `validate_final_primers.py`
- Machine-readable table: `FINAL_primer_pairs.tsv`
- Primer FASTA: `FINAL_primers.fasta`

---

## 1. Final primer list

| Gene | Acc. | Dir | Sequence (5'→3') | Len | Position (transcript) | Region | GC % | Tm (NN, °C) | Amplicon (bp) |
|---|---|---|---|---|---|---|---|---|---|
| CserDQD/SDH1.1 | PZ518666 | F | `GTTCACTGGTGGATTAAGTG` | 20 | 1570–1589 | CDS | 45.0 | **49.5** | **120** |
|  |  | R | `ATAAGCACCTTACTTGAAGC` | 20 | 1670–1689 | 3′UTR | 40.0 | **48.6** |  |
| CserDQD/SDH1.2 | PZ518667 | F | `GCAAGCAGATTTCAGAGCTG` | 20 | 2268–2287 | 3′UTR | 50.0 | **53.1** | **180** |
|  |  | R | `ACTCAGCAACCCAATAGTGC` | 20 | 2428–2447 | 3′UTR | 50.0 | **54.0** |  |
| CserDQD/SDH2 | PZ518668 | F | `GCACTTCGATTAGCCATGGA` | 20 | 492–511 | CDS | 50.0 | **54.0** | **240** |
|  |  | R | `TGCCACATCAGTGATATCCG` | 20 | 712–731 | CDS | 50.0 | **54.0** |  |
| CserDQD/SDH3 | PZ518669 | F | `TGGAGCCGACATTGTGAAGA` | 20 | 738–757 | CDS | 50.0 | **55.2** | **300** |
|  |  | R | `TGACCGACAGGCTTTCCAAT` | 20 | 1018–1037 | CDS | 50.0 | **55.1** |  |
| CserDQD/SDH4 | PZ518670 | F | `AGAGGAAGCACTAGGAGGTT` | 20 | 1216–1235 | CDS | 50.0 | **53.1** | **360** |
|  |  | R | `AAGAGTCTGGTCCATTTCGG` | 20 | 1556–1575 | CDS | 50.0 | **52.8** |  |
| CserDQD/SDH5 | PZ518671 | F | `ACTTGCTCAGTCCATAGGTG` | 20 | 1048–1067 | CDS | 50.0 | **53.0** | **420** |
|  |  | R | `ACCTCAGTCACAGGAATTCG` | 20 | 1448–1467 | CDS | 50.0 | **53.2** |  |

*Tm is the SantaLucia 1998 unified nearest-neighbour value at 250 nM primer,
50 mM Na⁺. If you order from an IDT/Thermo-style portal their NNSalt/Tm
tool will return something very similar. If you use a Taq with high-salt
buffer, add ~3–5 °C.*

### Amplicon spacing (multiplex-ready)

```
SDH1.1  |■■ 120 bp
SDH1.2  |■■■ 180 bp
SDH2    |■■■■ 240 bp
SDH3    |■■■■■ 300 bp
SDH4    |■■■■■■ 360 bp
SDH5    |■■■■■■■ 420 bp
         ↑ ~60 bp apart → all six bands resolve on 2% agarose or a TapeStation
```

---

## 2. In-silico specificity check

Exact-match PCR (both primers must match perfectly) against all six Cornus
transcripts:

| Pair | Target transcript | Non-target amplification |
|---|---|---|
| SDH1.1 F/R | PZ518666 – 120 bp | none |
| SDH1.2 F/R | PZ518667 – 180 bp | none |
| SDH2 F/R | PZ518668 – 240 bp | none |
| SDH3 F/R | PZ518669 – 300 bp | none |
| SDH4 F/R | PZ518670 – 360 bp | none |
| SDH5 F/R | PZ518671 – 420 bp | none |

- ✓ No pair amplifies a non-target Cornus DQD/SDH transcript.
- ✓ No 8-nt 3′-complementarity among the 12 primers (no cross-primer dimer).
- ✓ No primer has a self-dimer at ≥8 nt 3′ overlap.
- ✓ No homopolymer run ≥4 nt.
- ✓ All primers 20 nt, GC 40–50 %, Tm 48–55 °C (matched within each pair).

**Caveat** — this check is against the six annotated transcripts only.
Off-target hits in the wider *Cornus* transcriptome/genome still need to be
tested with BLAST (or `bowtie2-build` on the transcriptome) before ordering.

---

## 3. Design logic in plain terms

For every transcript the pipeline:

1. Slides a 20-nt window along the whole transcript and flags every window
   that does **not** occur in any of the other five transcripts.
2. Keeps only windows that pass basic primer quality:
   GC 40–60 %, no homopolymer ≥ 4 nt, no 3′ TTTT/AAAA, no GGG/CCC 3′ run,
   no ≥7-nt inverted self-complementarity.
3. Tries every forward/backward pair of *unique* windows 120–500 bp apart
   with matched NN-Tm (ΔTm ≤ 4 °C). That gives a per-gene candidate pool of
   a few thousand pairs.
4. Assigns each gene a target amplicon size (120 / 180 / 240 / 300 / 360 /
   420 bp), picks the highest-quality candidate nearest to that target, and
   rejects candidates that form an inter-primer dimer with anything already
   chosen.

## 4. Special case: SDH1.1 vs SDH1.2

`CserDQD/SDH1.1` (PZ518666) and `CserDQD/SDH1.2` (PZ518667) are **two splice
isoforms of the same locus**:

- Identical 5′UTR.
- Nearly identical CDS.
- Divergent 3′ ends.

Across the *whole* 1885 nt of SDH1.1 there are **only 19 unique 20-mers**, and
they all fall inside a ~19 nt block (transcript positions 1566–1584, right
around the CDS 3′ end). That block is too small to hold two 20-mers 120+ bp
apart, so no fully-both-unique pair exists for SDH1.1.

For SDH1.1 the pipeline falls back to a **one-sided-specific design**:

- Forward primer `GTTCACTGGTGGATTAAGTG` sits inside the 19-nt unique block
  → only SDH1.1 can prime extension in this reaction.
- Reverse primer `ATAAGCACCTTACTTGAAGC` sits in the 3′UTR and is *not*
  unique — it also appears in SDH1.2, but SDH1.2 has no matching forward
  primer, so it is not amplified.

The exact-match PCR above confirms: only PZ518666 gives a product.
If you want a *fully both-unique* pair for SDH1.1, we would have to widen
the unique block, which means switching to a genomic / exon-junction
strategy rather than transcript-based design.

---

## 5. Why only ONE primer pair per gene, not many alternatives?

This is what the earlier runs (`design_multiplex_primers.py`,
`flexible_multiplex_design.py`, `final_primer_design.py`) kept hitting.

Two different questions get conflated here:

### 5a. Why did the earlier scripts find only one *primer* (single oligo) per
gene 2–5, and never a working pair?

Because they had **three filters stacked together**, and any one of them was
too strict:

1. **Search region**: they looked for primers only inside the last ~60 bp of
   the 5′UTR and the first ~60 bp of the 3′UTR. Most paralogs share long
   patches in those regions, so almost no 20-mer survived the uniqueness
   filter, and the few that did often failed GC / homopolymer / 3′-end rules.
2. **Hairpin test**: the earlier `self_complementarity(..., threshold=4)`
   rejected a candidate if **any** 4-base inverted match appeared. In a
   random 20-mer that happens ~95 % of the time, so almost everything was
   rejected — including the only forward primer that would have paired with
   the found reverse primer.
3. **Product length**: they demanded 200–1500 bp from UTR-to-UTR pairs, but
   the CDS alone is 1500–1600 bp, so almost every surviving pair exceeded
   the length ceiling and got discarded.

Result: the pipeline kept printing things like
`CserDQD/SDH3: F = CCAAGAAACACAACCAGACC` — a single usable primer but no
partner that survived all three filters.

The current `design_final_pairs.py` fixes this by:

- widening the search from "ends of UTRs only" to the **whole transcript**;
- replacing the 4-nt hairpin rule with a realistic **7-nt inverted match**;
- allowing **120–500 bp** amplicons (qPCR / short multiplex products);
- separating the uniqueness constraint from the length constraint, so we
  never again reject a genuinely unique 20-mer just because its natural
  partner was too far.

### 5b. Why do we present only ONE pair per gene in the final set?

Not because there is only one option — there are thousands of passing
candidates for genes SDH2–5 (see the `candidate pairs` counts printed by
`design_final_pairs.py`: ~87 k to ~120 k per gene).

The reason is that for a **sixplex** PCR the six amplicons have to be spread
over the separable size range (here 120 / 180 / 240 / 300 / 360 / 420 bp).
Once the pipeline locks each gene onto one of those six size slots and
rejects any pair that forms a 3′ dimer with a pair already chosen, each gene
effectively has one best representative returned to the final table. If you
want backup / alternative pairs, just ask and we will dump, per gene, the
top 3–5 candidates that still pass all filters but are in neighbouring size
bins.

---

## 6. Suggested PCR set-up

**Conventional sixplex (gel readout)**

- 10× buffer with 1.5 mM MgCl₂, 0.2 mM each primer, 0.2 mM dNTPs, 1 U Taq.
- 95 °C 3 min; 30 cycles of 95 °C 15 s / **55 °C 20 s** / 72 °C 30 s;
  72 °C 5 min.
- 2 % agarose, 60 min at 120 V → six bands at 120 / 180 / 240 / 300 / 360 /
  420 bp.

**RT-qPCR (singleplex per gene)**

- SYBR / EvaGreen 2× master mix, 400 nM each primer, 72 °C 20 s extension.
- Use the same annealing window 54–56 °C; all six pairs are within 6 °C of
  each other so a single cycling programme works for all.

**Validation before ordering**

1. BLAST each primer pair against the full *Cornus* transcriptome
   (`annotated_final_nr.gff3` + matching cDNA).
2. Order as PAGE-pure or desalted, resuspend to 100 µM.
3. Test each singleplex on a mixed cDNA pool (all six templates together)
   and confirm only the expected band appears.

---

## 7. Files created

- `FINAL_primer_pairs.tsv` — the primer table above (tab-separated).
- `FINAL_primers.fasta` — twelve primers, one per line, ready for
  ordering-portal bulk upload.
- `FINAL_primers_report.md` — this report.
- `design_final_pairs.py` — reproducible generator.
- `validate_final_primers.py` — in-silico PCR and dimer checker.
