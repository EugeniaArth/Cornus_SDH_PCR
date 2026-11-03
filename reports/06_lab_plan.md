# 06 — Wet-lab plan

What to order, how to cycle, and what a correct result looks like on a gel. Everything here
follows from reports 03-05; nothing new is measured in this step.

## Ordering

* Desalted oligonucleotides are sufficient; no 5' tails, no modified bases.
* Resuspend to 100 µM in TE or nuclease-free water; store working dilutions at 20 µM.
* `data/results/FINAL_primers.fasta` is formatted for bulk upload (one sequence per record,
  name line = gene_direction).
* If using the backup SDH3 pair (report 05), swap the two SDH3 sequences before ordering.

## Six-plex, gel readout

| Component | Amount |
|---|---|
| 10× buffer with 1.5 mM MgCl₂ | 1× |
| Each primer | 0.2 µM |
| dNTPs | 0.2 mM |
| Taq | 1 U |

Cycling: 95 °C 3 min; 30× (95 °C 15 s / **55-58 °C 20 s** / 72 °C 30 s); 72 °C 5 min.
Run 2 % agarose, 60 min at 120 V.

Expected: six bands at **120, 180, 240, 300, 360, 420 bp**. A seventh band near 314 bp is the
known SDH1.2 forward × SDH1.1 reverse combination; it is not a new contaminant.

## RT-qPCR, singleplex

SYBR or EvaGreen 2× master mix, 400 nM each primer, 20 s extension at 72 °C. All six pairs sit
within 6 °C of each other, so 54-56 °C annealing works for screening; run a melt curve every time.

## Temperature choice, and why it is not arbitrary

The in-silico evidence gives one actionable constraint: **the SDH3 pair is the highest-Tm pair of
the set (55.2 / 55.1 °C) and it is also the pair with mismatch-tolerant sites on SDH4 and SDH2.**
Hotter annealing suppresses those sites, and SDH3 can take the heat. The weak link is the opposite
end: SDH1.1 (49.5 / 48.6 °C).

* Multiplex gel: start at 58 °C rather than 55 °C. If the 120 bp SDH1.1 band fades at 58 °C,
  that pair needs re-design (a longer anchor), not a lower temperature.
* SDH3 singleplex: 58-60 °C with melt-curve check. The SDH4 co-product has a different sequence,
  so it should melt apart from the true 300 bp band.
* If SDH3 must be unambiguous by size alone, use the backup pair from report 05.

## Validation sequence

1. Single-plex each pair on a cDNA pool containing all six templates; confirm only the expected band.
2. Multiplex at 55, 56, 57, 58 °C; pick the crispest six-band pattern.
3. Excise and sequence the 120 bp band (prove it is SDH1.1, not the ~1 kb SDH1.2 product) and the
   300 bp band (prove it is SDH3, not SDH4).
4. No-RT control for genomic DNA.

## Known failure modes

| Observation | Cause | Action |
|---|---|---|
| 120 bp band faint or absent | SDH1.1 has the lowest Tm and only a shared reverse site | confirm template present; re-design SDH1.1 with a longer unique anchor if it stays weak |
| extra band ~314 bp | SDH1.2 forward × SDH1.1 reverse on SDH1.2 (perfectly matched pair) | expected; ignore for presence/absence, do not use for quantification |
| extra band ~300/309 bp with SDH3 | mismatch-tolerant priming on SDH4 / SDH2 at soft annealing | raise annealing to 58-60 °C, or switch to the backup SDH3 pair |
| product on no-RT control | amplicons were designed on cDNA; intron coverage was not part of the design | redesign that pair across an exon junction using genomic data |

## Scope of these primers

The design space was the transcriptome, which is the right template for RT-PCR and RT-qPCR on
cDNA. Genomic DNA, introns and intergenic sequence are not covered. The assembly is clustered to
longest isoforms, so allelic variants and low-coverage isoforms are invisible to these checks;
a second *Cornus* genotype may carry SNPs under a primer site.
