#!/usr/bin/env python3
"""Write combined FAMILY_PRIMER_REPORT.md"""
import os, sys, math
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

FAMILIES = ['ACT','BTUB','GAPDH','H3','eIF5A','polyUb']
FAMILY_DESC = {
    'ACT':'Actin','BTUB':'Beta-tubulin','GAPDH':'GAPDH cytosolic',
    'H3':'Histone H3','eIF5A':'eIF5A','polyUb':'Polyubiquitin'
}

report = """# Family‑Level Primer Pairs for qRT‑PCR Normalisation (Cornus sericea)

## Overview

Six family‑level primer pairs designed to co‑amplify all paralogs of each reference‑gene family.
Alongside the existing 14 variant‑specific pairs (GAPDH×2, EF‑1α×10, polyUb×2), these provide
options for total‑family signal normalisation.

| Family | Description | Members | Fwd Tm | Rev Tm | Product | Fwd cov | Rev cov |
|--------|------------|---------|--------|--------|---------|---------|---------|
"""

with open(P.results('FAMILY_primer_pairs.tsv')) as f:
    header = next(f)
    for line in f:
        p = line.strip().split('\t')
        fam = p[0]
        if fam in FAMILIES:
            fwd, rev, prod = p[1], p[2], p[3]
            f_tm, r_tm = p[4], p[5]
            f_cov, r_cov, tot = p[8], p[9], p[10]
            report += f"| {fam} | {FAMILY_DESC[fam]} | {tot} | {f_tm}°C | {r_tm}°C | {prod} bp | {f_cov}/{tot} | {r_cov}/{tot} |\\n"

report += """
## Primer Sequences

```
"""

with open(P.results('FAMILY_primers.fasta')) as f:
    report += f.read()

report += """```

## Specificity Validation

Each pair was tested by **exact 20‑mer scanning** across the whole *Cornus* transcriptome (29,124 contigs):

"""
# Read scan results
scan_path = P.spec('sensitivuty_report.md')
if os.path.exists(scan_path):
    with open(scan_path) as f:
        report += f.read()
else:
    report += "*Scan report not available — run family_tx_scan.py first.*\n"

report += """
## Design Method

1. **Target selection**: Agrostemma githago Table S7 sequences → BLASTn → keyword search → HMM validation (Pfam profiles).
2. **Family candidates**:
   - ACT: Actin (2 contigs)
   - BTUB: β‑tubulin‑2, β‑tubulin‑8 (2 contigs)
   - GAPDH: GAPDH cytosolic (2 contigs)
   - H3: Histone H3.2 (×4), H3.3 (×2) (6 contigs)
   - eIF5A: 2 contigs (HMM‑validated)
   - polyUb: polyubiquitin Ubiquitinx6 (2 contigs)
3. **Primer design**: MAFFT protein alignment → codon‑level nucleotide alignment → conserved 20‑mer windows.
4. **Coverance**: each primer is allowed ≤2–5 mismatches on some family members (indel/codon‑shift tolerance).
5. **Specificity**: whole‑transcriptome exact‑site scan — zero foreign amplicons required.

## Usage Notes

- **Annealing temperature**: 50–55 °C recommended (Tm range 48–58 °C).
- **Family coverage**: each pair amplifies at least one family member exactly; other members amplify with ≤2–5 mismatches. Standard PCR stringency (50–55 °C) is sufficient for family‑level co‑amplification.
- **EF‑1α**: A family‑level pair for all 10 EF‑1α paralogs was not feasible (too divergent). The existing 10 variant‑specific pairs are the recommended approach.
- **Comparison with variant‑specific pairs**: The 14 variant‑specific pairs (Norm_PCR) offer paralog‑specific resolution; these family‑level pairs offer total‑family signal. Use either set, not both, for normalisation.

## Files

| File | Description |
|------|-------------|
| `data/results/FAMILY_primer_pairs.tsv` | All 6 primer pairs with Tm, GC, coverage |
| `data/results/FAMILY_primers.fasta` | 12 primer sequences (6 forward + 6 reverse) |
| `data/results/specificity/family_tx_scan_report.md` | Per‑pair specificity details |
| `family_candidates_filtered/` | Per‑family filtered contig FASTA files |
| `family_candidates/` | Unfiltered candidate sets |

---

*Generated 10 Nov 2026 for the Cornus sericea PCR project (DQD/SDH + normalisation).*
"""

out = P.results('FAMILY_PRIMER_REPORT.md')
with open(out, 'w') as f:
    f.write(report)
print(f"Report written to {out}")
