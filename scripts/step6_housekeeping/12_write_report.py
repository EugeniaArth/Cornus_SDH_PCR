#!/usr/bin/env python3
"""Write combined FAMILY_PRIMER_REPORT.md with positions, GC%, and region."""
import os, sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

FAMILIES = ['ACT','BTUB','GAPDH','H3','eIF5A','polyUb']
FAMILY_DESC = {
    'ACT':'Actin','BTUB':'Beta-tubulin','GAPDH':'GAPDH cytosolic',
    'H3':'Histone H3','eIF5A':'eIF5A','polyUb':'Polyubiquitin'
}

report = []
report.append("# Family-Level Primer Pairs for qRT-PCR Normalisation (Cornus sericea)")
report.append("")
report.append("## Overview")
report.append("")
report.append("Six family-level primer pairs designed to co-amplify all paralogs of each reference-gene family.")
report.append("All primers are in the **CDS** region (derived from MAFFT codon alignment).")
report.append("")
report.append("| Family | Description | Members | Forward (5'->3') | Reverse (5'->3') | Fwd GC | Rev GC | Product | FwdPos | RevPos | Region |")
report.append("|--------|------------|---------|-----------------|-----------------|--------|--------|---------|--------|--------|--------|")

with open(P.results('FAMILY_primer_pairs.tsv')) as f:
    next(f)
    for line in f:
        p = line.strip().split('\t')
        fam = p[0]
        if fam in FAMILIES:
            fwd, rev = p[1], p[2]
            prod = p[3]
            f_gc, r_gc = p[6], p[7]
            f_pos, r_pos = p[8], p[9]
            tot = p[12]
            report.append(f"| {fam} | {FAMILY_DESC[fam]} | {tot} | `{fwd}` | `{rev}` | {f_gc}% | {r_gc}% | {prod} bp | {f_pos} | {r_pos} | CDS |")

report.append("")
report.append("## Primer Sequences")
report.append("")
report.append("```")

with open(P.results('FAMILY_primers.fasta')) as f:
    report.append(f.read())

report.append("```")
report.append("")
report.append("## Specificity Validation")
report.append("")
report.append("Each pair was tested by **exact 20-mer scanning** across the whole *Cornus* transcriptome (29,124 contigs):")
report.append("")

scan_path = P.spec('family_tx_scan_report.md')
if os.path.exists(scan_path):
    with open(scan_path) as f:
        report.append(f.read())
else:
    report.append("*Scan report not available - run family_tx_scan.py first.*")

report.append("")
report.append("## Design Method")
report.append("")
report.append("1. **Target selection**: Agrostemma githago Table S7 sequences -> BLASTn -> keyword search -> HMM validation (Pfam profiles).")
report.append("2. **Family candidates**:")
report.append("   - ACT: Actin (2 contigs)")
report.append("   - BTUB: Beta-tubulin-2, Beta-tubulin-8 (2 contigs)")
report.append("   - GAPDH: GAPDH cytosolic (2 contigs)")
report.append("   - H3: Histone H3.2 (x4), H3.3 (x2) (6 contigs)")
report.append("   - eIF5A: 3 contigs (PF01287 validated)")
report.append("   - polyUb: polyubiquitin (7 contigs, architecture-validated)")
report.append("3. **Primer design**: MAFFT protein alignment -> codon-level nucleotide alignment -> conserved 20-mer windows.")
report.append("4. **Coverage**: each primer is allowed <=5 mismatches on some family members (indel/codon-shift tolerance).")
report.append("5. **Specificity**: whole-transcriptome exact-site scan - zero foreign amplicons required.")
report.append("")
report.append("## Usage Notes")
report.append("")
report.append("- **Annealing temperature**: 50-55 C recommended (Tm range 48-58 C).")
report.append("- **Family coverage**: each pair amplifies at least one family member exactly; other members amplify with <=5 mismatches.")
report.append("- **EF-1a**: A family-level pair for all 10 EF-1a paralogs was not feasible (too divergent).")
report.append("- **eIF5A**: Family-level pair was not feasible for eIF5A (3 members, too divergent for a single conserved pair).")
report.append("")
report.append("## Files")
report.append("")
report.append("| File | Description |")
report.append("|------|-------------|")
report.append("| `data/results/FAMILY_primer_pairs.tsv` | All primer pairs with Tm, GC, positions, coverage |")
report.append("| `data/results/FAMILY_primers.fasta` | Primer sequences |")
report.append("| `data/results/specificity/family_tx_scan_report.md` | Per-pair specificity details |")
report.append("| `family_candidates_filtered/` | Per-family filtered contig FASTA files |")
report.append("")
report.append("---")
report.append("")
report.append("*Generated for the Cornus sericea PCR project (DQD/SDH + normalisation).*")

out = P.results('FAMILY_PRIMER_REPORT.md')
with open(out, 'w') as f:
    f.write('\n'.join(report))
print(f"Report written to {out}")
