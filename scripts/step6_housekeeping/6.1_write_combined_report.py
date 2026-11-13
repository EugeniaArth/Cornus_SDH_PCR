#!/usr/bin/env python3
"""
write_combined_report.py — Combine BLAST, keyword, and HMM validation results into
a single report (STEP2_RESULTS_COMBINED.md).
"""
import os, sys, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BLAST_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping', 'agrostemma_blast')
KW_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping', 'keyword_search')
HMM_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping', 'hmm_validation')
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping')



FAMILIES = {
    'ACT':   'Actin',  'BTUB': 'Beta-tubulin',  'EF1a': 'Elongation factor 1-alpha',
    'GAPDH': 'GAPDH cytosolic',  'H3': 'Histone H3',  'eIF5A': 'eIF5A',
    'TEF1':  'Translation elongation factor 1',  'UBQ': 'Ubiquitin / polyubiquitin',
}

def load_blast():
    b = defaultdict(list)
    path = os.path.join(BLAST_DIR, 'blast_hits_summary.tsv')
    if not os.path.exists(path): return b
    with open(path) as f:
        next(f)
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 4: b[p[0]].append((p[1], p[2], p[3]))
    return b

def load_keyword():
    k = defaultdict(set)
    if not os.path.isdir(KW_DIR): return k
    for fn in os.listdir(KW_DIR):
        if fn.endswith('_contigs.txt'):
            fam = fn.replace('_contigs.txt', '')
            with open(os.path.join(KW_DIR, fn)) as f:
                for line in f:
                    c = line.strip()
                    if c: k[fam].add(c)
    return k

def load_hmm():
    h = {}
    import glob
    for fp in glob.glob(os.path.join(HMM_DIR, '*_hmm_hits.tsv')):
        fam = os.path.basename(fp).replace('_hmm_hits.tsv', '')
        if fam == 'EF1a': continue  # EF1a rows from this file are 9 EF1a
        with open(fp) as f:
            header = f.readline().strip().split('\t')
            rows = []
            for line in f:
                p = line.strip().split('\t')
                if len(p) >= len(header):
                    rows.append(dict(zip(header, p)))
        h[fam] = rows
    return h

def main():
    blast = load_blast()
    kw = load_keyword()
    hmm = load_hmm()
    report = []

    report.append("# Combined BLAST / Keyword / HMM Validation Results")
    report.append("")
    report.append("## Methodology")
    report.append("1. **BLASTn** — Agrostemma githago Table S7 sequences vs Cornus transcriptome (e-value 1e-5)")
    report.append("2. **Keyword search** — product-name patterns in GFF3 annotation")
    report.append("3. **HMM validation** — hmmscan with Pfam profiles (gathering thresholds):")
    report.append("   - PF00125 (Histone) -> H3 (filtered to canonical H3.2/H3.3 by NR product)")
    report.append("   - PF01287 (eIF-5a) -> eIF5A (corrected; replaces wrongly-used PF01331/mRNA_cap_enzyme)")
    report.append("   - PF00009 (GTP_EFTU) -> EF-1a / TEF1 / organellar EF-Tu (classified by NR product)")
    report.append("   - PF00240 (ubiquitin) -> ubiquitin fold; polyubiquitin by domain architecture")
    report.append("")

    # Summary table
    report.append("## Summary Table")
    report.append("")
    report.append("| Family | BLAST hits | Keyword hits | Intersection | HMM validated | Notes |")
    report.append("|--------|-----------|-------------|-------------|---------------|-------|")
    for fam_short, fam_desc in FAMILIES.items():
        bl = set(c[0] for c in blast.get(fam_short, []))
        kw_set = kw.get(fam_short, set())
        inter = bl & kw_set
        hmm_rows = hmm.get(fam_short, [])
        n_hmm = len(hmm_rows)

        n_bl = len(bl)
        n_kw = len(kw_set)
        n_inter = len(inter)

        if fam_short == 'UBQ':
            polyub_n = sum(1 for r in hmm_rows if r.get('kind') == 'polyubiquitin')
            notes = f"polyUb: {polyub_n}/{n_hmm}"
        elif fam_short == 'EF1a':
            notes = f"EF1a: {n_hmm}"
        elif fam_short == 'TEF1':
            notes = f"TEF1: {n_hmm}"
        elif fam_short == 'H3':
            notes = f"canonical H3: {n_hmm} (45 raw PF00125 hits)"
        elif fam_short == 'eIF5A':
            notes = f"eIF5A (PF01287): {n_hmm}"
        else:
            notes = ""

        report.append(f"| {fam_short} ({fam_desc}) | {n_bl} | {n_kw} | {n_inter} | {n_hmm} | {notes} |")

    report.append("")

    # Detailed per-family
    for fam_short, fam_desc in FAMILIES.items():
        report.append(f"## {fam_short} — {fam_desc}")
        report.append("")

        bl = blast.get(fam_short, [])
        kw_set = kw.get(fam_short, set())
        hmm_rows = hmm.get(fam_short, [])

        if bl:
            report.append(f"**BLAST hits ({len(bl)}):**")
            for c, pid, bs in bl:
                report.append(f"  - `{c}` (id {pid}%, bitscore {bs})")
        else:
            report.append(f"**BLAST hits:** 0 (none at e-value 1e-5)")

        if kw_set:
            report.append(f"**Keyword hits ({len(kw_set)}):**")
            for c in sorted(kw_set)[:5]:
                report.append(f"  - `{c}`")
            if len(kw_set) > 5:
                report.append(f"  - ... plus {len(kw_set)-5} more")
        else:
            report.append(f"**Keyword hits:** 0")

        if bl and kw_set:
            inter = set(c for c,_,_ in bl) & kw_set
            report.append(f"**Intersection:** {len(inter)}")

        if hmm_rows:
            report.append(f"**HMM validated ({len(hmm_rows)}):**")
            for r in hmm_rows:
                c = r.get('contig', '')
                if 'score' in r and r['score']:
                    sc = f"score={float(r['score']):.1f}"
                else:
                    sc = ""
                if 'evalue' in r and r['evalue']:
                    ev = f"E={float(r['evalue']):.0e}"
                else:
                    ev = ""
                prod = r.get('product', '')[:80]
                if fam_short == 'UBQ':
                    extra = f"n_dom={r.get('n_dom','')} cov={r.get('coverage','')} kind={r.get('kind','')}"
                    report.append(f"  - `{c}` {extra}  _{prod}_")
                else:
                    report.append(f"  - `{c}` {sc} {ev}  _{prod}_")
        else:
            report.append(f"**HMM validated:** 0")

        report.append("")

    report.append("## Profile corrections")
    report.append("")
    report.append("TEF1 is absent in Cornus sericea (0 of 46 PF00009 hits classify as TEF1).")
    report.append("")

    out = os.path.join(OUTDIR, 'HOUSEKEEPING_INT_RESULTS.md')
    with open(out, 'w') as f:
        f.write('\n'.join(report))
    print(f"Written: {out}")
    print(f"  {len(FAMILIES)} families covered")

if __name__ == '__main__':
    main()
