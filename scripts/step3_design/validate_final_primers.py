#!/usr/bin/env python3
"""
Cross-validate the FINAL primer pairs against all six Cornus transcripts.
Simulate a simple PCR:
  - For each pair, scan every transcript on both strands for exact-match
    forward binding and downstream reverse binding within 100-600 bp.
  - Report any amplicons found and the gene they came from.
  - Flag non-specific amplification.
Also flag cross-pair dimer potential with a more realistic k>=8 3' overlap test.
"""
import csv, os
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

def rc(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp[c] for c in reversed(s))

def load_fasta(path):
    seqs, sid, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sid: seqs[sid] = ''.join(buf)
                sid = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
        if sid: seqs[sid] = ''.join(buf)
    return seqs

def load_primers(path):
    primers = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            primers.append(row)
    return primers

def gc_content(s):
    return (s.count('G') + s.count('C')) / len(s) * 100

def tm_simple(s):
    return 64.9 + 41 * (s.count('G') + s.count('C') - 16.4) / len(s)

import math
_NN = {
    'AA': (-7.6, -21.7), 'TT': (-7.6, -21.7),
    'AT': (-7.2, -20.4),
    'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.1),
    'GC': (-9.8, -24.4),
    'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
}

def tm_nn(s, dna_conc=250e-9, na_conc=50e-3):
    dH = 0.0
    dS = 0.0
    for i in range(len(s) - 1):
        step = s[i:i+2]
        if step not in _NN:
            return tm_simple(s)
        h, ss = _NN[step]
        dH += h
        dS += ss
    for end_base in (s[0], s[-1]):
        if end_base in 'GC':
            dH += 0.1
            dS += -2.8
        elif end_base in 'AT':
            dH += 2.3
            dS += 4.1
        else:
            return tm_simple(s)
    dS += 0.368 * (len(s) - 1) * math.log(na_conc)
    R = 1.987
    return dH * 1000.0 / (dS + R * math.log(dna_conc / 4.0)) - 273.15

def amplicons_for_pair(fwd, rev, seq, max_len=600, min_len=80):
    """Return list of (start, end, len) of products from a single transcript."""
    hits = []
    rc_rev = rc(rev)
    n = len(seq)
    # exact search forward
    i = 0
    while True:
        i = seq.find(fwd, i)
        if i == -1: break
        # find rc_rev downstream within max_len
        j = seq.find(rc_rev, i)
        while j != -1 and j < i + max_len:
            amp_len = (j + len(rc_rev)) - i
            if amp_len >= min_len:
                hits.append((i+1, j+len(rc_rev), amp_len))
            j = seq.find(rc_rev, j+1)
        # reverse-strand product (fwd rc + rev): skip for standard orientation
        i += 1
    # also check reverse orientation: rev + rc(fwd) on the same strand
    j = 0
    while True:
        j = seq.find(rev, j)
        if j == -1: break
        k = seq.find(rc(fwd), j)
        while k != -1 and k < j + max_len:
            amp_len = (k + len(rc(fwd))) - j
            if amp_len >= min_len:
                hits.append((j+1, k+len(rc(fwd)), amp_len, 'rev_orientation'))
            k = seq.find(rc(fwd), k+1)
        j += 1
    return hits

def dimer_3p(a, b, k=8):
    """Return True if 3' end of a matches last k bases of rc(b) (primer dimer)."""
    if len(a) < k or len(b) < k:
        return False
    return a[-k:] == rc(b)[-k:] or b[-k:] == rc(a)[-k:]

def self_dimer(a, k=8):
    if len(a) < 2*k:
        return False
    r = rc(a)
    return any(a[i:i+k] in r for i in range(len(a)-k+1))

def main():
    seqs = load_fasta(P.require(P.raw('Cornus_DQD_SDHs.fasta'), 'curated cDNA, report 01'))
    primers = load_primers(P.require(P.results('FINAL_primer_pairs.tsv'), 'from step3_design/design_final_pairs.py'))

    # Group forward+reverse per gene
    pairs = {}
    for row in primers:
        g = row.get('Gene') or ''
        if not g:
            # continuation row (reverse line)
            if '_current' in locals() and _current:
                pairs[_current]['rev'] = row['Sequence (5\'→3\')']
            continue
        _current = g
        pairs[g] = {
            'acc': row['Accession'],
            'fwd': row['Sequence (5\'→3\')'],
            'mode': row['Mode'],
        }
        if row['Dir'] == 'R':
            pairs[g]['rev'] = row['Sequence (5\'→3\')']

    # The TSV stores R on next line with empty Gene, so we re-parse:
    rows = list(csv.DictReader(open(P.results('FINAL_primer_pairs.tsv')), delimiter='\t'))
    pairs = {}
    current_gene = None
    for r in rows:
        if r['Gene']:
            current_gene = r['Gene']
            pairs[current_gene] = {'acc': r['Accession'], 'mode': r['Mode']}
        if r['Dir'] == 'F':
            pairs[current_gene]['fwd'] = r["Sequence (5'→3')"]
        elif r['Dir'] == 'R':
            pairs[current_gene]['rev'] = r["Sequence (5'→3')"]

    print(f"Primer pairs loaded for genes: {list(pairs.keys())}\n")

    print("="*90)
    print("IN-SILICO PCR: exact-match amplicons on each Cornus transcript")
    print("="*90)

    non_specific = []
    for g, p in pairs.items():
        print(f"\n{g}  ({p['acc']})")
        print(f"  F: {p['fwd']}")
        print(f"  R: {p['rev']}")
        for gene, seq in seqs.items():
            hits = amplicons_for_pair(p['fwd'], p['rev'], seq)
            if hits:
                tag = "✓ TARGET" if gene == g else "⚠ CROSS"
                for h in hits:
                    print(f"    {tag}  {gene}: {h[2]} bp ({h[0]}-{h[1]})")
                if gene != g:
                    non_specific.append((g, gene, hits[0][2]))

    print("\n" + "="*90)
    print("Cross-amplification summary")
    print("="*90)
    if non_specific:
        for src, target, size in non_specific:
            print(f"  ⚠ {src} pair could also amplify {target} ({size} bp)")
    else:
        print("  ✓ No pair amplifies a non-target Cornus DQD/SDH transcript.")

    # Primer dimers
    print("\n" + "="*90)
    print("Cross-primer 3'-complementarity check (last 8 nt)")
    print("="*90)
    genes = list(pairs.keys())
    dimer_flags = []
    for i, g1 in enumerate(genes):
        for g2 in genes[i:]:
            for a_name, a in [('F', pairs[g1]['fwd']), ('R', pairs[g1]['rev'])]:
                for b_name, b in [('F', pairs[g2]['fwd']), ('R', pairs[g2]['rev'])]:
                    if a == b: continue
                    if dimer_3p(a, b, k=8):
                        dimer_flags.append((g1, a_name, g2, b_name))
    if dimer_flags:
        for d in dimer_flags:
            print(f"  ⚠ 3'-dimer: {d[0]}_{d[1]} vs {d[2]}_{d[3]}")
    else:
        print("  ✓ No 8-nt 3'-complementarity among all 12 primers.")

    # Individual primer quality
    print("\n" + "="*90)
    print("Individual primer QC")
    print("="*90)
    for g in genes:
        for d, key in [('F','fwd'), ('R','rev')]:
            s = pairs[g][key]
            gc = gc_content(s)
            md = tm_simple(s)
            nn = tm_nn(s)
            sd = "⚠self-dimer" if self_dimer(s) else "OK"
            print(f"  {g}_{d}: len={len(s)} GC={gc:.1f}% Tm(MD)={md:.1f}°C Tm(NN)={nn:.1f}°C {sd}")

    # Write a compact summary file (autogenerated; see FINAL_primers_report.md for the full write-up)
    with open(P.results('VALIDATION_auto_report.md'), 'w') as f:
        f.write("# Final Cornus DQD/SDH Primer Set (de novo)\n\n")
        f.write("All primers are 20 nt unless otherwise noted. Product sizes are spread\n")
        f.write("120 → 420 bp in ~60 bp steps so the six genes can be resolved together\n")
        f.write("on a 2% agarose gel or by fragment analyzer.\n\n")
        f.write("## Primer pairs\n\n")
        f.write("| Gene | Accession | Direction | Sequence (5'→3') | Region | Position | GC% | Tm(°C) | Amplicon bp |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write("| " + " | ".join([r['Gene'], r['Accession'], r['Dir'],
                                        r["Sequence (5'→3')"], r['Region'],
                                        r['Position'], r['GC%'], r['Tm(°C)'],
                                        r['Amplicon(bp)']]) + " |\n")
        f.write("\n## Specificity\n\n")
        f.write("In-silico exact-match PCR against all six Cornus DQD/SDH transcripts:\n\n")
        if non_specific:
            for src, target, size in non_specific:
                f.write(f"- ⚠ `{src}` pair may also amplify `{target}` ({size} bp)\n")
        else:
            f.write("- ✓ No pair amplifies a non-target Cornus DQD/SDH transcript.\n")
        f.write("\n## Primer-dimer check\n\n")
        if dimer_flags:
            for d in dimer_flags:
                f.write(f"- ⚠ 3'-dimer: `{d[0]}_{d[1]}` vs `{d[2]}_{d[3]}`\n")
        else:
            f.write("- ✓ No 8-nt 3'-complementarity among the 12 primers.\n")
        f.write("\n## Notes on design logic\n\n")
        f.write("- **Both-unique pairs** (genes SDH1.2, 2, 3, 4, 5): both forward and "
                "reverse 20-mers exist only in that transcript and in no other Cornus "
                "DQD/SDH transcript. This gives true gene-specific amplification.\n")
        f.write("- **SDH1.1 special case**: SDH1.1 owns only 38 unique nucleotides "
                "(1567-1604, 19 unique 20-mer starts; see "
                "`data/processed/unique_20mer_spans.tsv`) near the CDS 3' end. Because "
                "that window is too narrow to hold "
                "two 20-mers 120+ bp apart, the SDH1.1 pair uses a *unique forward "
                "anchor* (TGTTCACTGGTGGATTAAGT, inside the 38-nt unique window) plus a "
                "*shared reverse* that also exists in SDH1.2. Specificity still holds "
                "for SDH1.1 because extension cannot be initiated from SDH1.2 without "
                "the unique forward primer, but expect SDH1.1 to be the only clean "
                "amplicon in this pair.\n")
        f.write("\n## Why only ONE primer pair per gene (2, 3, 4, 5)?\n\n")
        f.write("- The paralogs are >90% identical across most of their CDS, so most "
                "20-mers are shared. Only short patches (~20-40 nt each) in the "
                "5'UTR, near variable codons in the CDS, and in the 3'UTR are truly "
                "unique per gene.\n")
        f.write("- Once a unique 20-mer is chosen as the gene-specific anchor, a second "
                "unique 20-mer must lie 120-500 bp away on the same transcript to "
                "complete a pair. That distance filter leaves a small set of options, "
                "usually only a handful per gene, so the pipeline returns one "
                "high-quality pair per gene rather than many alternatives.\n")
        f.write("- Additional pairs *could* be generated per gene by relaxing either "
                "the both-unique rule (allowing one shared partner) or the size "
                "spacing rule; if you want backup pairs, we can rerun with those "
                "relaxed settings.\n")

    print("\nSaved: data/results/VALIDATION_auto_report.md")

if __name__ == '__main__':
    main()
