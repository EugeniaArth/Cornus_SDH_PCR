#!/usr/bin/env python3
"""
Relaxed permissive window survey — count 20‑mers that cover ALL family members
at various mismatch thresholds.
"""
import os, sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

print("Imported paths.py from:", P.__file__)
print("Project root:", P.ROOT)
print("RAW directory:", P.RAW)
print("Expected CDS:", P.raw('family_cds.fasta'))
print("CDS exists:", os.path.exists(P.raw('family_cds.fasta')))


def load_fasta(path):
    seqs = OrderedDict()
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper())
        if name: seqs[name] = ''.join(buf)
    return seqs

def hamming(a, b):
    return sum(1 for x,y in zip(a,b) if x!=y)

def survey_family(fam, K=20):
    path = P.proc(f'{fam}_codon_aln.fasta')
    if not os.path.exists(path):
        return None
    align = load_fasta(path)
    labels = list(align.keys())
    if len(labels) < 2:
        return None
    # Use longest CDS as reference
    ref_label = max(labels, key=lambda l: len(align[l]))
    ref = align[ref_label]
    n = len(ref)
    n_mem = len(labels)
    others = [(l, align[l]) for l in labels if l != ref_label]

    results = {}
    for max_mm in range(6):
        count = 0
        for start in range(n - K + 1):
            ref_win = ref[start:start+K]
            if '-' in ref_win or 'N' in ref_win:
                continue
            covered = True
            for l, o in others:
                if start + K > len(o):
                    covered = False
                    break
                o_win = o[start:start+K]
                if '-' in o_win or 'N' in o_win:
                    covered = False
                    break
                if hamming(ref_win, o_win) > max_mm:
                    covered = False
                    break
            if covered:
                count += 1
        results[max_mm] = count
    return results, ref_label, n, n_mem

def main():
    families = ['ACT','BTUB','EF1a','GAPDH','H3','eIF5A','polyUb']
    print(f"{'Family':8s} {'Members':6s} {'RefLen':6s} {'mm0':5s} {'mm1':5s} {'mm2':5s} {'mm3':5s} {'mm4':5s} {'mm5':5s}")
    print("-"*60)
    for fam in families:
        r = survey_family(fam)
        if r is None:
            print(f"{fam:8s}  — codon alignment not found")
            continue
        stats, ref, n_cols, n_mem = r
        row = f"{fam:8s} {n_mem:3d}     {n_cols:4d}  "
        for mm in range(6):
            row += f"{stats.get(mm,0):4d} "
        print(row)

if __name__ == '__main__':
    main()
