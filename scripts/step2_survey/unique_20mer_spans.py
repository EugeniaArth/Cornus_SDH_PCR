#!/usr/bin/env python3
"""
Where is each transcript actually unique?

For every transcript, count the 20-mers that occur in no other of the six sequences (reverse
complement included, since a shared site can be on either strand), and report the span they
occupy plus the longest contiguous stretch of unique starts. This is what proves the SDH1.1
constraint quoted in reports 02 and 03: its only unique territory is a short window near the
3' end of the CDS, too narrow for two primers 120+ bp apart.

Output: data/processed/unique_20mer_spans.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths as P

K = 20
COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def main():
    seqs = read_fasta(P.require(P.raw('Cornus_DQD_SDHs.fasta'), 'curated cDNA, report 01'))
    out = []
    for gene in sorted(seqs):
        s = seqs[gene]
        other = set()
        for other_gene, o in seqs.items():
            if other_gene == gene:
                continue
            for i in range(len(o) - K + 1):
                w = o[i:i + K]
                other.add(w)
                other.add(rc(w))
        starts = [i for i in range(len(s) - K + 1) if s[i:i + K] not in other]
        if not starts:
            out.append({'gene': gene, 'len': len(s), 'unique_20mers': 0,
                        'first': '', 'last': '', 'span_nt': 0, 'longest_run': 0})
            continue
        run = best = 1
        for a, b in zip(starts, starts[1:]):
            run = run + 1 if b == a + 1 else 1
            best = max(best, run)
        out.append({'gene': gene, 'len': len(s), 'unique_20mers': len(starts),
                    'first': starts[0] + 1, 'last': starts[-1] + K,
                    'span_nt': starts[-1] + K - starts[0], 'longest_run': best})
        print('%-16s %4d unique 20-mers, span %d-%d (%d nt), longest contiguous run %d'
              % (gene, len(starts), starts[0] + 1, starts[-1] + K,
                 starts[-1] + K - starts[0], best))

    cols = ['gene', 'len', 'unique_20mers', 'first', 'last', 'span_nt', 'longest_run']
    with open(P.proc('unique_20mer_spans.tsv'), 'w', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for r in out:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print('\nwrote data/processed/unique_20mer_spans.tsv')


if __name__ == '__main__':
    main()
