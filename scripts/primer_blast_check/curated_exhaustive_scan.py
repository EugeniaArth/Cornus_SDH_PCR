#!/usr/bin/env python3
"""
Exhaustive (seed-independent) mismatch scan of the 12 final primers against the
SIX CURATED DQD/SDH transcripts (DQD_SDH_possible_nt.fasta = PZ518666-71),
i.e. the templates that will actually be used in the wet lab.

This complements the transcriptome scan: CserDQD/SDH1.1 has no separate contig
in the clustered assembly (it collapses into NODE_10980 with SDH1.2), so its
specificity can only be judged on the curated sequences.

For every primer x transcript combination the script reports the best 20-nt
window (full-length, no gaps) with its mismatch count, the length of the
perfectly matched run counting back from the primer 3' end, and the mismatch
positions.  Then it builds candidate amplicons under three rule sets.

Outputs: curated_exhaustive_sites.tsv, curated_exhaustive_products.tsv,
         TX_CURATED_EXHAUSTIVE.md
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
TX = os.path.join(PARENT, 'DQD_SDH_possible_nt.fasta')
PR = os.path.join(PARENT, 'FINAL_primers.fasta')
L = 20

def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if name: recs.append((name, ''.join(buf).upper()))
            name = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if name: recs.append((name, ''.join(buf).upper()))
    return recs

primers = read_fasta(PR)
txs = read_fasta(TX)
genes = [n for n, _ in txs]

# ---- exhaustive window scan ------------------------------------------------
MAXR = 6
sites = []
for pname, pseq in primers:
    gene = pname.rsplit('_', 1)[0]
    side = pname.rsplit('_', 1)[1]
    for tname, tseq in txs:
        arr = np.frombuffer(tseq.encode(), dtype=np.uint8)
        n = len(arr)
        if n < L: continue
        for orient in ('+', '-'):
            q = pseq if orient == '+' else rc(pseq)
            qb = np.frombuffer(q.encode(), dtype=np.uint8)
            match = np.zeros(n - L + 1, dtype=np.int16)
            for i in range(L):
                match += (arr[i:i + n - L + 1] == qb[i])
            mm = L - match
            best = int(mm.argmin())      # window with the FEWEST mismatches
            cand = set(int(p) for p in np.flatnonzero(mm <= MAXR))
            cand.add(int(best))          # always report the best window, even if mm > MAXR
            for pos in sorted(cand):
                tail = 0
                if orient == '+':
                    for k in range(L - 1, -1, -1):
                        if arr[pos + k] == qb[k]: tail += 1
                        else: break
                    pos_mm = [i + 1 for i in range(L) if arr[pos + i] != qb[i]]
                else:
                    for j in range(L):
                        if arr[pos + j] == qb[j]: tail += 1
                        else: break
                    pos_mm = [L - i for i in range(L) if arr[pos + i] != qb[i]]
                sites.append({'primer': pname, 'gene': gene, 'side': side, 'strand': orient,
                              'transcript': tname, 'start': int(pos) + 1, 'end': int(pos) + L,
                              'mm': int(mm[pos]), 'tail3': tail,
                              'is_best': int(pos == best),
                              'mm_positions': ','.join(map(str, pos_mm)) or '-'})

with open(os.path.join(HERE, 'curated_exhaustive_sites.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(sites[0].keys()), delimiter='\t')
    w.writeheader()
    for s in sorted(sites, key=lambda x: (x['primer'], x['transcript'], x['mm'], x['start'])):
        w.writerow(s)

# ---- candidate amplicons ---------------------------------------------------
def build(max_mm, need_tail, lo, hi):
    out = []
    pairs = {}
    for pname, pseq in primers:
        g = pname.rsplit('_', 1)[0]
        pairs.setdefault(g, {})[pname.rsplit('_', 1)[1]] = pname
    for g, d in sorted(pairs.items()):
        for tname in genes:
            fs = [s for s in sites if s['primer'] == d['Fwd'] and s['transcript'] == tname
                  and s['strand'] == '+' and s['mm'] <= max_mm and s['tail3'] >= need_tail]
            rs = [s for s in sites if s['primer'] == d['Rev'] and s['transcript'] == tname
                  and s['strand'] == '-' and s['mm'] <= max_mm and s['tail3'] >= need_tail]
            for a in fs:
                for b in rs:
                    size = b['end'] - a['start'] + 1
                    if not (lo <= size <= hi): continue
                    if a['mm'] + b['mm'] > max_mm: continue
                    own = (tname == g.replace('CserDQD_SDH', 'CserDQD/SDH'))
                    out.append({'pair_gene': g, 'transcript': tname, 'size': size,
                                'mm_fwd': a['mm'], 'mm_rev': b['mm'],
                                'tail_fwd': a['tail3'], 'tail_rev': b['tail3'],
                                'class': 'TARGET' if own else 'CROSS-PRODUCT'})
    return out

rules = {'strict  (<=3 mm, perfect 5-nt 3\' end, 80-600 bp)': build(3, 5, 80, 600),
         'medium  (<=4 mm, perfect 8-nt 3\' end, 60-1000 bp)': build(4, 8, 60, 1000),
         'loose   (<=6 mm, perfect 8-nt 3\' end, 60-1000 bp)': build(6, 8, 60, 1000)}
with open(os.path.join(HERE, 'curated_exhaustive_products.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['rule_set', 'pair_gene', 'transcript', 'size', 'mm_fwd', 'mm_rev',
                'tail_fwd', 'tail_rev', 'class'])
    for rule, ps in rules.items():
        for p in ps:
            w.writerow([rule] + [p[k] for k in ('pair_gene','transcript','size','mm_fwd',
                                                'mm_rev','tail_fwd','tail_rev','class')])

# ---- report ---------------------------------------------------------------
lines = ["# Exhaustive mismatch scan on the six curated DQD/SDH transcripts\n",
         f"All {len(primers)} primers x {len(txs)} transcripts x 2 strands; every 20-nt window was "
         f"compared (no BLAST seeding, no gaps). Every site with ≤{MAXR} mismatches is in "
         "`curated_exhaustive_sites.tsv`; the table below lists only the best window "
         "(fewest mismatches) per primer x transcript.\n",
         "## Best full-length window per primer x transcript\n",
         "| Primer | Transcript | mm (of 20) | 3'-perfect run | mismatch positions (primer 5'->3') |",
         "|---|---|---|---|---|"]
for s in sorted(sites, key=lambda x: (x['primer'], x['transcript'])):
    if not s['is_best']: continue
    star = ' **(designed site)**' if s['mm'] == 0 else ''
    lines.append(f"| {s['primer']} | {s['transcript']} | {s['mm']} | {s['tail3']} | "
                 f"{s['mm_positions']}{star} |")
lines.append("")
lines.append("## Candidate amplicons under three rule sets\n")
for rule, ps in rules.items():
    cross = [p for p in ps if p['class'] == 'CROSS-PRODUCT']
    lines.append(f"### {rule}")
    lines.append(f"- products: {len(ps)}  |  cross-products (pair on a non-target transcript): {len(cross)}")
    for p in ps:
        lines.append(f"  - {p['pair_gene']} → **{p['transcript']}** {p['size']} bp "
                     f"(mm {p['mm_fwd']}/{p['mm_rev']}, 3'-run {p['tail_fwd']}/{p['tail_rev']}) "
                     f"[{p['class']}]")
    lines.append("")
with open(os.path.join(HERE, 'TX_CURATED_EXHAUSTIVE.md'), 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"scanned {len(primers)} primers x {len(txs)} transcripts, best-window table written")
for rule, ps in rules.items():
    cross = [p for p in ps if p['class'] == 'CROSS-PRODUCT']
    print(f"  {rule}: products={len(ps)} cross={len(cross)}")
    for p in cross:
        print(f"      CROSS: {p['pair_gene']} -> {p['transcript']} {p['size']} bp "
              f"mm {p['mm_fwd']}/{p['mm_rev']} tail {p['tail_fwd']}/{p['tail_rev']}")
