#!/usr/bin/env python3
"""
Exhaustive (seed-independent) transcriptome scan of the 12 final primers.

Why this exists: BLAST needs an exact 7-mer seed, so a 20-mer primer site with
many scattered mismatches can be missed.  For short primers a brute-force
Hamming scan is both faster and strictly more sensitive, so this script counts
EVERY position in every transcript of the assembly.

Method (numpy, vectorised):
  * concatenate all 29,124 transcripts (48 Mb) into one uint8 buffer,
    non-ACGT characters become a separator that can never match;
  * for each primer, and for its reverse complement, compute the number of
    mismatches in the 20-nt window at every valid start position
    (windows overlapping a transcript boundary are discarded);
  * report the site-count distribution for <=1 ... <=6 mismatches;
  * for every site with <=6 mismatches, recompute the length of the perfectly
    matched run counting back from the primer 3' end (extension-competence);
  * pair forward/reverse sites on the same contig (60-1500 bp) and classify the
    resulting products against the contig annotation.

Outputs (inside primer_blast_check/):
  exhaustive_sites.tsv, exhaustive_products.tsv, TX_EXHAUSTIVE_SCAN.md
"""
import csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P
from collections import defaultdict

import numpy as np

TX_FASTA = P.ext('transcriptome.fasta')
GFF_TX = P.ext('transcriptome_annotation.gff3')
PRIMERS = P.results('FINAL_primers.fasta')
TARGETS = P.raw('Cornus_DQD_SDHs.fasta')

SEP = np.uint8(ord('-'))
MAX_MM = 6
MIN_LEN, MAX_LEN = 60, 1500
FAM = re.compile(r'dehydroquinate|shikimate dehydrogenase|SDH', re.I)

def rc(s):
    c = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return s.translate(c)[::-1]

def load_fasta_bytes(path):
    data = open(path, 'rb').read()
    recs = []
    for chunk in data.split(b'>')[1:]:
        nl = chunk.index(b'\n')
        name = chunk[:nl].split()[0].decode()          # keep header case
        seq = chunk[nl + 1:].replace(b'\n', b'').replace(b' ', b'').upper()
        recs.append((name, seq))
    return recs

def load_primers(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name: out.append((name, ''.join(buf)))
            name = line[1:]; buf = []
        elif line: buf.append(line)
    if name: out.append((name, ''.join(buf)))
    return out

def load_annotation():
    ann = {}
    for gff in (GFF_TX):
        if not os.path.isfile(gff): continue
        for line in open(gff):
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'mRNA': continue
            attrs = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
            prod = attrs.get('product', '')
            if p[0] not in ann or (prod and not ann[p[0]]):
                ann[p[0]] = prod
    return ann

def main():
    recs = load_fasta_bytes(TX_FASTA)
    print(f"transcripts: {len(recs)}")
    L = 20
    lens = np.array([len(s) for _, s in recs], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(lens)])
    total = int(starts[-1])
    pieces = []
    for i, (n_, s) in enumerate(recs):
        arr = np.frombuffer(s, dtype=np.uint8)
        ok = (arr == ord('A')) | (arr == ord('C')) | (arr == ord('G')) | (arr == ord('T'))
        clean = np.where(ok, arr, np.uint8(ord('-'))).astype(np.uint8)
        pieces.append(clean)
        lens[i] = len(clean)
    seq = np.concatenate(pieces)
    del pieces
    n = len(seq)
    print(f"total bases: {n/1e6:.1f} Mb")

    # windows that straddle a transcript boundary are not real sites: mask them
    invalid = np.zeros(n, dtype=bool)
    for b in starts[1:-1]:
        lo = max(0, int(b) - L + 1)
        invalid[lo:int(b)] = True
    npos = n - L + 1
    valid = np.flatnonzero(~invalid[:npos])
    print(f"valid start positions: {len(valid)/1e6:.1f} M")

    contig_starts = starts[:-1]
    def contig_of(pos):
        return int(np.searchsorted(contig_starts, pos, side='right') - 1)

    primers = load_primers(PRIMERS)
    plen = {nm: len(sq) for nm, sq in primers}
    assert all(l == L for l in plen.values()), "script assumes uniform 20-nt primers"

    # target contigs for the six curated genes (megablast-derived, saved by the main check)
    own_contigs = {}
    smap = P.spec('tx_selfmap.tsv')
    if os.path.isfile(smap):
        for row in csv.DictReader(open(smap), delimiter='\t'):
            own_contigs.setdefault(row['curated_transcript'], set()).add(row['contig'])

    ann = load_annotation()

    sites = []
    dist = defaultdict(int)
    for nm, sq in primers:
        gene = nm.rsplit('_', 1)[0]
        side = nm.rsplit('_', 1)[1]
        for orient, q in (('plus', sq), ('minus', rc(sq))):
            qb = np.frombuffer(q.encode(), dtype=np.uint8)
            matches = np.zeros(n - L + 1, dtype=np.uint8)
            for i in range(L):
                matches += (seq[i:i + n - L + 1] == qb[i])
            mm = L - matches
            cand = valid[mm[valid] <= MAX_MM]
            for pos in cand:
                m = int(mm[pos])
                dist[(nm, orient, m)] += 1
                # perfect run from the primer 3' end
                if orient == 'plus':
                    tail = 0
                    for k in range(L - 1, -1, -1):
                        if seq[pos + k] == ord(q[k]): tail += 1
                        else: break
                    strand, s5, s3 = '+', int(pos) + 1, int(pos) + L
                else:
                    # q is already rc(primer): it occupies [pos, pos+L) on the plus strand,
                    # and the primer's 3' end pairs with the LEFT edge (pos)
                    tail = 0
                    for j in range(L):
                        if seq[pos + j] == ord(q[j]): tail += 1
                        else: break
                    strand, s5, s3 = '-', int(pos) + L, int(pos) + 1
                cid = contig_starts[contig_of(int(pos))]
                tx = recs[contig_of(int(pos))][0]
                sites.append({'primer': nm, 'gene': gene, 'side': side, 'strand': strand,
                              'contig': tx, 'start': int(pos) + 1, 'end': int(pos) + L,
                              'mm': m, 'tail3': tail,
                              'own_target': tx in own_contigs.get(gene.replace('CserDQD_SDH', 'CserDQD/SDH'), set()),
                              'family': bool(FAM.search(ann.get(tx, '') or ''))})
            del matches, mm

    with open(P.spec('exhaustive_sites.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['primer','gene','side','strand','contig','start','end',
                                         'mm','tail3','own_target','family'], delimiter='\t')
        w.writeheader()
        for s in sorted(sites, key=lambda x: (x['primer'], x['mm'], x['contig'])):
            w.writerow(s)

    # ---- potential amplicons: F(+) and R(-) on the same contig, 60-1500 bp ----
    byctg = defaultdict(lambda: {'+': [], '-': []})
    for s in sites:
        if s['strand'] == '+': byctg[(s['gene'], s['side'], s['contig'])]['+'].append(s)
        else: byctg[(s['gene'], s['side'], s['contig'])]['-'].append(s)

    prods = []
    for gene in sorted({s['gene'] for s in sites}):
        fs = [s for s in sites if s['gene'] == gene and s['side'] == 'Fwd' and s['strand'] == '+']
        rs = [s for s in sites if s['gene'] == gene and s['side'] == 'Rev' and s['strand'] == '-']
        for a in fs:
            for b in rs:
                if a['contig'] != b['contig']: continue
                size = b['end'] - a['start'] + 1
                if not (MIN_LEN <= size <= MAX_LEN): continue
                mm = a['mm'] + b['mm']
                if mm > 8: continue
                klass = ('TARGET_CONTIG' if (a['own_target'] and b['own_target'])
                         else 'SDH_FAMILY_CONTIG' if (a['family'] or b['family'])
                         else 'UNRELATED_CONTIG')
                prods.append({'pair_gene': gene, 'contig': a['contig'], 'size': size,
                              'mm_fwd': a['mm'], 'mm_rev': b['mm'], 'mm_total': mm,
                              'tail3_fwd': a['tail3'], 'tail3_rev': b['tail3'],
                              'class': klass, 'annotation': (ann.get(a['contig']) or '-')[:100]})
    with open(P.spec('exhaustive_products.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['pair_gene','contig','size','mm_fwd','mm_rev','mm_total',
                                         'tail3_fwd','tail3_rev','class','annotation'], delimiter='\t')
        w.writeheader()
        for p in sorted(prods, key=lambda x: (x['class'], x['pair_gene'], x['size'])):
            w.writerow(p)

    # ---- report ----
    strict = [p for p in prods if p['mm_total'] <= 3 and p['tail3_fwd'] >= 5 and p['tail3_rev'] >= 5]
    ext = [p for p in prods if p['tail3_fwd'] >= 8 and p['tail3_rev'] >= 8 and p['mm_total'] <= 5]
    lines = []
    lines.append("# Exhaustive transcriptome scan (no BLAST seeding)\n")
    lines.append(f"Scanned every one of {n/1e6:.1f} Mb of the assembly at {len(valid)/1e6:.1f} M valid "
                 f"20-nt windows, for all 12 primers on both strands. "
                 f"Reported every site with ≤{MAX_MM} mismatches.\n")
    lines.append("## Sites per primer (mismatch distribution)\n")
    lines.append("| Primer | strand | 0 mm | ≤1 | ≤2 | ≤3 | ≤4 | ≤5 | ≤6 | best 3'-tail at ≤3 mm |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for nm, sq in primers:
        row = []
        for orient in ('plus', 'minus'):
            c = [sum(dist[(nm, orient, m)] for m in range(k + 1)) for k in range(MAX_MM + 1)]
            tails = [s['tail3'] for s in sites if s['primer'] == nm and s['strand'] == orient and s['mm'] <= 3]
            best = max(tails) if tails else 0
            row.append(f"| {nm} | {orient} | " + ' | '.join(str(x) for x in c) + f" | {best} |")
        lines.extend(row)
    lines.append("\n## Candidate products (exhaustive, 60-1500 bp)\n")
    lines.append(f"- all pairs (≤8 total mismatches): **{len(prods)}**")
    lines.append(f"- Primer-BLAST-strict (≤3 mm, perfect 5-nt 3' ends): **{len(strict)}**")
    lines.append(f"- extension-competent relaxed set (≤5 mm, perfect 8-nt 3' ends): **{len(ext)}**\n")
    lines.append("| Pair | Contig | Size | mm F/R | 3'-tail F/R | Class | Annotation |")
    lines.append("|---|---|---|---|---|---|---|")
    for p in sorted(prods, key=lambda x: (x['class'], x['pair_gene'], x['size'])):
        lines.append(f"| {p['pair_gene']} | {p['contig']} | {p['size']} | {p['mm_fwd']}/{p['mm_rev']} | "
                     f"{p['tail3_fwd']}/{p['tail3_rev']} | {p['class']} | {p['annotation'][:60]} |")
    lines.append("\n## Reading\n")
    lines.append("- `mm` counts mismatches over the full 20 nt, so a low-identity partial match is not a "
                 "priming site; the `3'-tail` column is what actually permits extension.\n")
    lines.append("- Every site is reported, including those BLAST could not seed.\n")
    with open(P.spec('TX_EXHAUSTIVE_SCAN.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"\nsites with <=6 mm: {len(sites)}")
    print(f"candidate products (<=8 mm, 60-1500 bp): {len(prods)}")
    print(f"  strict (<=3 mm, 5-nt perfect 3'): {len(strict)}")
    print(f"  extension-competent (<=5 mm, 8-nt perfect 3'): {len(ext)}")
    for p in prods:
        if p not in strict:
            print(f"  loose-only: {p['pair_gene']} {p['contig'].split('_cov')[0]} {p['size']}bp "
                  f"mm {p['mm_fwd']}/{p['mm_rev']} tail {p['tail3_fwd']}/{p['tail3_rev']} {p['class']}")

if __name__ == '__main__':
    main()
