#!/usr/bin/env python3
"""
Sensitivity / stress test of the transcriptome-wide specificity result.

The main check (run_transcriptome_check.py) used Primer-BLAST-default-like
cut-offs.  Here the cut-offs are deliberately RELAXED to see whether any
off-target product appears once the stringency is lowered:

  variant A: <=4 total mismatches, product 60-1000 bp
  variant B: <=6 total mismatches, product 60-1000 bp, 3' perfect match only 4 nt

Also reports, for every transcriptome contig annotated as a DQD/SDH-family
protein, how many primer hits it carries - i.e. whether other family members
exist in the assembly at all and whether primers come close to them.
"""
import csv, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

DB = P.blastdb('cornus_tx')
PRIMERS = P.results('FINAL_primers.fasta')
GFF_TX = P.ext('transcriptome_annotation.gff3')

def run(cmd, inp=None):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True)
    if r.returncode: sys.exit(r.stderr)
    return r.stdout

def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if name: recs.append((name, ''.join(buf)))
            name = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if name: recs.append((name, ''.join(buf)))
    return recs

primers = read_fasta(PRIMERS)
plen = {n: len(s) for n, s in primers}
cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
raw = run(['blastn', '-task', 'blastn-short', '-word_size', '7', '-evalue', '1000',
           '-perc_identity', '60', '-max_target_seqs', '5000', '-db', DB,
           '-outfmt', f'6 {cols}'],
          '\n'.join(f'>{n}\n{s}' for n, s in primers))

aln = []
for line in raw.strip().splitlines():
    if not line: continue
    p = line.split('\t')
    rec = {'primer': p[0], 'tx': p[1], 'alen': int(p[3]), 'mm': int(p[4]),
           'qstart': int(p[6]), 'qend': int(p[7]), 'sstart': int(p[8]),
           'send': int(p[9]), 'qseq': p[12], 'sseq': p[13]}
    aln.append(rec)
print(f"raw hits at perc_identity>=60: {len(aln)}")

def is_full(r, name): return (r['qstart'] == 1 and r['qend'] == plen[name]
                              and r['alen'] == plen[name])
def tail_perfect(r, need):
    g = 0
    for qb, sb in reversed(list(zip(r['qseq'], r['sseq']))):
        if qb.upper() == sb.upper(): g += 1
        else: break
    return g >= need

pairs = {}
for n, s in primers:
    pairs.setdefault(n.rsplit('_', 1)[0], {})[n.rsplit('_', 1)[1]] = n

def products(max_mm, lo, hi, need3):
    out = []
    for gene, d in sorted(pairs.items()):
        fh = [r for r in aln if r['primer'] == d['Fwd'] and is_full(r, d['Fwd'])
              and r['sstart'] < r['send'] and r['mm'] <= max_mm and tail_perfect(r, need3)]
        rh = [r for r in aln if r['primer'] == d['Rev'] and is_full(r, d['Rev'])
              and r['sstart'] > r['send'] and r['mm'] <= max_mm and tail_perfect(r, need3)]
        for a in fh:
            for b in rh:
                if a['tx'] != b['tx']: continue
                size = b['sstart'] - a['sstart'] + 1
                if not (lo <= size <= hi): continue
                if a['mm'] + b['mm'] > max_mm: continue
                out.append((gene, a['tx'], size, a['mm'], b['mm']))
    return out

# annotation lookup
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

FAM = re.compile(r'dehydroquinate|shikimate|dehydrogenase, chloroplastic|'
                 r'SDH', re.I)

lines = []
for tag, (mm, lo, hi, n3) in {
        'A (<=4 mm, 60-1000 bp, 3-prime 5 nt perfect)': (4, 60, 1000, 5),
        'B (<=6 mm, 60-1000 bp, 3-prime 4 nt perfect)': (6, 60, 1000, 4),
        'C (<=8 mm, 60-1500 bp, 3-prime 3 nt perfect)': (8, 60, 1500, 3),
        }.items():
    ps = products(mm, lo, hi, n3)
    fam = [p for p in ps if FAM.search(ann.get(p[1], ''))]
    lines.append(f"## Variant {tag}\n")
    lines.append(f"total candidate products: {len(ps)}  "
                 f"(on DQD/SDH-annotated contigs: {len(fam)}; on other contigs: {len(ps)-len(fam)})\n")
    for g, tx, size, m1, m2 in ps:
        cls = 'SDH-family' if FAM.search(ann.get(tx, '')) else 'OTHER'
        lines.append(f"  - {g} | {tx} | {size} bp | mm {m1}+{m2} | {cls} | "
                     f"{(ann.get(tx) or '-')[:80]}\n")
    lines.append("")

# family content of the assembly
fam_contigs = {k: v for k, v in ann.items() if FAM.search(v or '')}
hit_counts = {}
for r in aln:
    if r['tx'] in fam_contigs and is_full(r, r['primer']) and r['mm'] <= 4:
        hit_counts.setdefault(r['tx'], 0)
        hit_counts[r['tx']] += 1
lines.append("## DQD/SDH-family-annotated contigs in the assembly\n")
lines.append(f"contigs with a family annotation: {len(fam_contigs)}\n")
lines.append("contigs carrying a near-full-length (20/20 aligned, <=4 mm) primer site: "
             f"{len(hit_counts)}\n")
for tx, n in sorted(hit_counts.items(), key=lambda x: -x[1]):
    lines.append(f"  - {tx}: {n} sites | {(fam_contigs.get(tx) or '')[:70]}\n")

with open(P.spec('TX_SENSITIVITY.md'), 'w') as f:
    f.write("# Sensitivity test of the transcriptome-wide specificity check\n\n")
    f.write("Thresholds deliberately relaxed vs the main check to see if off-target products "
            "appear below the cut-offs. blastn-short, word_size 7, perc_identity 60, "
            "full-length (20/20) primer alignment required.\n\n")
    f.write('\n'.join(lines))
print('\n'.join(lines[:6]))
print(f"\nWrote TX_SENSITIVITY.md ({len(products(4,60,1000,5))} / "
      f"{len(products(6,60,1000,4))} / {len(products(8,60,1500,3))} products in variants A/B/C)")
