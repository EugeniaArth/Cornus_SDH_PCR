#!/usr/bin/env python3
"""
Verification of the recommended backup CserDQD/SDH3 pairs.

Four configurations are compared on identical terms:
  * in-service  : the six pairs of FINAL_primers.fasta (baseline)
  * #1          : top-ranked backup
  * #2          : backup with zero foreign sites even at <=8 mm / >=3-nt 3' run
  * #3          : highest-Tm Tier A backup (the option used when the mitigation is
                  a hotter annealing temperature)

Stage 1 - six curated transcripts, all 36 forward x reverse combinations, five rules.
Stage 2 - whole transcriptome (29,124 contigs, 48.1 Mb), seed-free scan of every
          20-nt window on both strands, sites recorded up to 6 mismatches per primer,
          then product pairing (60-1000 bp) under the same rules used in
          TX_EXHAUSTIVE_SCAN.md so the numbers are directly comparable.

Products are classified as
  designed   : both primers perfect (0 mm) on the contig that carries their own gene
  on-family  : at least one mismatch, contig maps to one of the six curated DQD/SDH
               transcripts (i.e. a paralog cross-product)
  off-family : at least one mismatch, any other contig in the assembly

Full product lists go to backup_tx_products.tsv / backup_curated_products.tsv;
the markdown holds the counts and the worst examples only.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'step3_design'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from design_final_pairs import load_fasta, rc
from design_sdh3_backup import scan, read_final, L, GENE, ALT_LO, ALT_HI, TXOME, tail_plus, tail_minus

SEQS = load_fasta(P.raw('Cornus_DQD_SDHs.fasta'))
ARR = {g: np.frombuffer(s.encode(), dtype=np.uint8) for g, s in SEQS.items()}
GENES = sorted(ARR)
SHORT = {g: g.split('/')[-1] for g in GENES}
NAMES = {g: "CserDQD_SDH" + g.split('/')[-1].replace('SDH', '') for g in GENES}
DESIGNED_SIZE = {'CserDQD/SDH1.1': 120, 'CserDQD/SDH1.2': 180, 'CserDQD/SDH2': 240,
                 'CserDQD/SDH3': 300, 'CserDQD/SDH4': 360, 'CserDQD/SDH5': 420}
RULES = [("designed        (both primers 0 mm)", 0, 0),
         ("strict          (total ≤3 mm, 3′ runs ≥5 nt)", 3, 5),
         ("competent       (total ≤6 mm, 3′ runs ≥5 nt)", 6, 5),
         ("long anchor     (total ≤8 mm, 3′ runs ≥8 nt)", 8, 8),
         ("any             (total ≤8 mm, no 3′ requirement)", 8, 0)]
CUR_RULES = [("strict          (total ≤3 mm, 3′ runs ≥5 nt)", 3, 5),
             ("competent       (total ≤6 mm, 3′ runs ≥5 nt)", 6, 5),
             ("long anchor     (total ≤8 mm, 3′ runs ≥8 nt)", 8, 8),
             ("permissive      (each ≤6 mm, 3′ runs ≥5 nt)", 12, 5),
             ("very permissive (each ≤8 mm, 3′ runs ≥3 nt)", 16, 3)]

fam = list(csv.DictReader(open(P.spec('tx_selfmap.tsv')), delimiter='\t'))
FAMILY = {r['contig'] for r in fam}
GENE_OF = defaultdict(set)
for r in fam:
    GENE_OF[r['contig']].add(r['curated_transcript'].split('/')[-1])
SHORT2SIZE = {v: DESIGNED_SIZE[k] for k, v in SHORT.items()}


def load_txome():
    recs, name, buf = [], None, []
    for line in open(TXOME):
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if name: recs.append((name, ''.join(buf)))
            name = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if name: recs.append((name, ''.join(buf)))
    seqall = np.frombuffer(''.join(s for _, s in recs).encode(), dtype=np.uint8)
    starts = np.zeros(len(recs), dtype=np.int64)
    acc = 0
    for k, (_, s) in enumerate(recs):
        starts[k] = acc; acc += len(s)
    n = len(seqall)
    invalid = np.zeros(n, dtype=bool)
    for b in starts[1:]:
        invalid[max(0, int(b) - L + 1):int(b)] = True
    valid = np.flatnonzero(~invalid[:n - L + 1])
    return recs, seqall, starts, valid, n

CACHE = {}
def tx_sites(p, seqall, starts, valid, n, recs, maxmm=6):
    if (p, maxmm) in CACHE: return CACHE[(p, maxmm)]
    qp = np.frombuffer(p.encode(), np.uint8)
    qr = np.frombuffer(rc(p).encode(), np.uint8)
    hits = []
    for orient, qb in (('+', qp), ('-', qr)):
        m = np.zeros(n - L + 1, dtype=np.int16)
        for i in range(L):
            m += (seqall[i:i + n - L + 1] == qb[i])
        mm = L - m
        for pos in valid[mm[valid] <= maxmm]:
            pos = int(pos)
            tail = (tail_plus(seqall, pos, qp) if orient == '+' else tail_minus(seqall, pos, qr))
            ci = int(np.searchsorted(starts, pos, side='right') - 1)
            hits.append({'tx': recs[ci][0], 'strand': orient, 'start': pos + 1,
                         'end': pos + L, 'mm': int(mm[pos]), 'tail': tail})
    CACHE[(p, maxmm)] = hits
    return hits

def all_pairings(pairset, sites, max_total, min_run):
    """Every forward x reverse combination present in the tube (36 for a 6-plex):
    F of pair x with R of pair y, on every transcript, total mismatches <= max_total,
    both perfect 3' runs >= min_run, product ALT_LO-ALT_HI bp."""
    genes = list(pairset)
    prep = {}
    for g in genes:
        fwd, rev = defaultdict(list), defaultdict(list)
        for s in sites[pairset[g][0]]:
            if s['strand'] == '+' and s['tail'] >= min_run: fwd[s['tx']].append(s)
        for s in sites[pairset[g][1]]:
            if s['strand'] == '-' and s['tail'] >= min_run: rev[s['tx']].append(s)
        prep[g] = (fwd, rev)
    out = []
    for x in genes:
        for y in genes:
            bf, br = prep[x][0], prep[y][1]
            for t in set(bf) & set(br):
                for a in bf[t]:
                    for b in br[t]:
                        tot = a['mm'] + b['mm']
                        sz = b['end'] - a['start'] + 1
                        if tot > max_total or not (ALT_LO <= sz <= ALT_HI): continue
                        out.append({'pair': SHORT[x], 'rpair': SHORT[y],
                                    'target': SHORT[y] if x == y else SHORT[x] + '×' + SHORT[y],
                                    'contig': t, 'size': sz, 'f_mm': a['mm'], 'r_mm': b['mm'],
                                    'total': tot, 'f_tail': a['tail'], 'r_tail': b['tail'],
                                    'cross': x != y})
    return out

def classify(pr):
    if pr['f_mm'] == 0 and pr['r_mm'] == 0:
        return 'designed' if not pr['cross'] else 'cross-exact'
    return 'on-family' if pr['contig'] in FAMILY else 'off-family'

def main():
    cands = list(csv.DictReader(open(P.spec('SDH3_BACKUP_PAIRS.tsv')), delimiter='\t'))
    A = [c for c in cands if c['tier'] == 'A']
    picks = {'in-service (baseline)': None,
             '#1 top-ranked backup': A[0]}
    zero = [c for c in A if int(c['very_permissive_foreign']) == 0]
    if zero: picks['#2 zero foreign sites even at ≤8 mm/3-nt run'] = zero[0]
    picks['#3 highest-Tm backup (hot annealing)'] = sorted(
        A, key=lambda c: -(float(c['fwd_tm']) + float(c['rev_tm'])))[0]

    finals = read_final(P.results('FINAL_primers.fasta'))
    in_service = {g: (finals[NAMES[g] + '_Fwd'], finals[NAMES[g] + '_Rev']) for g in GENES}
    print("loading transcriptome ...")
    recs, seqall, starts, valid, n = load_txome()
    print(f"  {len(recs)} contigs, {n/1e6:.1f} Mb, {len(FAMILY)} DQD/SDH-family contigs")

    md = ["# Backup CserDQD/SDH3 pair — verification against the in-service set\n",
          "The same exhaustive machinery is run on four configurations: the six pairs in "
          "service and three SDH3 replacements. All 36 forward×reverse combinations of each "
          "6-plex are examined, on the six curated transcripts and on all 48.1 Mb of the "
          "assembly (seed-free, every 20-nt window, both strands).\n",
          "Note: only **five** of the six genes have their own contig in the clustered "
          "assembly (SDH1.1 is collapsed into NODE_10980 with SDH1.2), so the maximum number "
          "of detectable designed amplicons transcriptome-wide is 5, not 6.\n"]
    all_prod_rows, cur_rows = [], []

    for label, c in picks.items():
        six = dict(in_service)
        if c is not None:
            six[GENE] = (c['fwd_seq'], c['rev_seq'])
        md.append(f"\n## {label}\n")
        if c is not None:
            md.append(f"- SDH3 primers under test: F `{c['fwd_seq']}` ({c['fwd_pos']}, "
                      f"{c['fwd_tm']} °C) / R `{c['rev_seq']}` ({c['rev_pos']}, {c['rev_tm']} °C), "
                      f"{c['product']} bp")
        else:
            md.append("- SDH3 primers in service: F `%s` (%s) / R `%s` (%s), 300 bp"
                      % (in_service[GENE][0], '738-757', in_service[GENE][1], '1018-1037'))

        # ---------------- stage 1: curated ----------------
        sites = {}
        for g in GENES:
            for p in six[g]:
                if p not in sites: sites[p] = scan(p, ARR)
        md.append("\n### Stage 1 — six curated transcripts, all 36 forward×reverse combinations\n")
        md.append("| rule | designed amplicons | other products | foreign-paralog / cross-pair detail |")
        md.append("|---|---|---|---|")
        for rname, maxtot, minrun in CUR_RULES:
            pr = all_pairings(six, sites, maxtot, minrun)
            des, off = set(), set()
            for q in pr:
                cur_rows.append({'config': label, 'rule': rname.split('(')[0].strip(),
                                 'fwd_of_pair': q['pair'], 'rev_of_pair': q['rpair'],
                                 'transcript': q['contig'].split('/')[-1],
                                 'cross': 'yes' if q['cross'] else 'no', 'size': q['size'],
                                 'f_mm': q['f_mm'], 'r_mm': q['r_mm'],
                                 'f_tail': q['f_tail'], 'r_tail': q['r_tail']})
                tag = (f"F:{q['pair']}+R:{q['rpair']} on {q['contig'].split('/')[-1]} "
                       f"{q['size']} bp (F {q['f_mm']} mm/3′{q['f_tail']}, "
                       f"R {q['r_mm']} mm/3′{q['r_tail']})" + (' *[cross]*' if q['cross'] else ''))
                own = (not q['cross'] and q['f_mm'] == 0 and q['r_mm'] == 0
                       and q['size'] == SHORT2SIZE.get(q['contig'].split('/')[-1]))
                (des if own else off).add(tag)
            md.append(f"| {rname} | {len(des)} | {len(off)} | "
                      f"{'; '.join(sorted(off)[:4]) if off else '—'}"
                      + (f" (+{len(off)-4} more)" if len(off) > 4 else "") + " |")
            print(f"  {label} | curated {rname.split('(')[0].strip()}: designed={len(des)} "
                  f"other={len(off)}")

        # ---------------- stage 2: transcriptome ----------------
        primers = sorted({p for g in GENES for p in six[g]})
        tsites = {p: tx_sites(p, seqall, starts, valid, n, recs, 6) for p in primers}
        md.append("\n### Stage 2 — whole transcriptome, every window, same rules as TX_EXHAUSTIVE_SCAN.md\n")
        md.append("| primer | sites ≤2 mm | ≤4 mm | ≤6 mm | longest 3′ run at ≤4 mm |")
        md.append("|---|---|---|---|---|")
        for p in primers:
            s = tsites[p]
            tag = ' ← **backup**' if p in six[GENE] else ''
            t4 = [x['tail'] for x in s if x['mm'] <= 4]
            md.append(f"| `{p}`{tag} | {sum(1 for x in s if x['mm'] <= 2)} | "
                      f"{sum(1 for x in s if x['mm'] <= 4)} | {len(s)} | {max(t4) if t4 else 0} |")
        md.append("\n| rule | designed | cross-pair exact | on-family (paralog) | off-family (unrelated contigs) | worst example |")
        md.append("|---|---|---|---|---|---|")
        for rname, maxtot, minrun in RULES:
            pr = all_pairings(six, tsites, maxtot, minrun)
            d = [q for q in pr if classify(q) == 'designed']
            xe = [q for q in pr if classify(q) == 'cross-exact']
            on = [q for q in pr if classify(q) == 'on-family']
            off = [q for q in pr if classify(q) == 'off-family']
            for q in pr:
                all_prod_rows.append({'config': label, 'rule': rname.split('(')[0].strip(),
                                      'fwd_of_pair': q['pair'], 'rev_of_pair': q['rpair'],
                                      'contig': q['contig'], 'size': q['size'],
                                      'f_mm': q['f_mm'], 'r_mm': q['r_mm'], 'f_tail': q['f_tail'],
                                      'r_tail': q['r_tail'], 'class': classify(q),
                                      'family_gene': '/'.join(sorted(GENE_OF.get(q['contig'], ['?'])))})
            worst_ex = '—'
            if on or off:
                w = sorted(on + off, key=lambda z: (z['total'], -min(z['f_tail'], z['r_tail'])))[0]
                worst_ex = (f"F:{w['pair']}+R:{w['rpair']}→{w['contig'][:20]}… {w['size']} bp, "
                            f"{w['total']} mm, 3′ {w['f_tail']}/{w['r_tail']}")
            uniq = lambda lst: len(set((q['pair'], q['rpair'], q['contig'], q['size']) for q in lst))
            md.append(f"| {rname} | {uniq(d)} | {uniq(xe)} | {uniq(on)} | {uniq(off)} | {worst_ex} |")
            print(f"  {label} | tx {rname.split('(')[0].strip()}: designed={uniq(d)} "
                  f"cross-exact={uniq(xe)} on-family={uniq(on)} off-family={uniq(off)}")
        md.append("")

    with open(P.spec('backup_tx_products.tsv'), 'w', newline='') as fh:
        cols = ['config', 'rule', 'fwd_of_pair', 'rev_of_pair', 'contig', 'family_gene',
                'size', 'f_mm', 'r_mm', 'f_tail', 'r_tail', 'class']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        seen = set()
        for r in sorted(all_prod_rows, key=lambda x: (x['config'], x['rule'], x['class'], x['fwd_of_pair'], x['rev_of_pair'])):
            k = tuple(r.values())
            if k in seen: continue
            seen.add(k); w.writerow(r)
    with open(P.spec('backup_curated_products.tsv'), 'w', newline='') as fh:
        cols = ['config', 'rule', 'fwd_of_pair', 'rev_of_pair', 'transcript', 'cross',
                'size', 'f_mm', 'r_mm', 'f_tail', 'r_tail']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        seen = set()
        for r in sorted(cur_rows, key=lambda x: (x['config'], x['rule'], x['fwd_of_pair'], x['rev_of_pair'],
                                             x['size'], tuple(str(v) for v in x.values()))):  # full tiebreak = reproducible bytes
            k = tuple(r.values())
            if k in seen: continue
            seen.add(k); w.writerow(r)
    md.append("\nFull product tables: `backup_curated_products.tsv` (six transcripts) and "
              "`backup_tx_products.tsv` (whole transcriptome).")
    open(P.spec('SDH3_BACKUP_VERIFY.md'), 'w').write("\n".join(md) + "\n")
    print("wrote SDH3_BACKUP_VERIFY.md, backup_curated_products.tsv, backup_tx_products.tsv")

if __name__ == '__main__':
    main()
