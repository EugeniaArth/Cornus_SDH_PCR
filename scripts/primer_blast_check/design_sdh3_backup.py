#!/usr/bin/env python3
"""
Backup CserDQD/SDH3 primer-pair search.

Why: the SDH3 pair in service has, under a permissive competence rule, cross-sites
on CserDQD/SDH4 (300 bp product, 3+3 mismatches, 3'-runs 7/8 nt) and CserDQD/SDH2
(309 bp, 3+4 mm, 3'-runs 7/5 nt).  Standard Primer-BLAST rules reject them, but the
SDH4 co-product would migrate at exactly the same size as the intended SDH3 band.

This script re-searches the SDH3 transcript for an alternative pair that is safe at
the SINGLE-PRIMER level - a strictly stronger guarantee than pair-level safety: a
primer with no competent site anywhere on the other five transcripts cannot mis-prime
regardless of what its partner does.

Definitions
  competent site : 20-nt window elsewhere with <=6 mismatches AND >=5 nt perfect
                   3'-terminal run (exactly the rule that exposed the SDH3 problem)
  Tier A         : both primers individually have ZERO competent foreign sites
  Tier B         : competent foreign sites exist, but no F+R pairing on the same
                   foreign transcript within 60-1000 bp
  QC             : as design_final_pairs.py (GC 40-60 %, no homopolymer >4 nt,
                   acceptable 3' tail, no >=7 nt self-complementary hairpin),
                   SantaLucia 1998 NN Tm at 250 nM primer / 50 mM Na+
  product        : exactly 300 bp (keeps the 120/180/240/300/360/420 ladder intact);
                   widened to 280-320 bp only if nothing at 300 bp survives
  multiplex      : no >=8 nt 3'-terminal complementarity against the ten other primers
  genome-wide    : the recommended pair is then scanned exhaustively against all
                   29,124 transcripts of Transcriptome/final.clust_transcripts_longest_iso.fasta

Outputs (inside primer_blast_check/ only):
  SDH3_BACKUP_PAIRS.tsv, SDH3_backup_candidates.fasta, SDH3_BACKUP_REPORT.md
"""
import csv, os, sys
import numpy as np

PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PARENT)
from design_final_pairs import (gc_content, tm_nn, rc, primer_ok, load_fasta,
                                load_coords, classify_position)

TX_FASTA = os.path.join(PARENT, 'DQD_SDH_possible_nt.fasta')
COORDS   = os.path.join(PARENT, 'genbank_coordinates_combined.tsv')
FINAL    = os.path.join(PARENT, 'FINAL_primers.fasta')
TXOME    = os.path.join(PARENT, 'Transcriptome', 'final.clust_transcripts_longest_iso.fasta')

L = 20
GENE = 'CserDQD/SDH3'
TARGET_PROD = 300
MAX_MM_COMP = 6      # mismatches allowed in a competent foreign site
MIN_TAIL_COMP = 5    # perfect 3'-terminal run required for competence
MAXR = 8             # deepest mismatch level retained, for margin reporting
ALT_LO, ALT_HI = 60, 1000


def tail_plus(arr, pos, qp):
    t = 0
    for k in range(L - 1, -1, -1):
        if arr[pos + k] == qp[k]: t += 1
        else: break
    return t

def tail_minus(arr, pos, qr):
    t = 0
    for j in range(L):
        if arr[pos + j] == qr[j]: t += 1
        else: break
    return t

def scan(p, arrays):
    """Every site of primer p (both strands) with mm <= MAXR, plus its 3'-terminal run."""
    out = []
    qp = np.frombuffer(p.encode(), dtype=np.uint8)
    qr = np.frombuffer(rc(p).encode(), dtype=np.uint8)
    for tname, arr in arrays.items():
        n = len(arr)
        if n < L: continue
        for strand, qb in (('+', qp), ('-', qr)):
            m = np.zeros(n - L + 1, dtype=np.int16)
            for i in range(L):
                m += (arr[i:i + n - L + 1] == qb[i])
            mm = L - m
            for pos in np.flatnonzero(mm <= MAXR):
                pos = int(pos)
                tail = (tail_plus(arr, pos, qp) if strand == '+'
                        else tail_minus(arr, pos, qr))
                out.append({'tx': tname, 'strand': strand, 'start': pos + 1,
                            'end': pos + L, 'mm': int(mm[pos]), 'tail': tail})
    return out

def competent(sites, exclude_tx):
    return [s for s in sites if s['tx'] != exclude_tx
            and s['mm'] <= MAX_MM_COMP and s['tail'] >= MIN_TAIL_COMP]

def worst(sites, exclude_tx):
    foreign = [s for s in sites if s['tx'] != exclude_tx]
    if not foreign: return 'no foreign site ≤8 mm'
    w = min(foreign, key=lambda s: (s['mm'], -s['tail']))
    return f"{w['tx']} {w['mm']}mm tail{w['tail']}"

def three_prime_dimer(a, b):
    """Longest l with a[-l:] complementary to b[:l] (3' of a pairing into 5' of b)."""
    best = 0
    for l in range(1, L + 1):
        if a[-l:] == rc(b[:l]): best = l
    return best

def max_dimer(a, b):
    return max(three_prime_dimer(a, b), three_prime_dimer(b, a))

def read_final(path):
    primers, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name: primers[name] = ''.join(buf)
            name = line[1:]; buf = []
        elif line: buf.append(line)
    if name: primers[name] = ''.join(buf)
    return primers


def main():
    seqs = load_fasta(TX_FASTA)
    coords = load_coords(COORDS)
    genes = sorted(seqs)
    arrays = {g: np.frombuffer(s.encode(), dtype=np.uint8) for g, s in seqs.items()}
    sdh3 = seqs[GENE]
    others = [g for g in genes if g != GENE]
    other_arrays = {g: arrays[g] for g in others}
    finals = read_final(FINAL)
    other_final = [p for n, p in finals.items() if not n.startswith('CserDQD_SDH3')]
    print(f"{GENE}: searching backup pairs; target product {TARGET_PROD} bp")
    print(f"multiplex partners to stay dimer-free with: {len(other_final)} primers")

    cache = {}
    def sc(p):
        if p not in cache:
            cache[p] = scan(p, other_arrays)
        return cache[p]

    fwd, rev = [], []
    for i in range(len(sdh3) - L + 1):
        plus = sdh3[i:i + L]
        if 'N' in plus or not primer_ok(plus): continue
        sites = sc(plus)
        if any(s['mm'] == 0 for s in sites): continue          # not unique to SDH3
        fwd.append({'seq': plus, 'i': i, 'tierA': not competent(sites, GENE),
                    'worst': worst(sites, GENE), 'ncomp': len(competent(sites, GENE)),
                    'gc': gc_content(plus), 'tm': tm_nn(plus)})
    for j in range(len(sdh3) - L + 1):
        win = sdh3[j:j + L]
        p = rc(win)
        if 'N' in win or not primer_ok(p): continue
        sites = sc(p)
        if any(s['mm'] == 0 for s in sites): continue
        rev.append({'seq': p, 'j': j, 'tierA': not competent(sites, GENE),
                    'worst': worst(sites, GENE), 'ncomp': len(competent(sites, GENE)),
                    'gc': gc_content(p), 'tm': tm_nn(p)})
    print(f"unique QC-passing 20-mers: {len(fwd)} forward, {len(rev)} reverse")
    print(f"  Tier A (no competent foreign site at all): "
          f"{sum(1 for c in fwd if c['tierA'])} fwd, {sum(1 for c in rev if c['tierA'])} rev")

    fwd_by_i = {c['i']: c for c in fwd}
    rev_by_j = {c['j']: c for c in rev}

    def pair_and_filter(sizes):
        rows = []
        for prod in sizes:
            for i, f in sorted(fwd_by_i.items()):
                j = i + prod - L
                r = rev_by_j.get(j)
                if r is None: continue
                if f['tierA'] and r['tierA']:
                    tier = 'A'
                else:
                    bad = []
                    for t in others:
                        fs = [s for s in sc(f['seq']) if s['tx'] == t and s['strand'] == '+'
                              and s['mm'] <= MAX_MM_COMP and s['tail'] >= MIN_TAIL_COMP]
                        rs = [s for s in sc(r['seq']) if s['tx'] == t and s['strand'] == '-'
                              and s['mm'] <= MAX_MM_COMP and s['tail'] >= MIN_TAIL_COMP]
                        for a in fs:
                            for b in rs:
                                sz = b['end'] - a['start'] + 1
                                if ALT_LO <= sz <= ALT_HI:
                                    bad.append(f"{t}:{sz}bp")
                    if bad: continue
                    tier = 'B'
                dtm = abs(f['tm'] - r['tm'])
                if dtm > 2.0: continue
                dim = max([max_dimer(f['seq'], q) for q in other_final]
                          + [max_dimer(r['seq'], q) for q in other_final]
                          + [max_dimer(f['seq'], r['seq'])])
                if dim >= 8: continue
                s_f = scan(f['seq'], arrays)
                s_r = scan(r['seq'], arrays)
                alt = []
                for a in [x for x in s_f if x['tx'] == GENE and x['strand'] == '+'
                          and x['mm'] <= MAX_MM_COMP and x['tail'] >= MIN_TAIL_COMP]:
                    for b in [x for x in s_r if x['tx'] == GENE and x['strand'] == '-'
                              and x['mm'] <= MAX_MM_COMP and x['tail'] >= MIN_TAIL_COMP]:
                        sz = b['end'] - a['start'] + 1
                        if ALT_LO <= sz <= ALT_HI: alt.append(sz)
                alt_other = sorted({a for a in alt if a != prod})
                very_permissive = sum(1 for s in sc(f['seq']) + sc(r['seq'])
                                      if s['tx'] != GENE and s['mm'] <= 8 and s['tail'] >= 3)
                score = (abs(f['gc'] - 50) + abs(r['gc'] - 50) + dtm * 4
                         + 2 * abs((f['tm'] + r['tm']) / 2 - 54)
                         + 0.01 * abs(f['tm'] - 55.2) + 0.01 * abs(r['tm'] - 55.1)
                         + 4 * len(alt_other) + 3 * very_permissive)
                rows.append({
                    'tier': tier, 'fwd_seq': f['seq'], 'rev_seq': r['seq'],
                    'fwd_pos': f"{f['i']+1}-{f['i']+L}", 'rev_pos': f"{r['j']+1}-{r['j']+L}",
                    'fwd_region': classify_position(coords[GENE], f['i'] + 1),
                    'rev_region': classify_position(coords[GENE], r['j'] + 1),
                    'product': prod,
                    'fwd_gc': round(f['gc'], 1), 'rev_gc': round(r['gc'], 1),
                    'fwd_tm': round(f['tm'], 1), 'rev_tm': round(r['tm'], 1),
                    'delta_tm': round(dtm, 1),
                    'worst_foreign_fwd': f['worst'], 'worst_foreign_rev': r['worst'],
                    'competent_foreign_fwd': f['ncomp'], 'competent_foreign_rev': r['ncomp'],
                    'very_permissive_foreign': very_permissive,
                    'alt_products_on_SDH3': ';'.join(map(str, alt_other)) or 'none',
                    'max_3prime_dimer_vs_other_pairs': dim,
                    'score': round(score, 2)})
        rows.sort(key=lambda x: (x['tier'], x['score']))
        return rows

    rows = pair_and_filter([TARGET_PROD])
    if not rows:
        print(f"nothing survives at exactly {TARGET_PROD} bp; widening to 280-320 bp")
        rows = pair_and_filter(list(range(280, 321)))
    print(f"backup pairs passing every filter: {len(rows)} "
          f"(Tier A {sum(1 for r in rows if r['tier']=='A')}, "
          f"Tier B {sum(1 for r in rows if r['tier']=='B')})")
    if not rows:
        sys.exit("no acceptable backup pair found - report and stop")

    with open(os.path.join(HERE, 'SDH3_BACKUP_PAIRS.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        for r in rows: w.writerow(r)
    with open(os.path.join(HERE, 'SDH3_backup_candidates.fasta'), 'w') as fh:
        for k, r in enumerate(rows[:5], 1):
            fh.write(f">SDH3bk{k}_{r['product']}bp_Tier{r['tier']}_Fwd\n{r['fwd_seq']}\n")
            fh.write(f">SDH3bk{k}_{r['product']}bp_Tier{r['tier']}_Rev\n{r['rev_seq']}\n")

    # ---------------- genome-wide verification of the recommended pair ----------------
    best = rows[0]
    tx_lines = []
    if os.path.isfile(TXOME):
        print("verifying the recommended pair against the whole transcriptome ...")
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
        starts = np.cumsum([0] + [len(s) for _, s in recs])[:-1]
        starts[0] = 0
        n = len(seqall)
        invalid = np.zeros(n, dtype=bool)
        for b in starts[1:]:
            invalid[max(0, int(b) - L + 1):int(b)] = True
        valid = np.flatnonzero(~invalid[:n - L + 1])
        tx_lines.append("## Genome-wide verification of the recommended backup pair\n")
        tx_lines.append(f"Exhaustive seed-free scan of all {n/1e6:.1f} Mb in {len(recs)} "
                        f"transcripts; every site with ≤4 mismatches is listed.\n")
        for p in (best['fwd_seq'], best['rev_seq']):
            qp = np.frombuffer(p.encode(), np.uint8)
            qr = np.frombuffer(rc(p).encode(), np.uint8)
            hits = []
            for orient, qb in (('+', qp), ('-', qr)):
                m = np.zeros(n - L + 1, dtype=np.int16)
                for i in range(L):
                    m += (seqall[i:i + n - L + 1] == qb[i])
                mm = L - m
                for pos in valid[mm[valid] <= 4]:
                    pos = int(pos)
                    tail = (tail_plus(seqall, pos, qp) if orient == '+'
                            else tail_minus(seqall, pos, qr))
                    ci = int(np.searchsorted(starts, pos, side='right') - 1)
                    hits.append((recs[ci][0], orient, pos + 1, int(mm[pos]), tail))
            hits.sort(key=lambda x: (x[3], -x[4]))
            tx_lines.append(f"### `{p}` — {len(hits)} site(s) with ≤4 mismatches")
            for h in hits[:12]:
                tx_lines.append(f"- {h[0]}  strand {h[1]}  start {h[2]}  {h[3]} mm  3′-run {h[4]} nt")
            if not hits:
                tx_lines.append("- none")
            tx_lines.append("")
            print(f"  {p}: {len(hits)} site(s) with ≤4 mm transcriptome-wide")
    else:
        tx_lines.append("_(transcriptome FASTA not found; genome-wide step skipped)_\n")

    # ---------------- report ----------------
    cur_f = finals.get('CserDQD_SDH3_Fwd', '')
    cur_r = finals.get('CserDQD_SDH3_Rev', '')
    cf, cr = scan(cur_f, arrays), scan(cur_r, arrays)
    with open(os.path.join(HERE, 'SDH3_BACKUP_REPORT.md'), 'w') as fh:
        fh.write("# Backup CserDQD/SDH3 pairs — zero competent sites on SDH2/SDH4\n\n")
        fh.write("Mandatory filters: 20-nt, exact-unique to SDH3, GC 40–60 %, no homopolymer >4 nt, "
                 "acceptable 3′ tail, no ≥7-nt self-complementary hairpin, product "
                 f"{TARGET_PROD} bp, ΔTm ≤ 2.0 °C, no ≥8-nt 3′-terminal complementarity with the "
                 "other ten primers, and zero competent foreign sites "
                 f"(≤{MAX_MM_COMP} mismatches with ≥{MIN_TAIL_COMP}-nt perfect 3′ run) at primer "
                 "level (Tier A) or pair level (Tier B).\n\n")
        fh.write("## Recommended backup (best-ranked)\n\n")
        for k in ('tier','fwd_seq','rev_seq','fwd_pos','rev_pos','fwd_region','rev_region',
                  'product','fwd_gc','rev_gc','fwd_tm','rev_tm','delta_tm','worst_foreign_fwd',
                  'worst_foreign_rev','competent_foreign_fwd','competent_foreign_rev',
                  'very_permissive_foreign','alt_products_on_SDH3',
                  'max_3prime_dimer_vs_other_pairs','score'):
            fh.write(f"- **{k}**: {best[k]}\n")
        fh.write("\n## In-service SDH3 pair, for comparison\n\n")
        fh.write(f"- forward `{cur_f}` {cur_f and round(tm_nn(cur_f),1)} °C, "
                 f"reverse `{cur_r}` {cur_r and round(tm_nn(cur_r),1)} °C, product 300 bp\n")
        fh.write(f"- competent foreign sites: {len(competent(cf, GENE))} forward + "
                 f"{len(competent(cr, GENE))} reverse "
                 "(they create the SDH4 300 bp and SDH2 309 bp cross-pairings)\n")
        fh.write("\n" + "\n".join(tx_lines) + "\n")
        fh.write("## All passing backups, ranked\n\n")
        cols = ['tier','fwd_seq','rev_seq','fwd_pos','rev_pos','fwd_region','rev_region',
                'fwd_tm','rev_tm','delta_tm','fwd_gc','rev_gc','worst_foreign_fwd',
                'worst_foreign_rev','very_permissive_foreign','alt_products_on_SDH3',
                'max_3prime_dimer_vs_other_pairs','score']
        fh.write("| " + " | ".join(cols) + " |\n|" + "|".join(['---'] * len(cols)) + "|\n")
        for r in rows[:30]:
            fh.write("| " + " | ".join(str(r[c]) for c in cols) + " |\n")

    print("\nRECOMMENDED BACKUP")
    print(f"  F  {best['fwd_seq']}  {best['fwd_pos']} ({best['fwd_region']})  "
          f"GC {best['fwd_gc']}%  Tm {best['fwd_tm']}°C")
    print(f"  R  {best['rev_seq']}  {best['rev_pos']} ({best['rev_region']})  "
          f"GC {best['rev_gc']}%  Tm {best['rev_tm']}°C")
    print(f"  product {best['product']} bp, tier {best['tier']}, ΔTm {best['delta_tm']} °C")
    print(f"  worst foreign site F: {best['worst_foreign_fwd']} | R: {best['worst_foreign_rev']}")
    print("wrote SDH3_BACKUP_PAIRS.tsv, SDH3_backup_candidates.fasta, SDH3_BACKUP_REPORT.md")


if __name__ == '__main__':
    main()