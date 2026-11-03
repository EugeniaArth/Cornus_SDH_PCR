#!/usr/bin/env python3
"""
Final gene-specific primer pair design for six Cornus DQD/SDH transcripts.

Strategy:
  * Extract every 20-mer unique to a given transcript (not found in any other transcript).
  * Preferred pair: BOTH forward and reverse unique anchors, facing each other, 120-400 bp apart.
  * Fallback (needed for SDH1.1, which owns only a 38 nt unique window vs SDH1.2; see data/processed/unique_20mer_spans.tsv):
    one unique anchor primer + one shared partner primer nearby. Specificity is
    guaranteed by the unique anchor (only that transcript can prime extension).
  * Report amplicons spread over the 120-400 bp range so products can be told apart.
"""
import csv
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

def gc_content(s):
    return (s.count('G') + s.count('C')) / len(s) * 100

def tm_basic(s):
    # Marmur-Doty / salt-adjusted formula for 14-100 nt (Na+ ~50 mM)
    gc = s.count('G') + s.count('C')
    return 64.9 + 41 * (gc - 16.4) / len(s)

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
    """SantaLucia 1998 unified NN Tm in °C, for non-self-complementary oligos."""
    import math
    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/(mol*K)
    for i in range(len(s) - 1):
        step = s[i:i+2]
        if step not in _NN:
            return tm_basic(s)
        h, ss = _NN[step]
        dH += h
        dS += ss
    # Initiation (per terminal base pair)
    for end_base in (s[0], s[-1]):
        if end_base in 'GC':
            dH += 0.1
            dS += -2.8
        elif end_base in 'AT':
            dH += 2.3
            dS += 4.1
        else:
            return tm_basic(s)
    # Salt correction
    dS += 0.368 * (len(s) - 1) * math.log(na_conc)
    R = 1.987  # cal/(mol*K)
    tm_K = dH * 1000.0 / (dS + R * math.log(dna_conc / 4.0))
    return tm_K - 273.15

def homopolymer(s, n=4):
    return any(b * n in s for b in 'ATGC')

def rc(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp[c] for c in reversed(s))

def three_prime_ok(s):
    tail = s[-5:]
    if tail.endswith('TTTT') or tail.endswith('AAAA'):
        return False
    if tail.count('G') >= 4 or tail.count('C') >= 4:
        return False
    return True

def hairpin(s, k=7):
    # Require at least k complementary bases in inverted orientation
    r = rc(s)
    for i in range(len(s) - k + 1):
        if s[i:i+k] in r:
            return True
    return False

def primer_ok(s):
    return (40 <= gc_content(s) <= 60
            and not homopolymer(s)
            and three_prime_ok(s)
            and not hairpin(s))

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

def load_coords(path):
    coords = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if not row['Gene']: continue
            coords[row['Gene']] = {
                'acc': row['Accession'],
                'len': int(row['Length']),
                'five':  (int(row['5UTR_start']),  int(row['5UTR_end'])),
                'cds':   (int(row['CDS_start']),   int(row['CDS_end'])),
                'three': (int(row['3UTR_start']),  int(row['3UTR_end'])),
            }
    return coords

def classify_position(c, pos_1based):
    if c['five'][0] <= pos_1based <= c['five'][1]:
        return '5UTR'
    if c['cds'][0] <= pos_1based <= c['cds'][1]:
        return 'CDS'
    if c['three'][0] <= pos_1based <= c['three'][1]:
        return '3UTR'
    return 'other'

def build_candidates(seqs, gene, uniq, K=20,
                     min_prod=120, max_prod=500, prefer_prod=220, max_dtm=4):
    s = seqs[gene]
    n = len(s)
    uniq_set = {p for p, _ in uniq}
    kmer_at = {p: k for p, k in uniq}

    both_unique = []
    anchor_only = []

    # ---- Pass 1: both primers unique ----
    for i, k_i in uniq:
        if not primer_ok(k_i):
            continue
        for j in range(i + (min_prod - K), min(i + max_prod, n - K + 1)):
            if j not in uniq_set:
                continue
            rbind = kmer_at[j]
            if not primer_ok(rbind):
                continue
            product = (j + K) - i
            if not (min_prod <= product <= max_prod):
                continue
            dtm = abs(tm_nn(k_i) - tm_nn(rbind))
            if dtm > max_dtm:
                continue
            both_unique.append({
                'fwd': k_i,
                'rev': rc(rbind),
                'fwd_pos': i + 1,
                'rev_pos': j + 1,
                'product': product,
                'fwd_gc': gc_content(k_i),
                'rev_gc': gc_content(rbind),
                'fwd_tm': tm_nn(k_i),
                'rev_tm': tm_nn(rbind),
                'score': abs(gc_content(k_i) - 50) + abs(gc_content(rbind) - 50)
                         + dtm * 3 + abs(product - prefer_prod) * 0.05,
                'mode': 'both primers unique',
            })

    if both_unique:
        both_unique.sort(key=lambda x: x['score'])
        return both_unique, 'both_unique'

    # ---- Pass 2: one unique anchor, one shared partner ----
    for i, k_i in uniq:
        if not primer_ok(k_i):
            continue
        # Anchor on forward strand, partner downstream (becomes reverse primer)
        for j in range(i + (min_prod - K), min(i + max_prod, n - K + 1)):
            if j in uniq_set:
                continue  # already tried in pass 1
            partner = s[j:j+K]
            if not primer_ok(partner):
                continue
            product = (j + K) - i
            if not (min_prod <= product <= max_prod):
                continue
            dtm = abs(tm_nn(k_i) - tm_nn(partner))
            if dtm > max_dtm:
                continue
            anchor_only.append({
                'fwd': k_i,
                'rev': rc(partner),
                'fwd_pos': i + 1,
                'rev_pos': j + 1,
                'product': product,
                'fwd_gc': gc_content(k_i),
                'rev_gc': gc_content(partner),
                'fwd_tm': tm_nn(k_i),
                'rev_tm': tm_nn(partner),
                'score': abs(gc_content(k_i) - 50) + abs(gc_content(partner) - 50)
                         + dtm * 3 + abs(product - prefer_prod) * 0.05,
                'mode': 'unique F + shared R (specificity from unique 5\' anchor)',
            })
        # Anchor on forward strand, partner upstream (anchor becomes reverse primer)
        for j in range(max(0, i - (max_prod - K)), i - (min_prod - K)):
            if j in uniq_set:
                continue
            partner = s[j:j+K]
            if not primer_ok(partner):
                continue
            product = (i + K) - j
            if not (min_prod <= product <= max_prod):
                continue
            dtm = abs(tm_nn(k_i) - tm_nn(partner))
            if dtm > max_dtm:
                continue
            anchor_only.append({
                'fwd': partner,
                'rev': rc(k_i),
                'fwd_pos': j + 1,
                'rev_pos': i + 1,
                'product': product,
                'fwd_gc': gc_content(partner),
                'rev_gc': gc_content(k_i),
                'fwd_tm': tm_nn(partner),
                'rev_tm': tm_nn(k_i),
                'score': abs(gc_content(partner) - 50) + abs(gc_content(k_i) - 50)
                         + dtm * 3 + abs(product - prefer_prod) * 0.05,
                'mode': 'shared F + unique R (specificity from unique 3\' anchor)',
            })

    anchor_only.sort(key=lambda x: x['score'])
    return anchor_only, 'anchor_only'

def main():
    seqs = load_fasta(P.require(P.raw('Cornus_DQD_SDHs.fasta'), 'curated cDNA, report 01'))
    coords = load_coords(P.require(P.proc('genbank_coordinates.tsv'), 'from step1_inputs/extract_genbank_coords.py'))
    genes = sorted(coords.keys())
    print(f"Transcripts: {genes}\n")

    K = 20
    unique_kmers = {}
    for g in genes:
        others = [seqs[o] for o in genes if o != g]
        uniq = []
        for i in range(len(seqs[g]) - K + 1):
            km = seqs[g][i:i+K]
            if any(km in o for o in others):
                continue
            uniq.append((i, km))
        unique_kmers[g] = uniq
        print(f"{g}: {len(uniq)} unique {K}-mers")

    print()

    per_gene = {}
    per_gene_mode = {}
    for g in genes:
        cands, mode = build_candidates(seqs, g, unique_kmers[g], K=K)
        per_gene[g] = cands
        per_gene_mode[g] = mode
        print(f"{g}: {len(cands)} candidate pairs (mode={mode})")

    # Greedy selection: spread sizes AND avoid dimers
    # We will assign DIFFERENT target product sizes per gene so multiplex gel
    # separation works, then pick a candidate closest to that target.
    TARGETS = {
        'CserDQD/SDH1.1': 120,
        'CserDQD/SDH1.2': 180,
        'CserDQD/SDH2':   240,
        'CserDQD/SDH3':   300,
        'CserDQD/SDH4':   360,
        'CserDQD/SDH5':   420,
    }

    def has_dimer(p1, p2, k=7):
        """Return True if any pair among p1/p2 has >=k bp 3' complementarity."""
        seqs_list = [(p1['fwd'], 'F1'), (p1['rev'], 'R1'), (p2['fwd'], 'F2'), (p2['rev'], 'R2')]
        for i, (s1, n1) in enumerate(seqs_list):
            for s2, n2 in seqs_list[i+1:]:
                if s1[-k:] == rc(s2)[-k:]:
                    return True
                if s2[-k:] == rc(s1)[-k:]:
                    return True
        return False

    def clashes_with_chosen(c, chosen_list):
        for other in chosen_list:
            if has_dimer(c, other):
                return True
        return False

    chosen = {}
    chosen_list = []
    # Process genes in fixed order
    for g in genes:
        cands = per_gene[g]
        if not cands:
            chosen[g] = None
            continue
        target = TARGETS.get(g, 220)
        # Sort candidates by |product - target| + tiny score
        scored = sorted(cands, key=lambda c: (abs(c['product'] - target), c['score']))
        pick = None
        for c in scored[:400]:
            if not clashes_with_chosen(c, chosen_list):
                pick = c
                break
        if pick is None:
            pick = scored[0]
        chosen[g] = pick
        chosen_list.append(pick)

    # Report dimers if any remain
    all_pairs = [(g, chosen[g]) for g in genes if chosen[g]]
    dimer_notes = []
    for i, (g1, p1) in enumerate(all_pairs):
        for g2, p2 in all_pairs[i+1:]:
            if has_dimer(p1, p2):
                dimer_notes.append((g1, g2))

    # Output table
    print("\n" + "=" * 92)
    print("FINAL PRIMER PAIRS (de novo, Cornus-specific, unique vs the other 5 transcripts)")
    print("=" * 92)

    header = ['Gene', 'Accession', 'Dir', 'Sequence (5\'→3\')', 'Len',
              'Position', 'Region', 'GC%', 'Tm(°C)', 'Amplicon(bp)', 'Mode']
    rows = []
    for g in genes:
        p = chosen[g]
        c = coords[g]
        if not p:
            rows.append([g, c['acc'], '-', 'NO PAIR AVAILABLE', '-', '-', '-',
                         '-', '-', '-', '-'])
            continue
        freg = classify_position(c, p['fwd_pos'])
        rreg = classify_position(c, p['rev_pos'])
        rows.append([g, c['acc'], 'F', p['fwd'], len(p['fwd']),
                     f"{p['fwd_pos']}-{p['fwd_pos']+19}", freg,
                     f"{p['fwd_gc']:.1f}", f"{p['fwd_tm']:.1f}", p['product'], p['mode']])
        rows.append(['', '', 'R', p['rev'], len(p['rev']),
                     f"{p['rev_pos']}-{p['rev_pos']+19}", rreg,
                     f"{p['rev_gc']:.1f}", f"{p['rev_tm']:.1f}", '', ''])

    with open(P.results('FINAL_primer_pairs.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(header)
        for r in rows:
            w.writerow(r)

    for r in rows:
        print('\t'.join(str(x) for x in r))

    print("\n--- Amplicon sizes ---")
    sizes = sorted([(g, chosen[g]['product']) for g in genes if chosen[g]],
                   key=lambda x: x[1])
    for g, s in sizes:
        print(f"  {g}: {s} bp")
    if len(sizes) >= 2:
        diffs = [sizes[i+1][1] - sizes[i][1] for i in range(len(sizes)-1)]
        print(f"  Minimum spacing: {min(diffs)} bp")
        if min(diffs) >= 40:
            print("  ✓ Spacing OK for agarose-gel / multiplex discrimination.")
        else:
            print("  ⚠ Consider fragment analyzer or separate reactions for close pairs.")

    if dimer_notes:
        print("\n⚠ Possible inter-primer 3' complementarity detected:")
        for d in dimer_notes[:10]:
            print("   ", d)
    else:
        print("\n✓ No strong inter-primer 3' complementarity among chosen set.")

    # FASTA
    with open(P.results('FINAL_primers.fasta'), 'w') as f:
        for g in genes:
            p = chosen[g]
            if not p: continue
            name = g.replace('/', '_')
            f.write(f">{name}_Fwd\n{p['fwd']}\n")
            f.write(f">{name}_Rev\n{p['rev']}\n")

    print("\nSaved: FINAL_primer_pairs.tsv, FINAL_primers.fasta")

if __name__ == '__main__':
    main()
