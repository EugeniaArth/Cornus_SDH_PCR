#!/usr/bin/env python3
"""
Design one family‑level primer pair per family.
Uses real transcript sequences for hamming matching (tolerates codon indels).
"""
import os, sys, math
from collections import OrderedDict


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

def rc(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp[c] for c in reversed(s))

def gc_content(s):
    return (s.count('G')+s.count('C'))/len(s)*100

_NN = {
    'AA':(-7.6,-21.7),'TT':(-7.6,-21.7),'AT':(-7.2,-20.4),'TA':(-7.2,-21.3),
    'CA':(-8.5,-22.7),'TG':(-8.5,-22.7),'GT':(-8.4,-22.4),'AC':(-8.4,-22.4),
    'CT':(-7.8,-21.0),'AG':(-7.8,-21.0),'GA':(-8.2,-22.2),'TC':(-8.2,-22.2),
    'CG':(-10.6,-27.1),'GC':(-9.8,-24.4),'GG':(-8.0,-19.9),'CC':(-8.0,-19.9),
}

def tm_nn(s, dna_conc=250e-9, na_conc=50e-3):
    dH=dS=0.0
    for i in range(len(s)-1):
        step=s[i:i+2]
        if step not in _NN: return None
        h,ss=_NN[step]; dH+=h; dS+=ss
    for eb in (s[0],s[-1]):
        if eb in 'GC': dH+=0.1; dS+=-2.8
        elif eb in 'AT': dH+=2.3; dS+=4.1
        else: return None
    dS+=0.368*(len(s)-1)*math.log(na_conc)
    return dH*1000.0/(dS+1.987*math.log(dna_conc/4.0))-273.15

def homopolymer(s, n=4):
    return any(b*n in s for b in 'ATGC')

def three_prime_ok(s):
    tail=s[-5:]
    if tail.endswith('TTTT') or tail.endswith('AAAA'): return False
    if tail.count('G')>=4 or tail.count('C')>=4: return False
    return True

def hairpin(s, k=7):
    r=rc(s)
    return any(s[i:i+k] in r for i in range(len(s)-k+1))

def primer_ok(s):
    return (40 <= gc_content(s) <= 60 and not homopolymer(s)
            and three_prime_ok(s) and not hairpin(s))

def load_fasta(path):
    seqs=OrderedDict()
    name,buf=None,[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if name: seqs[name]=''.join(buf)
                name,buf=line[1:].split()[0],[]
            else: buf.append(line.upper())
        if name: seqs[name]=''.join(buf)
    return seqs

def hamming(a,b):
    return sum(1 for x,y in zip(a,b) if x!=y)

def best_match(primer, tseq, K=20):
    """Return best hamming distance for primer in tseq."""
    best=20
    for wi in range(len(tseq)-K+1):
        mm=hamming(primer,tseq[wi:wi+K])
        if mm<best: best=mm
    return best

FAMILIES = ['ACT','BTUB','EF1a','GAPDH','H3','eIF5A','polyUb']
FAMILY_DESC = {
    'ACT':'Actin','BTUB':'Beta-tubulin','EF1a':'Elongation factor 1-alpha',
    'GAPDH':'Glyceraldehyde-3-phosphate dehydrogenase','H3':'Histone H3',
    'eIF5A':'Eukaryotic translation initiation factor 5A','polyUb':'Polyubiquitin'
}
FAMILY_COVERAGE_RULES = {
    'ACT': {'min_cov': 2, 'max_mm': 0},
    'BTUB': {'min_cov': 2, 'max_mm': 0},
    'EF1a': {'min_cov': 8, 'max_mm': 5},
    'GAPDH': {'min_cov': 2, 'max_mm': 0},
    'H3': {'min_cov': 6, 'max_mm': 3},
    'eIF5A': {'min_cov': 2, 'max_mm': 0},
    'polyUb': {'min_cov': 2, 'max_mm': 2},
}

def design_for_family(fam, cds_all, tx_seqs):
    print(f"\n{'='*50}")
    print(f"Designing family-level pair for {fam} ({FAMILY_DESC[fam]})")
    aln_path = P.proc(f'{fam}_codon_aln.fasta')
    if not os.path.exists(aln_path):
        print(f"  SKIP: no codon alignment at {aln_path}")
        return None
    align = load_fasta(aln_path)
    aln_labels = sorted([l for l in align if l.strip()])
    if len(aln_labels) < 2:
        print(f"  SKIP: alignment has <2 members")
        return None
    ref_label = max(aln_labels, key=lambda l: len(align[l]))
    ref = align[ref_label]
    n_cols = len(ref)
    rule = FAMILY_COVERAGE_RULES[fam]
    K = 20; MIN_P = 70; MAX_P = 250
    print(f"  Members: {len(aln_labels)}, ref={ref_label}, cols={n_cols}, rule={rule}")

    # Real transcript sequences
    real_seqs = {}
    for l in aln_labels:
        if l in tx_seqs: real_seqs[l] = tx_seqs[l]
        elif l in cds_all: real_seqs[l] = cds_all[l]

    # Survey forward candidates
    fwd_candidates = []
    for start in range(n_cols - K + 1):
        ref_win = ref[start:start+K]
        if '-' in ref_win or 'N' in ref_win: continue
        if not primer_ok(ref_win): continue
        fwd_mms = [best_match(ref_win, real_seqs.get(l,''), K) for l in aln_labels]
        cov = sum(1 for m in fwd_mms if m <= rule['max_mm'])
        if cov >= rule['min_cov']:
            fwd_candidates.append({'pos':start+1, 'seq':ref_win, 'cov':cov,
                                   'mms':fwd_mms})
    print(f"  Forward candidates: {len(fwd_candidates)}")
    if not fwd_candidates:
        # Relax primer_ok to GC 35-65, no homopolymer 5
        for start in range(n_cols - K + 1):
            ref_win = ref[start:start+K]
            if '-' in ref_win or 'N' in ref_win or homopolymer(ref_win,5): continue
            gc=gc_content(ref_win)
            if gc<35 or gc>65: continue
            fwd_mms = [best_match(ref_win, real_seqs.get(l,''), K) for l in aln_labels]
            cov = sum(1 for m in fwd_mms if m <= rule['max_mm'])
            if cov >= max(rule['min_cov']-1, 1):
                fwd_candidates.append({'pos':start+1, 'seq':ref_win, 'cov':cov, 'mms':fwd_mms})
        print(f"    After relax: {len(fwd_candidates)}")
    if not fwd_candidates: print("  FAIL: no forward candidate"); return None

    # Search pairs
    pairs = []
    for fw in fwd_candidates:
        f_pos = fw['pos']
        for r_start in range(f_pos + MIN_P - K, min(f_pos + MAX_P, n_cols - K + 1)):
            ref_win = ref[r_start:r_start+K]
            if '-' in ref_win or 'N' in ref_win: continue
            rev_seq_rc = rc(ref_win)
            if not primer_ok(rev_seq_rc):
                gc_r = gc_content(rev_seq_rc)
                if gc_r < 35 or gc_r > 65 or homopolymer(rev_seq_rc, 5): continue

            # Compute coverage on real sequences
            fwd_mms = [best_match(fw['seq'], real_seqs.get(l,''), K) for l in aln_labels]
            rev_mms = [best_match(rev_seq_rc, rc(real_seqs.get(l,'')), K) for l in aln_labels]
            fwd_cov = sum(1 for m in fwd_mms if m <= rule['max_mm'])
            rev_cov = sum(1 for m in rev_mms if m <= rule['max_mm'])

            if fwd_cov < rule['min_cov']: continue
            if rev_cov < max(1, rule['min_cov'] - 1): continue  # allow 1 less

            # Check that both F members and R members don't exceed max_mm+2 (indel tolerance)
            if any(m > rule['max_mm'] + 2 for m in rev_mms): continue
            if any(m > rule['max_mm'] for m in fwd_mms): continue

            tm_f = tm_nn(fw['seq'])
            tm_r = tm_nn(rev_seq_rc)
            if tm_f is None or tm_r is None or abs(tm_f - tm_r) > 6: continue

            # Check dimer (5-nt 3' complementarity)
            r_rc = rc(rev_seq_rc)
            def has_dimer(f, r, k=5):
                return any(r[i:i+k] in f for i in range(len(r)-k+1)) or \
                       any(f[i:i+k] in r_rc for i in range(len(f)-k+1))
            if has_dimer(fw['seq'], rev_seq_rc): continue

            # Real product on reference
            ref_t = real_seqs.get(ref_label, '')
            fi = ref_t.find(fw['seq'])
            ri = ref_t.rfind(rc(rev_seq_rc))
            prod = (ri+K-fi) if fi>=0 and ri>fi else (r_start+K-1)-f_pos+1

            score = (10 - fwd_cov/len(aln_labels)*10) + abs(tm_f-60) + abs(tm_r-60)
            score += abs(gc_content(fw['seq'])-50)/2 + abs(gc_content(rev_seq_rc)-50)/2 + abs(prod-130)*0.02

            pairs.append({'fwd_seq':fw['seq'],'rev_seq':rev_seq_rc,'fwd_pos':f_pos,'rev_pos':r_start+K-1,
                          'fwd_tm':tm_f,'rev_tm':tm_r,'fwd_gc':gc_content(fw['seq']),'rev_gc':gc_content(rev_seq_rc),
                          'product':prod,'fwd_cov':fwd_cov,'rev_cov':rev_cov,'total':len(aln_labels),'score':score})

    print(f"  Candidate pairs: {len(pairs)}")
    if not pairs: print(f"  FAIL: no valid pair found"); return None
    pairs.sort(key=lambda x: x['score'])
    best = pairs[0]
    print(f"  ** {best['fwd_seq']} / {best['rev_seq']}  Tm={best['fwd_tm']:.1f}/{best['rev_tm']:.1f}  "
          f"GC={best['fwd_gc']:.0f}/{best['rev_gc']:.0f}  cov={best['fwd_cov']}/{best['rev_cov']}/{best['total']}  "
          f"prod={best['product']}bp  score={best['score']:.1f}")
    return best

def main():
    print("=== Step 3: Family-level primer pairs ===\n")
    cds_all = load_fasta(P.raw('family_cds.fasta'))
    tx_seqs = load_fasta(P.raw('family_transcripts.fasta'))
    results = {}
    for fam in FAMILIES:
        results[fam] = design_for_family(fam, cds_all, tx_seqs)
    out = P.results('FAMILY_primer_pairs.tsv')
    with open(out, 'w') as f:
        f.write("Family\tForward\tReverse\tProduct\tFwdTm\tRevTm\tFwdGC\tRevGC\tFwdCov\tRevCov\tTotal\n")
        for fam in FAMILIES:
            b=results[fam]
            if b: f.write(f"{fam}\t{b['fwd_seq']}\t{b['rev_seq']}\t{b['product']}\t"
                          f"{b['fwd_tm']:.1f}\t{b['rev_tm']:.1f}\t{b['fwd_gc']:.1f}\t{b['rev_gc']:.1f}\t"
                          f"{b['fwd_cov']}\t{b['rev_cov']}\t{b['total']}\n")
    fa = P.results('FAMILY_primers.fasta')
    with open(fa, 'w') as f:
        for fam in FAMILIES:
            b=results[fam]
            if b: f.write(f">{fam}_F\n{b['fwd_seq']}\n>{fam}_R\n{b['rev_seq']}\n")
    print(f"\nWritten: {out}, {fa}")

if __name__ == '__main__':
    main()
