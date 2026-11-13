#!/usr/bin/env python3
"""
/conserved_blocks.py – Family-level conservation survey using MAFFT alignment.

For each family:
  1. MAFFT-align the protein sequences.
  2. Build codon-level CDS alignment from the protein alignment.
  3. Find conserved 20-mer windows in the CDS (100% identity across all members, no gaps).
  4. Also find windows with ≤2 mismatches and perfect 3′-anchor (last 8 nt).
  5. Report: number of conserved candidates, conserved territory, reference windows.

Output: data/processed/family_conservation.tsv, per-family stats.
"""
import os, sys, subprocess, re, math
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

print("Imported paths.py from:", P.__file__)
print("Project root:", P.ROOT)
print("RAW directory:", P.RAW)
print("Expected CDS:", P.raw('family_cds.fasta'))
print("CDS exists:", os.path.exists(P.raw('family_cds.fasta')))

def rc(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp[c] for c in reversed(s))

def gc_content(s):
    return (s.count('G')+s.count('C'))/len(s)*100

CODON = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GCC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C',
    'TGT':'C','TGA':'*','TGG':'W',
}

def load_fasta(path):
    seqs = OrderedDict()
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.upper())
        if name: seqs[name] = ''.join(buf)
    return seqs

def write_fasta(seqs, path):
    with open(path,'w') as f:
        for k,v in seqs.items():
            f.write(f'>{k}\n')
            for i in range(0,len(v),60): f.write(v[i:i+60]+'\n')

def run_mafft(seqs, outpath):
    """Run MAFFT on an OrderedDict {label: sequence}, return list of aligned labels."""
    tmppath = outpath + '.in.fasta'
    write_fasta(seqs, tmppath)
    in_order = list(seqs.keys())
    cmd = ['mafft', '--auto', tmppath]
    with open(outpath, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    os.remove(tmppath)
    aligned = load_fasta(outpath)
    return aligned, list(aligned.keys())

def build_codon_alignment(prot_align, protein_seqs, cds_seqs):
    """Given a protein alignment and the original unaligned CDS+protein dicts,
    build a codon-level nucleotide alignment by inserting 3-nt gaps where the protein has gaps."""
    result = OrderedDict()
    if not prot_align:
        return result
    # Get aligned sequences
    align_labels = list(prot_align.keys())
    if not align_labels:
        return result
    # Build a codon map for each sequence
    for label in align_labels:
        p_aln = prot_align[label].upper()
        raw_prot = protein_seqs.get(label, '')
        raw_cds = cds_seqs.get(label, '')
        # Map aligned protein columns to CDS positions
        cds_aln = []
        pi = 0  # position in raw_prot (ungapped)
        for i, aa in enumerate(p_aln):
            if aa == '-':
                cds_aln.extend(['-','-','-'])
            else:
                if pi*3 + 2 < len(raw_cds):
                    cds_aln.append(raw_cds[pi*3])
                    cds_aln.append(raw_cds[pi*3+1])
                    cds_aln.append(raw_cds[pi*3+2])
                else:
                    cds_aln.extend(['N','N','N'])
                pi += 1
        result[label] = ''.join(cds_aln)
    return result

def find_conserved_windows(cds_align, K=20):
    """
    Find windows of K consecutive columns in the CDS alignment that are:
      - gap-free in all sequences
      - identical (100% match) across all sequences
    Returns list of (start_column_1based, consensus_sequence).
    """
    if not cds_align:
        return []
    seqs = list(cds_align.values())
    n = len(seqs[0])
    K = min(K, n)
    windows = []
    for start in range(n - K + 1):
        col_block = [seq[start:start+K] for seq in seqs]
        cons = col_block[0]
        # Check no gaps and all identical
        ok = True
        for s in col_block:
            if '-' in s or 'N' in s:
                ok = False
                break
            if s != cons:
                ok = False
                break
        if ok:
            windows.append((start+1, cons))
    return windows

def find_permissive_windows(cds_align, reference_label, K=20, max_mm=2):
    """
    Find windows where reference windows appear in ≥N-1 members with ≤max_mm mismatches
    and perfect 3′ anchor (last 8 nt match = 0 mismatches).
    """
    if not cds_align or reference_label not in cds_align:
        return []
    ref_seq = cds_align[reference_label]
    ncols = len(ref_seq)
    other_seqs = [cds_align[l] for l in cds_align if l != reference_label]
    windows = []
    for start in range(ncols - K + 1):
        ref_win = ref_seq[start:start+K]
        if '-' in ref_win or 'N' in ref_win:
            continue
        cons = ref_win
        anchor = ref_win[-8:]
        passed = 0
        for o in other_seqs:
            o_win = o[start:start+K]
            if '-' in o_win or 'N' in o_win:
                continue
            mm = sum(1 for x,y in zip(ref_win, o_win) if x!=y)
            if mm <= max_mm and o_win[-8:] == anchor:
                passed += 1
        if passed >= len(other_seqs):  # all others pass
            windows.append((start+1, cons, passed, len(other_seqs)))
    return windows

FAMILIES = ['ACT','BTUB','EF1a','GAPDH','H3','eIF5A','polyUb']
FAMILY_DESCRIPTION = {
    'ACT': 'Actin', 'BTUB': 'Beta-tubulin', 'EF1a': 'Elongation factor 1-alpha',
    'GAPDH': 'Glyceraldehyde-3-phosphate dehydrogenase', 'H3': 'Histone H3',
    'eIF5A': 'Eukaryotic translation initiation factor 5A', 'polyUb': 'Polyubiquitin'
}

def main():
    P.require(P.raw('family_cds.fasta'), 'Run build_targets.py first')
    P.require(P.raw('family_protein.fasta'), 'Run build_targets.py first')
    P.require(P.raw('family_targets.tsv'), 'Run build_targets.py first')

    cds_all = load_fasta(P.raw('family_cds.fasta'))
    prot_all = load_fasta(P.raw('family_protein.fasta'))

    print("=== Step 2: Family conservation survey with MAFFT alignment ===\n")

    summary_rows = []

    for fam in FAMILIES:
        print(f"\n--- {fam} ({FAMILY_DESCRIPTION[fam]}) ---")
        # Filter CDS/protein sequences for this family
        fam_labels = sorted([l for l in cds_all if l.startswith(fam+'_') or l.startswith(fam+' ')])
        if len(fam_labels) < 2:
            print(f"  WARNING: Only {len(fam_labels)} member, skipping pairwise")
            continue

        fam_prot = OrderedDict((l, prot_all[l]) for l in fam_labels)
        fam_cds = OrderedDict((l, cds_all[l]) for l in fam_labels)

        print(f"  Members: {len(fam_labels)}")
        for l in fam_labels:
            print(f"    {l}: CDS={len(cds_all[l])} bp, prot={len(prot_all[l])} aa")

        # Align proteins
        aln_out = P.proc(f'{fam}_prot_aln.fasta')
        prot_aln, aln_labels = run_mafft(fam_prot, aln_out)
        print(f"  Protein alignment: {len(prot_aln)} seqs, length {len(list(prot_aln.values())[0]) if prot_aln else 0} columns")

        # Build codon alignment
        cds_aln = build_codon_alignment(prot_aln, fam_prot, fam_cds)
        if not cds_aln:
            print(f"  SKIP: codon alignment failed")
            continue
        aln_len = len(list(cds_aln.values())[0])
        print(f"  Codon alignment: {aln_len} columns")

        # Save per-family codon alignment
        cds_aln_path = P.proc(f'{fam}_codon_aln.fasta')
        with open(cds_aln_path, 'w') as f:
            for k, v in cds_aln.items():
                f.write(f'>{k}\n')
                for i in range(0, len(v), 60):
                    f.write(v[i:i+60] + '\n')

        # Find conserved windows (100% exact)
        exact_wins = find_conserved_windows(cds_aln, K=20)
        print(f"  Exact conserved 20-mer windows: {len(exact_wins)}")

        # Find permissive windows
        ref_label = fam_labels[0]  # first member as reference
        perm_wins = find_permissive_windows(cds_aln, ref_label, max_mm=2)
        print(f"  Permissive 20-mer windows (≤2 mm, perfect 3′ anchor): {len(perm_wins)}")

        # Report first few exact windows
        for col, seq in exact_wins[:3]:
            print(f"    column {col:5d}: {seq} GC={gc_content(seq):.1f}%")
        if len(exact_wins) > 3:
            print(f"    ... and {len(exact_wins)-3} more")

        summary_rows.append((fam, len(fam_labels), aln_len, len(exact_wins), len(perm_wins)))

    # Write summary
    print("\n=== Conservation Summary ===")
    print(f"{'Family':8s} {'Members':7s} {'AlnLen':6s} {'Exact':5s} {'Permissive':10s}")
    print("-"*40)
    for fam, n_mem, aln_len, n_exact, n_perm in summary_rows:
        print(f"{fam:8s} {n_mem:3d}   {aln_len:5d}   {n_exact:3d}    {n_perm:4d}")

    # Also write to TSV
    with open(P.proc('family_conservation.tsv'), 'w') as f:
        f.write("Family\tMembers\tAlignedColumns\tExact20mers\tPermissive20mers\n")
        for r in summary_rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print(f"\nConservation summary written to {P.proc('family_conservation.tsv')}")

if __name__ == '__main__':
    main()
