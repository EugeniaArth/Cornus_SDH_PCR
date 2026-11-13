#!/usr/bin/env python3
"""
step1_inputs/build_targets.py — Build target CDS / protein sets from filtered family FASTA files.

Reads family_candidates_filtered/<FAMILY>_filtered.fasta (produced by step0), finds the longest
clean ORF to determine sense orientation and CDS, writes:

  data/raw/family_transcripts.fasta  — sense-oriented full transcripts
  data/raw/family_cds.fasta          — CDS nucleotide sequences
  data/raw/family_protein.fasta      — translated proteins
  data/raw/family_targets.tsv        — family × contig × ORF coordinates

No transcriptome FASTA or GFF needed — sequences are already in the filtered files from step0.
"""
import os, sys, re, bisect
from collections import OrderedDict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

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
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C',
    'TGT':'C','TGA':'*','TGG':'W',
}

def rc(s):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    return ''.join(comp[c] for c in reversed(s))

def longest_clean_orf(seq):
    """Return (1-based_start, end, aa_length) of longest clean ATG→stop ORF."""
    n = len(seq)
    starts = [[] for _ in range(3)]
    for i in range(n - 2):
        if seq[i:i+3] == 'ATG':
            starts[i%3].append(i)
    stops = [[] for _ in range(3)]
    for i in range(n - 2):
        if CODON.get(seq[i:i+3], '') == '*':
            stops[i%3].append(i)
    best = (0, 0, 0)
    for fr in range(3):
        st = sorted(starts[fr])
        sp = sorted(stops[fr])
        for a in st:
            j = bisect.bisect_left(sp, a)
            end = (sp[j] + 3) if j < len(sp) else (n - ((n - a) % 3))
            aa = (end - a) // 3
            if aa > best[2]:
                best = (a + 1, end, aa)
    return best

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

FAMILIES = ['ACT', 'BTUB', 'EF1a', 'GAPDH', 'H3', 'eIF5A', 'polyUb']
FAMILY_DESC = {
    'ACT': 'Actin', 'BTUB': 'Beta-tubulin', 'EF1a': 'Elongation factor 1-alpha',
    'GAPDH': 'Glyceraldehyde-3-phosphate dehydrogenase', 'H3': 'Histone H3',
    'eIF5A': 'Eukaryotic translation initiation factor 5A', 'polyUb': 'Polyubiquitin'
}

def main():
    print("=== Step 1: Build target sequences from filtered family FASTA files ===\n")

    # Load filtered FASTA files directly — no transcriptome or GFF needed
    all_contigs = {}
    for fam in FAMILIES:
        fp = os.path.join(P.FAMILY_CANDIDATES, f'{fam}_filtered.fasta')
        if not os.path.exists(fp):
            print(f"WARNING: {fp} not found — run step0 first")
            continue
        seqs = load_fasta(fp)
        for c, s in seqs.items():
            all_contigs[c] = (fam, s)
        print(f"  {fam}: {len(seqs)} contigs loaded")

    print(f"\nTotal: {len(all_contigs)} contigs")

    # Find ORF on both strands, pick the longest ≥50 aa
    records = []
    transcripts = OrderedDict()
    cdss = OrderedDict()
    proteins = OrderedDict()

    for contig, (fam, seq) in all_contigs.items():
        op = longest_clean_orf(seq)
        om = longest_clean_orf(rc(seq))

        if op[2] >= om[2] and op[2] >= 50:
            start, end, aa_len = op
            orf_seq = seq[start-1:end]
            transcript = seq
            strand = '+'
        elif om[2] >= 50:
            start, end, aa_len = om
            orf_seq = rc(seq[start-1:end])
            transcript = rc(seq)
            strand = '-'
        else:
            print(f"  WARNING: {contig} ({fam}) no ORF ≥50 aa, using full sequence as-is")
            start, end, aa_len = 1, len(seq), 0
            orf_seq = ''
            transcript = seq
            strand = '+'

        prot = ''.join(CODON.get(orf_seq[i:i+3], 'X') for i in range(0, len(orf_seq)-2, 3)) if orf_seq else ''

        contig_short = re.sub(r'^NODE_(\d+).*', r'NODE_\1', contig)
        label = f"{fam}_{contig_short}"
        transcripts[label] = transcript
        if orf_seq: cdss[label] = orf_seq
        if prot: proteins[label] = prot

        records.append(dict(family=fam, contig=contig, label=label,
                            strand=strand, ORF_start=start, ORF_end=end,
                            ORF_aa=aa_len, length=len(seq)))
        print(f"  {label}: {len(seq)} nt, ORF={start}-{end} ({aa_len} aa) {strand}")

    # Write outputs
    def write_fasta(d, path):
        with open(path, 'w') as f:
            for k, v in d.items():
                f.write(f'>{k}\n')
                for i in range(0, len(v), 60):
                    f.write(v[i:i+60] + '\n')
        print(f"  → {path} ({len(d)} seqs)")

    write_fasta(transcripts, P.raw('family_transcripts.fasta'))
    write_fasta(cdss, P.raw('family_cds.fasta'))
    write_fasta(proteins, P.raw('family_protein.fasta'))

    with open(P.raw('family_targets.tsv'), 'w') as f:
        f.write("Family\tLabel\tContig\tLength\tStrand\tORF_start\tORF_end\tORF_aa\n")
        for r in records:
            f.write(f"{r['family']}\t{r['label']}\t{r['contig']}\t{r['length']}\t"
                    f"{r['strand']}\t{r['ORF_start']}\t{r['ORF_end']}\t{r['ORF_aa']}\n")
    print(f"  → {P.raw('family_targets.tsv')}")

    counts = Counter(r['family'] for r in records)
    print("\nFinal target set:")
    for fam in FAMILIES:
        print(f"  {fam} ({FAMILY_DESC[fam]}): {counts.get(fam, 0)} contigs")
    print(f"Total: {len(records)} contigs")

if __name__ == '__main__':
    main()
