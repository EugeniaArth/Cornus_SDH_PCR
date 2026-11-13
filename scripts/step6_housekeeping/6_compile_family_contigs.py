#!/usr/bin/env python3
"""
Compile candidate contigs for each reference‑gene family.
Inputs:
  - BLAST intersection hits (high‑confidence)
  - Keyword‑search contigs (broad)
  - HMM‑validated contigs (after hmmscan, optional)
Output: per‑family FASTA files and a consolidated TSV.
"""
import os, sys, re
from collections import defaultdict

# Paths
HERE = os.path.dirname(os.path.abspath(__file__))
BASEDIR =  os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping')
TRANSCRIPTOME  = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'external', 'transcriptome.fasta')
OUTDIR = os.path.join(BASEDIR, 'family_candidates')


def load_blast_hits():
    """Load contigs from BLAST intersection (high‑confidence)."""
    blast_file = os.path.join(BASEDIR, 'agrostemma_blast', 'blast_hits_summary.tsv')
    families = defaultdict(set)
    if os.path.exists(blast_file):
        with open(blast_file, 'r') as f:
            next(f)  # header
            for line in f:
                if line.strip():
                    fam, contig, *_ = line.strip().split('\t')
                    families[fam].add(contig)
    return families

def load_keyword_contigs():
    """Load contigs from keyword search."""
    families = defaultdict(set)
    keyword_dir = os.path.join(BASEDIR, 'keyword_search')
    for fname in os.listdir(keyword_dir):
        if fname.endswith('_contigs.txt'):
            fam = fname.replace('_contigs.txt', '')
            with open(os.path.join(keyword_dir, fname), 'r') as f:
                for line in f:
                    contig = line.strip()
                    if contig:
                        families[fam].add(contig)
    return families

def load_hmm_validated():
    """Placeholder: will be filled after hmmscan."""
    return defaultdict(set)

def extract_sequences(contig_ids, out_fasta):
    """Extract sequences for given contig IDs from transcriptome FASTA."""
    # Simple linear scan (48 MB file is fine)
    wanted = set(contig_ids)
    found = set()
    with open(TRANSCRIPTOME, 'r') as fin, open(out_fasta, 'w') as fout:
        write = False
        seq = []
        for line in fin:
            if line.startswith('>'):
                # Save previous sequence
                if write and seq:
                    fout.write(''.join(seq))
                # Parse header
                header = line[1:].strip()
                # Extract contig ID (remove possible description)
                contig = header.split()[0]
                if contig in wanted:
                    write = True
                    fout.write(f'>{contig}\n')
                    found.add(contig)
                    seq = []
                else:
                    write = False
            elif write:
                seq.append(line)
        # Last sequence
        if write and seq:
            fout.write(''.join(seq))
    print(f"  Extracted {len(found)}/{len(wanted)} sequences to {out_fasta}")
    return found

def main():
    print("Compiling candidate contigs for each family...")

    blast = load_blast_hits()
    keyword = load_keyword_contigs()

    # Define families of interest (short names)
    families_def = {
        'ACT': 'Actin',
        'BTUB': 'Beta‑tubulin',
        'EF1a': 'Elongation factor 1‑alpha',
        'GAPDH': 'Glyceraldehyde‑3‑phosphate dehydrogenase',
        'H3': 'Histone H3',
        'eIF5A': 'Eukaryotic translation initiation factor 5A',
        'polyUb': 'Polyubiquitin'
    }

    # For each family, decide which contig sets to use
    # For now, take union of BLAST and keyword (except for polyUb where we have our own list)
    all_contigs = defaultdict(set)

    for fam in families_def:
        if fam in ['ACT', 'BTUB', 'EF1a', 'GAPDH']:
            # Use BCAST hits (high confidence)
            all_contigs[fam] = blast.get(fam, set())
            if not all_contigs[fam]:
                # fallback to keyword
                all_contigs[fam] = keyword.get(fam, set())
        elif fam in ['H3', 'eIF5A']:
            # Use keyword (will be filtered by HMM later)
            all_contigs[fam] = keyword.get(fam, set())
        elif fam == 'polyUb':
            # Use our existing Ubiquitinx6 contigs
            polyub_file = 'hmm_search/final_targets.tsv'
            if os.path.exists(polyub_file):
                with open(polyub_file, 'r') as f:
                    for line in f:
                        if 'PolyUb' in line:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                contig = parts[1]
                                all_contigs['polyUb'].add(contig)
            else:
                print("WARNING: final_targets.tsv not found; using keyword UBQ contigs")
                all_contigs[fam] = keyword.get('UBQ', set())

    # Write summary TSV
    with open(os.path.join(OUTDIR, 'family_contigs_summary.tsv'), 'w') as f:
        f.write("Family\tDescription\tContigID\tSource\n")
        for fam, desc in families_def.items():
            for contig in sorted(all_contigs[fam]):
                source = 'BLAST' if contig in blast.get(fam, set()) else 'keyword'
                f.write(f"{fam}\t{desc}\t{contig}\t{source}\n")

    # Extract sequences for each family
    for fam in families_def:
        contigs = all_contigs[fam]
        if contigs:
            out_fasta = os.path.join(OUTDIR, f'{fam}_candidates.fasta')
            found = extract_sequences(contigs, out_fasta)
            missing = contigs - found
            if missing:
                print(f"    WARNING: {len(missing)} contigs not found in transcriptome")
                with open(os.path.join(OUTDIR, f'{fam}_missing.txt'), 'w') as f:
                    for c in sorted(missing):
                        f.write(c + '\n')
        else:
            print(f"  {fam}: no contigs collected")

    print(f"\nResults saved to {OUTDIR}/")

if __name__ == '__main__':
    main()
