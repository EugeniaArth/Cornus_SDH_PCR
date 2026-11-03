#!/usr/bin/env python3
import sys
import csv
import os
import re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P
from collections import defaultdict

def read_cds_coords(csv_file):
    """Return dict gene -> (cds_start, cds_end, transcript_length)"""
    coords = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row['Gene']
            if gene not in coords:
                try:
                    cds_start = int(row['CDS_start'])
                    cds_end = int(row['CDS_end'])
                    tx_len = int(row['Length'])
                    coords[gene] = (cds_start, cds_end, tx_len)
                except (ValueError, KeyError):
                    print(f"Warning: missing CDS coordinates for {gene}")
    return coords

def parse_genbank_sequence(gb_file):
    """Extract the nucleotide sequence from a GenBank file.
    Returns uppercase sequence string.
    """
    with open(gb_file, 'r') as f:
        content = f.read()
    # Find ORIGIN block
    origin_start = content.find('ORIGIN')
    if origin_start == -1:
        raise ValueError(f'No ORIGIN found in {gb_file}')
    # Find the terminating '//' after ORIGIN
    seq_text = content[origin_start:]
    # We'll parse lines after ORIGIN up to '//'
    lines = seq_text.split('\n')
    seq_lines = []
    for line in lines[1:]:          # skip the 'ORIGIN' line itself
        line = line.strip()
        if line == '//':
            break
        # remove leading numbers and spaces, keep only letters
        # typical line: '    181 ccttgtagaa ...'
        # keep only a-z and A-Z
        letters = re.sub(r'[^a-zA-Z]', '', line)
        seq_lines.append(letters.upper())
    return ''.join(seq_lines)

def load_cDNA_from_genbank(coord_file):
    """Load cDNA sequences from GenBank files.
    Reads coordinate TSV, for each accession loads the corresponding .gb file,
    extracts sequence, maps gene name.
    Returns dict gene -> sequence.
    """
    # mapping from accession to gene (and length for validation)
    acc_to_gene = {}
    with open(coord_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            acc = row['Accession']
            gene = row['Gene']
            length = int(row['Length'])
            acc_to_gene[acc] = (gene, length)

    seqs = {}
    for acc, (gene, expected_len) in acc_to_gene.items():
        gb_path = os.path.join(P.RAW, f'{acc}.gb')
        if not os.path.exists(gb_path):
            print(f'ERROR: missing GenBank file {gb_path}')
            continue
        seq = parse_genbank_sequence(gb_path)
        if len(seq) != expected_len:
            print(f'Warning: {gene} length mismatch: GB {len(seq)} vs coordinate table {expected_len}')
        seqs[gene] = seq
    return seqs

def main():
    csv_file = P.proc('genbank_coordinates.tsv')

    seqs = load_cDNA_from_genbank(csv_file)
    print(f"Loaded {len(seqs)} cDNA sequences from GenBank records")
    # Write combined cDNA FASTA for downstream steps (replaces DQD_SDH_possible_nt.fasta)
    with open(P.raw('Cornus_DQD_SDHs.fasta'), 'w') as f:
        for gene, seq in sorted(seqs.items()):
            f.write(f'>{gene}\n{seq}\n')
    print("Written combined cDNA to", P.raw('Cornus_DQD_SDHs.fasta'))

    coords = read_cds_coords(csv_file)
    print(f"Loaded CDS coordinates for {len(coords)} genes")

    # Check all genes present
    missing = [g for g in seqs if g not in coords]
    if missing:
        print(f"Warning: missing CDS coordinates for {missing}")
        # maybe skip
        for g in missing:
            seqs.pop(g, None)

    # Extract UTRs
    five_utrs = {}
    three_utrs = {}
    cds_seqs = {}
    for gene, seq in seqs.items():
        if gene not in coords:
            continue
        start, end, tx_len = coords[gene]
        # Validate
        if start < 1 or end > len(seq) or start > end:
            print(f"Invalid coordinates for {gene}: start={start}, end={end}, len={len(seq)}")
            continue
        # 5' UTR: bases 1 to start-1 (since start is 1-indexed)
        five_utr = seq[:start-1]
        # CDS: bases start to end inclusive
        cds = seq[start-1:end]
        # 3' UTR: bases end+1 to end of sequence
        three_utr = seq[end:]

        five_utrs[gene] = five_utr
        three_utrs[gene] = three_utr
        cds_seqs[gene] = cds

        print(f"{gene}: transcript {len(seq)} bp, CDS {len(cds)} bp ({start}-{end}), 5'UTR {len(five_utr)} bp, 3'UTR {len(three_utr)} bp")

    # Write UTRs to files
    with open(P.proc('cornus_5utr.fasta'), 'w') as f:
        for gene, utr in five_utrs.items():
            f.write(f">{gene}\n{utr}\n")
    with open(P.proc('cornus_3utr.fasta'), 'w') as f:
        for gene, utr in three_utrs.items():
            f.write(f">{gene}\n{utr}\n")
    with open(P.proc('cornus_cds.fasta'), 'w') as f:
        for gene, cds in cds_seqs.items():
            f.write(f">{gene}\n{cds}\n")

    print("\nWritten cornus_5utr.fasta, cornus_3utr.fasta, cornus_cds.fasta")

    # Generate summary table
    with open(P.proc('cornus_utr_summary.tsv'), 'w') as f:
        f.write("Gene\tTranscript_len\tCDS_len\t5UTR_len\t3UTR_len\t5UTR_seq\t3UTR_seq\n")
        for gene in sorted(seqs.keys()):
            if gene in coords:
                start, end, tx_len = coords[gene]
                f.write(f"{gene}\t{tx_len}\t{len(cds_seqs[gene])}\t{len(five_utrs[gene])}\t{len(three_utrs[gene])}\t{five_utrs[gene]}\t{three_utrs[gene]}\n")

    # Now look for unique k-mers in UTRs (optional)
    print("\nSearching for unique 20-mers in 5' UTRs...")
    k = 20
    five_kmers = {}
    for gene, utr in five_utrs.items():
        kmers = set()
        for i in range(len(utr) - k + 1):
            kmers.add(utr[i:i+k])
        five_kmers[gene] = kmers

    # Count occurrences across genes
    kmer_counts = defaultdict(int)
    for gene, kmers in five_kmers.items():
        for kmer in kmers:
            kmer_counts[kmer] += 1

    # Find unique kmers (appear only in one gene)
    unique_kmers = {}
    for gene, kmers in five_kmers.items():
        uniq = sorted(kmer for kmer in kmers if kmer_counts[kmer] == 1)  # sorted: reproducible output
        unique_kmers[gene] = uniq
        print(f"{gene}: {len(uniq)} unique {k}-mers in 5' UTR")
        if uniq:
            # Show first few
            for kmer in uniq[:3]:
                gc = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
                print(f"   {kmer} GC={gc:.1f}%")

    # Same for 3' UTR
    print("\nSearching for unique 20-mers in 3' UTRs...")
    three_kmers = {}
    for gene, utr in three_utrs.items():
        kmers = set()
        for i in range(len(utr) - k + 1):
            kmers.add(utr[i:i+k])
        three_kmers[gene] = kmers

    kmer_counts3 = defaultdict(int)
    for gene, kmers in three_kmers.items():
        for kmer in kmers:
            kmer_counts3[kmer] += 1

    unique_kmers3 = {}
    for gene, kmers in three_kmers.items():
        uniq = sorted(kmer for kmer in kmers if kmer_counts3[kmer] == 1)  # sorted: reproducible output
        unique_kmers3[gene] = uniq
        print(f"{gene}: {len(uniq)} unique {k}-mers in 3' UTR")
        if uniq:
            for kmer in uniq[:3]:
                gc = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
                print(f"   {kmer} GC={gc:.1f}%")

    # Write unique kmers to file for primer design
    with open(P.proc('cornus_unique_5utr_kmers.tsv'), 'w') as f:
        f.write("Gene\tKmer\tGC\tPosition\n")
        for gene, utr in five_utrs.items():
            uniq = unique_kmers.get(gene, [])
            for kmer in uniq:
                # find position(s)
                pos = utr.find(kmer)
                gc = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
                f.write(f"{gene}\t{kmer}\t{gc:.1f}\t{pos}\n")

    with open(P.proc('cornus_unique_3utr_kmers.tsv'), 'w') as f:
        f.write("Gene\tKmer\tGC\tPosition\n")
        for gene, utr in three_utrs.items():
            uniq = unique_kmers3.get(gene, [])
            for kmer in uniq:
                pos = utr.find(kmer)
                gc = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
                f.write(f"{gene}\t{kmer}\t{gc:.1f}\t{pos}\n")

    print("\nUnique kmer tables written.")

if __name__ == '__main__':
    main()
