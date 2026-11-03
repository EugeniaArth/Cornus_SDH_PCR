#!/usr/bin/env python3
import csv
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths as P
from collections import defaultdict

def main():
    csv_file = P.align('Cser_transcript_nucleotide_positions.csv')
    # We'll group by gene and AtSDH_position
    rows = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    
    # Build dict gene -> {position -> residue}
    gene_data = defaultdict(dict)
    for row in rows:
        gene = row['Sequence']
        pos = int(row['AtSDH_position'])
        residue = row['residue']
        gene_data[gene][pos] = residue
    
    # List of all positions
    all_positions = sorted(set(pos for data in gene_data.values() for pos in data.keys()))
    genes = sorted(gene_data.keys())
    
    print(f"Genes: {genes}")
    print(f"Number of positions: {len(all_positions)}")
    
    # Find positions where residues differ across genes
    variable_positions = []
    for pos in all_positions:
        residues = [gene_data[g].get(pos, '-') for g in genes]
        if len(set(residues)) > 1:
            variable_positions.append((pos, residues))
    
    print(f"\nFound {len(variable_positions)} variable positions")
    # Show top 20 most variable
    for pos, residues in variable_positions[:20]:
        print(f"AtSDH position {pos}: {dict(zip(genes, residues))}")
    
    # Now extract nucleotide windows around variable positions
    # We need mapping from AtSDH position to transcript_codon_start
    pos_to_codon = {}
    for row in rows:
        pos = int(row['AtSDH_position'])
        gene = row['Sequence']
        if gene not in pos_to_codon:
            pos_to_codon[gene] = {}
        try:
            codon_start = int(row['transcript_codon_start'])
            pos_to_codon[gene][pos] = codon_start
        except (ValueError, KeyError):
            pass
    
    # For each variable position, get codon start for each gene
    print("\nCodon start positions for variable residues:")
    for pos, residues in variable_positions[:10]:
        print(f"\nAtSDH {pos}:")
        for g in genes:
            codon_start = pos_to_codon[g].get(pos)
            if codon_start:
                print(f"  {g}: codon start {codon_start}, AA {gene_data[g].get(pos, '-')}")
    
    # Load cDNA sequences
    seqs = {}
    with open(P.raw('DQD_SDH_possible_nt.fasta'), 'r') as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    seqs[seq_id] = ''.join(seq)
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if seq_id is not None:
            seqs[seq_id] = ''.join(seq)
    
    # Extract windows of 60 nt around codon start
    window_size = 60
    print(f"\nExtracting {window_size} nt windows around variable codons:")
    for pos, residues in variable_positions[:15]:
        print(f"\nAtSDH {pos}:")
        for g in genes:
            codon_start = pos_to_codon[g].get(pos)
            if codon_start and g in seqs:
                seq = seqs[g]
                start = max(0, codon_start - window_size//2)
                end = min(len(seq), codon_start + window_size//2 + 3)  # include codon
                window = seq[start:end]
                print(f"  {g}: {window}")
                # Check uniqueness: compare with windows from other genes
                # We'll compute later
    
    # Now find unique windows: for each gene, find variable positions where the surrounding window is unique
    print("\nSearching for unique windows (>=10 bp difference) across genes...")
    unique_windows = []
    for pos, residues in variable_positions:
        windows = {}
        for g in genes:
            codon_start = pos_to_codon[g].get(pos)
            if codon_start and g in seqs:
                seq = seqs[g]
                start = max(0, codon_start - 30)  # 30 bp upstream
                end = min(len(seq), codon_start + 33)   # include codon + 30 downstream
                windows[g] = seq[start:end]
        # Check if any window is unique (no overlap with others)
        # Simple approach: compare sequences
        for g1, w1 in windows.items():
            is_unique = True
            for g2, w2 in windows.items():
                if g1 == g2:
                    continue
                # Check if w1 is substring of w2 or vice versa
                if w1 in w2 or w2 in w1:
                    is_unique = False
                    break
                # Compute similarity
                # For simplicity, just check exact equality of first 20 bp
                if w1[:20] == w2[:20]:
                    is_unique = False
                    break
            if is_unique and len(w1) >= 20:
                print(f"  AtSDH {pos}: {g1} has unique window length {len(w1)}")
                unique_windows.append((pos, g1, w1))
    
    print(f"\nFound {len(unique_windows)} unique windows")
    # Write to file
    with open(P.proc('variable_unique_windows.tsv'), 'w') as f:
        f.write("AtSDH_position\tGene\tWindowSeq\tWindowLength\n")
        for pos, g, w in unique_windows:
            f.write(f"{pos}\t{g}\t{w}\t{len(w)}\n")
    print("Saved to variable_unique_windows.tsv")

if __name__ == '__main__':
    main()