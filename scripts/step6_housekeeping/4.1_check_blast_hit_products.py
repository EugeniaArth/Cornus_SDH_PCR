#!/usr/bin/env python3
import re, sys

gff_path = '/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/external/transcriptome.gff3'

# Read blast hits summary
blast_hits = {}
with open('/Users/eugenianikonorova/Desktop/Work/Plant_sequencing//PCR_Cornus_SDH/data/housekeeping/agrostemma_blast/blast_hits_summary.tsv', 'r') as f:
    next(f)  # skip header
    for line in f:
        if line.strip():
            fam, contig, pident, bitscore, evalue = line.strip().split('\t')
            blast_hits.setdefault(fam, []).append(contig)

# Parse GFF for product info
contig_to_product = {}
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        contig = parts[0]
        attrs = parts[8]
        product_match = re.search(r'product=([^;]+)', attrs)
        if product_match:
            product = product_match.group(1)
            if contig not in contig_to_product:
                contig_to_product[contig] = product
            elif product != contig_to_product[contig]:
                # Keep the first one
                pass

print("=== BLAST hit contigs and their annotated products ===")
for fam in sorted(blast_hits.keys()):
    print(f"\n{fam}:")
    for contig in blast_hits[fam]:
        product = contig_to_product.get(contig, 'NO_PRODUCT')
        print(f"  {contig}: {product}")
