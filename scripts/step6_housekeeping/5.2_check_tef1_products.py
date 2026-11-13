#!/usr/bin/env python3
import re, sys

gff_path = '/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/external/transcriptome.gff3'

# Load TEF1 contigs
with open('/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/housekeeping/keyword_search/TEF1_contigs.txt', 'r') as f:
    tef1_contigs = [line.strip() for line in f if line.strip()]

# Parse GFF for these contigs
contig_to_product = {}
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        contig = parts[0]
        if contig in tef1_contigs:
            attrs = parts[8]
            product_match = re.search(r'product=([^;]+)', attrs)
            if product_match:
                product = product_match.group(1)
                if contig not in contig_to_product:
                    contig_to_product[contig] = product

print("TEF1 keyword contigs and their products:")
for contig in tef1_contigs:
    product = contig_to_product.get(contig, 'NO_PRODUCT')
    print(f"{contig}: {product}")
