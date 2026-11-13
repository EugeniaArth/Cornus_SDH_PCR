#!/usr/bin/env python3
import gzip, re, sys, os

# Path to GFF3 (assume uncompressed)
gff_path = 'PCR_Cornus_SDH/data/external/transcriptome.gff3'

# Families and their keyword patterns
families = {
    'ACT': ['actin', 'Actin'],  # but exclude 'interacting'
    'BTUB': ['tubulin', 'Tubulin'],
    'EF1a': ['elongation factor 1', 'EF-1', 'EF1'],
    'GAPDH': ['glyceraldehyde-3-phosphate dehydrogenase', 'GAPDH'],
    'H3': ['histone H3', 'Histone H3', 'H3.2', 'H3.3', 'CENH3'],
    'TEF1': ['translation elongation factor 1', 'TEF1'],
    'eIF5A': ['translation initiation factor 5A', 'eIF5A', 'IF5A'],
    'UBQ': ['ubiquitin', 'Ubiquitin'],

}

# Compile regex patterns for each family
patterns = {}
for fam, keywords in families.items():
    # Build regex that matches any keyword (case-insensitive)
    pat = re.compile('|'.join(keywords), re.IGNORECASE)
    patterns[fam] = pat

def parse_gff(gff_file):
    """Parse GFF3 and return dict of contig->set(product strings)"""
    contig_products = {}
    line_count = 0
    with open(gff_file, 'r') as f:
        for line in f:
            line_count += 1
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            attrs = parts[8]
            # Parse attributes for product
            product_match = re.search(r'product=([^;]+)', attrs)
            if product_match:
                product = product_match.group(1)
                contig = parts[0]
                contig_products.setdefault(contig, set()).add(product)
    print(f"Parsed {line_count} lines, found {len(contig_products)} contigs with product annotation")
    return contig_products

print("Searching for candidate reference genes in Cornus transcriptome annotation...")
contig_products = parse_gff(gff_path)

# For each family, count contigs matching pattern
family_hits = {}
for fam, pat in patterns.items():
    matching_contigs = set()
    for contig, products in contig_products.items():
        for prod in products:
            if pat.search(prod):
                # For ACT, exclude 'interacting'
                if fam == 'ACT' and 'interacting' in prod.lower():
                    continue
                matching_contigs.add(contig)
                break  # count contig once per family
    family_hits[fam] = matching_contigs
    print(f"{fam}: {len(matching_contigs)} contigs")

# Write results
os.makedirs('/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/housekeeping/keyword_search', exist_ok=True)
with open('/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/housekeeping/keyword_search/family_contig_counts.tsv', 'w') as f:
    f.write("Family\tContigCount\n")
    for fam in sorted(family_hits.keys()):
        f.write(f"{fam}\t{len(family_hits[fam])}\n")

# Write contig lists per family
for fam, contigs in family_hits.items():
    if contigs:
        with open(f'/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/housekeeping/keyword_search/{fam}_contigs.txt', 'w') as f:
            for c in sorted(contigs):
                f.write(f"{c}\n")

print("\nDetailed breakdown for ambiguous families:")
# For ACT, show product strings to verify
print("\nACT products (excluding 'interacting'):")
act_products = set()
for contig, products in contig_products.items():
    for prod in products:
        if 'actin' in prod.lower() and 'interacting' not in prod.lower():
            act_products.add(prod)
for p in sorted(act_products)[:20]:
    print(f"  {p}")

print("\nTEF1 vs EF1a distinction:")
tef1_products = set()
ef1a_products = set()
for contig, products in contig_products.items():
    for prod in products:
        if 'translation elongation factor 1' in prod.lower():
            tef1_products.add(prod)
        elif 'elongation factor 1 alpha' in prod.lower() or 'EF-1α' in prod.lower():
            ef1a_products.add(prod)
print(f"TEF1 products: {len(tef1_products)}")
for p in sorted(tef1_products)[:5]:
    print(f"  {p}")
print(f"EF1a products: {len(ef1a_products)}")
for p in sorted(ef1a_products)[:5]:
    print(f"  {p}")

print(f"\nResults written to keyword_search/")
