#!/usr/bin/env python3
import sys
import os
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths as P
import re

def parse_genbank(gb_file):
    """Return dict with accession, length, and features."""
    with open(gb_file, 'r') as f:
        content = f.read()
    
    # Extract accession
    acc_match = re.search(r'ACCESSION\s+(\S+)', content)
    accession = acc_match.group(1) if acc_match else os.path.basename(gb_file).replace('.gb', '')
    
    # Extract locus length
    locus_match = re.search(r'LOCUS\s+\S+\s+(\d+)\s+bp', content)
    length = int(locus_match.group(1)) if locus_match else None
    
    # Find FEATURES section
    features_start = content.find('FEATURES')
    if features_start == -1:
        return {'accession': accession, 'length': length, 'features': []}
    
    # Get lines after FEATURES until ORIGIN or //
    lines = content[features_start:].split('\n')
    features = []
    current_feature = None
    for line in lines:
        line = line.rstrip()
        if line.startswith('ORIGIN') or line.startswith('//'):
            break
        if line.startswith('     ') and not line.startswith('     ' * 2):
            # Feature line like "5'UTR           1..130"
            if current_feature:
                features.append(current_feature)
            # Parse feature line
            # Example: "5'UTR           1..130"
            parts = line.strip().split()
            if len(parts) >= 2:
                feat_type = parts[0]
                location = parts[1]
                current_feature = {'type': feat_type, 'location': location, 'quals': {}}
            else:
                current_feature = None
        elif current_feature and line.strip().startswith('/'):
            # Qualifier line
            qual_line = line.strip()
            if '=' in qual_line:
                key, val = qual_line[1:].split('=', 1)
                current_feature['quals'][key] = val.strip('"')
            else:
                # continuation line? ignore
                pass
    if current_feature:
        features.append(current_feature)
    
    return {'accession': accession, 'length': length, 'features': features}

def main():
    gb_files = P.genbank_records('PZ5186')
    print(f"Found {len(gb_files)} GenBank files")
    
    all_data = []
    for gb in sorted(gb_files):
        data = parse_genbank(gb)
        print(f"\n{gb}: length {data['length']}")
        for feat in data['features']:
            print(f"  {feat['type']} {feat['location']}")
            if 'note' in feat['quals']:
                print(f"    note: {feat['quals']['note']}")
            if 'product' in feat['quals']:
                print(f"    product: {feat['quals']['product']}")
        all_data.append(data)
    
    # Write summary table
    with open(P.proc('genbank_features.tsv'), 'w') as f:
        f.write("Accession\tLength\tGene\t5UTR_start\t5UTR_end\tCDS_start\tCDS_end\t3UTR_start\t3UTR_end\n")
        for data in all_data:
            acc = data['accession']
            length = data['length']
            # Find features
            five_utr = None
            cds = None
            three_utr = None
            gene_name = None
            for feat in data['features']:
                if feat['type'] == "5'UTR":
                    five_utr = feat['location']
                elif feat['type'] == 'CDS':
                    cds = feat['location']
                    if 'note' in feat['quals']:
                        gene_name = feat['quals']['note']
                    elif 'product' in feat['quals']:
                        # extract gene name maybe
                        pass
                elif feat['type'] == "3'UTR":
                    three_utr = feat['location']
            
            # Parse location strings like "1..130" or "131..1627"
            def parse_loc(loc):
                if not loc:
                    return None, None
                m = re.match(r'(\d+)\.\.(\d+)', loc)
                if m:
                    return int(m.group(1)), int(m.group(2))
                return None, None
            
            f5s, f5e = parse_loc(five_utr) if five_utr else (None, None)
            cds_s, cds_e = parse_loc(cds) if cds else (None, None)
            t3s, t3e = parse_loc(three_utr) if three_utr else (None, None)
            
            f.write(f"{acc}\t{length}\t{gene_name}\t{f5s}\t{f5e}\t{cds_s}\t{cds_e}\t{t3s}\t{t3e}\n")
    
    print("\nSaved genbank_features.tsv")

if __name__ == '__main__':
    main()