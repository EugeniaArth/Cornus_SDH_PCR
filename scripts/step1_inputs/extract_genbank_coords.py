#!/usr/bin/env python3
import sys, os, re, csv
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths as P

def parse_genbank(gb_file):
    with open(gb_file, 'r') as f:
        content = f.read()

    acc_match = re.search(r'ACCESSION\s+(\S+)', content)
    accession = acc_match.group(1) if acc_match else os.path.basename(gb_file).replace('.gb', '')

    locus_match = re.search(r'LOCUS\s+\S+\s+(\d+)\s+bp', content)
    length = int(locus_match.group(1)) if locus_match else None

    def_match = re.search(r'DEFINITION\s+(.+?)\n', content, re.DOTALL)
    definition = def_match.group(1).replace('\n', ' ').strip() if def_match else ''

    features_start = content.find('FEATURES')
    if features_start == -1:
        return {'accession': accession, 'length': length, 'definition': definition, 'features': []}

    lines = content[features_start:].split('\n')
    features = []
    current = None
    for line in lines:
        line = line.rstrip()
        if line.startswith('ORIGIN') or line.startswith('//'):
            break
        if line.startswith('     ') and not line.startswith('     ' * 2):
            if current:
                features.append(current)
            parts = line.strip().split()
            if len(parts) >= 2:
                ftype = parts[0]
                loc = parts[1]
                current = {'type': ftype, 'location': loc, 'quals': {}}
            else:
                current = None
        elif current and line.strip().startswith('/'):
            qual = line.strip()
            if '=' in qual:
                key, val = qual[1:].split('=', 1)
                current['quals'][key] = val.strip('"')
            else:
                pass
    if current:
        features.append(current)

    return {'accession': accession, 'length': length, 'definition': definition, 'features': features}

def parse_loc(loc):
    if not loc:
        return None, None
    m = re.match(r'(\d+)\.\.(\d+)', loc)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None

def main():
    gb_paths = P.genbank_records('PZ5186')          # full paths, so CWD does not matter
    if len(gb_paths) < 6:
        sys.exit('FATAL: only %d of the 6 GenBank records are in data/raw/ (%s).\n'
                 '       Stopping rather than writing a short table over the verified one.\n'
                 '       The verified 6-row table is committed in data/processed/.'
                 % (len(gb_paths), ', '.join(os.path.basename(p) for p in gb_paths)))
    print(f"Processing {len(gb_paths)} GenBank files from {P.RAW}")

    # Parse each GB file
    gb_records = []
    for gb in sorted(gb_paths):
        data = parse_genbank(gb)
        accession = data['accession']

        # Extract gene name from note qualifier
        gene_name = None
        for feat in data['features']:
            if feat['type'] == 'CDS' and 'note' in feat['quals']:
                note = feat['quals']['note']
                if 'CserDQD' in note:
                    gene_name = note.split(';')[0].strip()
                    gene_name = gene_name.replace('CserDQD_SDH', 'CserDQD/SDH')
                break

        if not gene_name:
            m = re.search(r'CserDQD/SDH[\d\.]+', data['definition'])
            if m:
                gene_name = m.group(0)

        # Extract UTR/CDS coordinates
        five_start = five_end = None
        cds_start = cds_end = None
        three_start = three_end = None
        for feat in data['features']:
            if feat['type'] == "5'UTR":
                five_start, five_end = parse_loc(feat['location'])
            elif feat['type'] == 'CDS':
                cds_start, cds_end = parse_loc(feat['location'])
            elif feat['type'] == "3'UTR":
                three_start, three_end = parse_loc(feat['location'])

        gb_records.append({
            'accession': accession,
            'gene': gene_name,
            'length': data['length'],
            '5utr_start': five_start,
            '5utr_end': five_end,
            'cds_start': cds_start,
            'cds_end': cds_end,
            '3utr_start': three_start,
            '3utr_end': three_end
        })

        # Report extracted coordinates
        print(f"\n{gene_name} ({accession}):")
        print(f"  GB: CDS {cds_start}-{cds_end}, length {data['length']}")
        if five_start:
            print(f"  5'UTR: {five_start}-{five_end}")
        if three_start:
            print(f"  3'UTR: {three_start}-{three_end}")

    # Write combined table (only GenBank coordinates, no CSV columns)
    with open(P.proc('genbank_coordinates.tsv'), 'w') as f:
        f.write("Accession\tGene\tLength\t5UTR_start\t5UTR_end\tCDS_start\tCDS_end\t3UTR_start\t3UTR_end\n")
        for rec in gb_records:
            f.write(f"{rec['accession']}\t{rec['gene']}\t{rec['length']}\t{rec['5utr_start']}\t{rec['5utr_end']}\t{rec['cds_start']}\t{rec['cds_end']}\t{rec['3utr_start']}\t{rec['3utr_end']}\n")

    print("\nSaved genbank_coordinates.tsv")

    # Summary
    print("\n=== GenBank Feature Summary ===")
    for rec in gb_records:
        if rec['gene']:
            five_len = rec['5utr_end'] - rec['5utr_start'] + 1 if rec['5utr_start'] else 0
            cds_len = rec['cds_end'] - rec['cds_start'] + 1 if rec['cds_start'] else 0
            three_len = rec['3utr_end'] - rec['3utr_start'] + 1 if rec['3utr_start'] else 0
            print(f"{rec['gene']}: {rec['length']} bp, 5'UTR {five_len} bp, CDS {cds_len} bp, 3'UTR {three_len} bp")

if __name__ == '__main__':
    main()
