#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import os, sys, re

def extract_sequences():
    ns = {'a': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

    # Load shared strings
    sst = ET.parse('PCR_Cornus_SDH/data/housekeeping/external/sharedStrings.xml').getroot()
    shared = []
    for si in sst.findall('a:si', ns):
        txt = ''.join((t.text or '') for t in si.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t'))
        shared.append(txt)

    # Load sheet 7 (Table S7)
    sh = ET.parse('PCR_Cornus_SDH/data/housekeeping/external/sheet7.xml').getroot()
    rows = {}
    for row in sh.findall('a:sheetData/a:row', ns):
        d = {}
        for cell in row.findall('a:c', ns):
            v = cell.find('a:v', ns)
            if v is None:
                continue
            val = shared[int(v.text)] if cell.get('t') == 's' else v.text
            d[cell.get('r')] = val
        rows[int(row.get('r'))] = d

    # Extract families and sequences
    families = {}
    current_family = None

    for r in sorted(rows):
        a_val = None
        # Find column A value
        for cell_ref in rows[r]:
            if cell_ref.startswith('A'):
                a_val = rows[r][cell_ref]
                break

        if a_val and not a_val.replace('.', '').isdigit() and not 'Isotig' in a_val and not 'Table' in a_val and not 'According' in a_val:
            current_family = a_val.strip()
            families[current_family] = []
            print(f"Found family: {current_family}")
        elif current_family:
            # Find nucleotide sequence in column H
            seq = None
            for cell_ref in rows[r]:
                if cell_ref.startswith('H'):
                    seq = rows[r][cell_ref]
                    break

            if seq and len(seq) >= 50:
                # Clean up sequence: remove whitespace, keep only ACTG letters and N
                seq_clean = re.sub(r'[^ACGTUNacgtun]', '', seq.upper())
                if len(seq_clean) >= 50:
                    # Get isotig ID (column A)
                    isotig = None
                    for cell_ref in rows[r]:
                        if cell_ref.startswith('A'):
                            isotig = rows[r][cell_ref]
                            break
                    if isotig:
                        families[current_family].append((isotig, seq_clean))

    # Write FASTA files
    os.makedirs('PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs', exist_ok=True)

    family_map = {
        'Actin (ACT)': 'ACT',
        'β-Tubulin (βTUB)': 'BTUB',
        'Elongation factor 1 α (EF1α)': 'EF1a',
        'Glyceraldehyde-3-phosphate dehydrogenase (GAPDH)': 'GAPDH',
        'Histone H3 (H3)': 'H3',
        'Translation elongation factor 1 (TEF1)': 'TEF1',
        'Eukaryotic translation initiation factor 5A1': 'eIF5A',
        'Squalene synthase (SQS)': 'SQS',
        'β-amyrin synthase (BAS)': 'BAS'
    }

    summary = []
    for full_name, seqs in families.items():
        if not seqs:
            continue
        short = family_map.get(full_name, full_name.split()[0])
        outfile = f'PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs/{short}.fasta'
        with open(outfile, 'w') as f:
            for i, (isotig, seq) in enumerate(seqs, 1):
                f.write(f'>{short}_{i:02d}_isotig{isotig}\n')
                # Write sequence in lines of 60 characters
                for j in range(0, len(seq), 60):
                    f.write(seq[j:j+60] + '\n')

        # Also write a concatenated file with all families
        concat_file = f'PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs/ALL_AGROSTEMMA.fasta'
        with open(concat_file, 'a') as cf:
            for i, (isotig, seq) in enumerate(seqs, 1):
                cf.write(f'>{short}_{i:02d}_isotig{isotig}\n{seq}\n')

        summary.append(f"{short}: {len(seqs)} sequences")
        print(f"  {short}: {len(seqs)} sequences written to {outfile}")

    # Write summary
    with open('PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs/SUMMARY.txt', 'w') as sf:
        sf.write("Agrostemma githago sequences extracted from Table S7\n")
        sf.write("="*60 + "\n")
        for line in summary:
            sf.write(line + "\n")

    total_seqs = sum(len(seqs) for seqs in families.values())
    print(f"\nTotal {total_seqs} sequences extracted")

    # Also create a mapping file for BLAST queries
    with open('/Users/eugenianikonorova/Desktop/Work/Plant_sequencing/PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs/family_info.tsv', 'w') as tf:
        tf.write("Family\tShortName\tIsotig\tSequenceLength\n")
        for full_name, seqs in families.items():
            if not seqs:
                continue
            short = family_map.get(full_name, full_name.split()[0])
            for isotig, seq in seqs:
                tf.write(f"{full_name}\t{short}\t{isotig}\t{len(seq)}\n")

    return families

if __name__ == '__main__':
    extract_sequences()
