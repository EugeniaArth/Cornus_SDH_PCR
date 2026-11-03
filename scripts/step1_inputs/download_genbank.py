#!/usr/bin/env python3
"""Download the six C. sericea DQD/SDH GenBank records (PZ518666–PZ518671) from NCBI.

If a record already exists in data/raw/, it is skipped. Requires network access.
"""
import os
import sys
import subprocess
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths as P

ACCESSIONS = ['PZ518666', 'PZ518667', 'PZ518668', 'PZ518669', 'PZ518670', 'PZ518671']
URL_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
PARAMS = 'db=nuccore&rettype=gb&retmode=text'

def download_one(acc):
    out = os.path.join(P.RAW, f'{acc}.gb')
    if os.path.exists(out):
        print(f'{acc}: already present')
        return False
    cmd = ['curl', '-sS', '--fail', '--max-time', '60',
           f'{URL_BASE}?{PARAMS}&id={acc}', '-o', out]
    print(f'{acc}: downloading...', end='', flush=True)
    try:
        subprocess.run(cmd, check=True)
        # quick sanity check: file non‑empty and contains “LOCUS”
        with open(out, 'r') as f:
            head = f.read(200)
            if 'LOCUS' not in head:
                print('WARNING: downloaded file does not look like a GenBank record')
        print(' done')
        return True
    except subprocess.CalledProcessError as e:
        print(f' FAILED ({e})')
        if os.path.exists(out):
            os.remove(out)
        return False

def main():
    os.makedirs(P.RAW, exist_ok=True)
    any_new = False
    for acc in ACCESSIONS:
        if download_one(acc):
            any_new = True
            time.sleep(1)   # polite delay between requests
    if any_new:
        print('\nNew records saved to', P.RAW)
    else:
        print('\nAll six GenBank records are already present.')

if __name__ == '__main__':
    main()