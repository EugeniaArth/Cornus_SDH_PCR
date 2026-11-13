#!/usr/bin/env python3
"""
step0_inputs/extract_filtered_candidates.py — Build filtered family FASTA files
from combined BLAST, keyword, and HMM validation evidence.

For each family, determines the authoritative candidate contig set:
  ACT, BTUB, GAPDH    — high‑confidence BLAST intersection (BLAST ∩ keyword)
  EF1a                — HMM PF00009 (GTP_EFTU) EF1a-class (9 contigs)
  H3                  — HMM PF00125, filtered to canonical H3.2/H3.3 by NR product
  eIF5A               — HMM PF01287 (eIF-5a, CORRECTED profile)
  polyUb              — HMM PF00240, polyubiquitin by architecture rule
  TEF1                — none (confirmed absent in Cornus)

Outputs:
  family_candidates_filtered/<FAMILY>_filtered.fasta  — extracted transcript sequences
  family_candidates_filtered/filtered_summary.tsv      — contig × product × evidence

Run from family_primer_design/:  python3 scripts/step0_inputs/extract_filtered_candidates.py
"""
import os, re, sys
from collections import OrderedDict

# Resolve paths relative to script location
HERE = os.path.dirname(os.path.abspath(__file__))

# base directory
_PROJECT = os.path.abspath(os.path.join(HERE, '..', '..'))

# PATHs
HOUSEKEEPING = os.path.join(_PROJECT, 'data', 'housekeeping')

TRANSCRIPTOME = os.path.join(_PROJECT, 'data','external', 'transcriptome.fasta')
GFF_PATH = os.path.join(_PROJECT, 'data','external', 'transcriptome.gff3')
HMM_DIR = os.path.join(HOUSEKEEPING, 'hmm_validation')
BLAST_FILE = os.path.join(HOUSEKEEPING,'agrostemma_blast', 'blast_hits_summary.tsv')
KW_DIR = os.path.join(HOUSEKEEPING, 'keyword_search')
OUTDIR = os.path.join(HOUSEKEEPING, 'family_candidates_filtered')
os.makedirs(OUTDIR, exist_ok=True)

FAMILIES = OrderedDict([
    ('ACT',   {'source': 'blast_kw_intersect', 'desc': 'Actin'}),
    ('BTUB',  {'source': 'blast_kw_intersect', 'desc': 'Beta-tubulin'}),
    ('EF1a',  {'source': 'hmm_ef1a',           'desc': 'Elongation factor 1-alpha'}),
    ('GAPDH', {'source': 'blast_kw_intersect', 'desc': 'GAPDH cytosolic'}),
    ('H3',    {'source': 'hmm_h3',             'desc': 'Histone H3 (canonical)'}),
    ('eIF5A', {'source': 'hmm_eif5a',          'desc': 'eIF5A'}),
    ('polyUb',{'source': 'hmm_polyub',         'desc': 'Polyubiquitin'}),
])

def load_fasta(path):
    seqs, name, buf = OrderedDict(), None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper())
        if name: seqs[name] = ''.join(buf)
    return seqs

def extract_seqs(contig_ids, tx):
    """Return OrderedDict {contig_id: sequence} for all requested IDs."""
    result = OrderedDict()
    for c in contig_ids:
        if c in tx:
            result[c] = tx[c]
        else:
            # maybe the contig ID doesn't have the full annotation suffix
            full = [k for k in tx if k.startswith(c + '_') or k.startswith(c)]
            if full:
                result[full[0]] = tx[full[0]]
    return result

def load_products():
    """contig -> product string from GFF."""
    products = {}
    if not os.path.exists(GFF_PATH): return products
    with open(GFF_PATH) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            m = re.search(r'product=([^;]+)', parts[8])
            if m and parts[0] not in products:
                products[parts[0]] = m.group(1)
    return products

def read_hmm_hits(tsv_path):
    """Read HMM hit TSV, return list of dicts."""
    if not os.path.exists(tsv_path): return []
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        rows = []
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= len(header):
                rows.append(dict(zip(header, p)))
    return rows

def read_blast():
    """return dict family -> set of contig IDs (BLAST hits)."""
    b = {}
    if not os.path.exists(BLAST_FILE): return b
    with open(BLAST_FILE) as f:
        next(f)
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                fam = p[0]
                b.setdefault(fam, []).append(p[1])
    return b

def read_keyword(fam_short):
    """return set of contig IDs for a family from keyword search."""
    path = os.path.join(KW_DIR, f'{fam_short}_contigs.txt')
    if not os.path.exists(path): return set()
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

def main():
    print("=== Step 0: Extract filtered family candidates from evidence ===\n")

    # 1. Load transcriptome
    print(f"Loading transcriptome from {TRANSCRIPTOME}...")
    tx = load_fasta(TRANSCRIPTOME)
    print(f"  {len(tx)} contigs loaded")
    products = load_products()
    print(f"  {len(products)} product annotations loaded")

    # 2. Read evidence
    h3_hits = read_hmm_hits(os.path.join(HMM_DIR, 'H3_hmm_hits.tsv'))
    eif_hits = read_hmm_hits(os.path.join(HMM_DIR, 'eIF5A_hmm_hits.tsv'))
    ubq_hits = read_hmm_hits(os.path.join(HMM_DIR, 'UBQ_hmm_hits.tsv'))
    ef1a_hits = read_hmm_hits(os.path.join(HMM_DIR, 'EF1a_hmm_hits.tsv'))
    blast = read_blast()

    print(f"\nEvidence summary:")
    print(f"  H3 (canonical): {len(h3_hits)} contigs")
    print(f"  eIF5A (PF01287): {len(eif_hits)} contigs")
    print(f"  EF1a (PF00009 EF1a-class): {len(ef1a_hits)} contigs")
    print(f"  Ubiquitin total: {len(ubq_hits)}, polyubiquitin: {sum(1 for r in ubq_hits if r.get('kind') == 'polyubiquitin')}")

    # 3. Determine per-family contig sets
    family_contigs = {}

    # ACT: BLAST ∩ keyword intersection
    bl_act = set(blast.get('ACT', []))
    kw_act = read_keyword('ACT')
    inter_act = bl_act & kw_act
    family_contigs['ACT'] = list(inter_act) if inter_act else bl_act
    print(f"\nACT: BLAST={len(bl_act)} kw={len(kw_act)} ∩={len(inter_act)} → {len(family_contigs['ACT'])} contigs")

    # BTUB: BLAST ∩ keyword
    bl_btub = set(blast.get('BTUB', []))
    kw_btub = read_keyword('BTUB')
    inter_btub = bl_btub & kw_btub
    family_contigs['BTUB'] = list(inter_btub) if inter_btub else [c for c in bl_btub if 'tubulin beta' in (products.get(c,'').lower())]
    print(f"BTUB: BLAST={len(bl_btub)} kw={len(kw_btub)} ∩={len(inter_btub)} → {len(family_contigs['BTUB'])} contigs")

    # GAPDH: BLAST ∩ keyword
    bl_gapdh = set(blast.get('GAPDH', []))
    kw_gapdh = read_keyword('GAPDH')
    inter_gapdh = bl_gapdh & kw_gapdh
    family_contigs['GAPDH'] = list(inter_gapdh) if inter_gapdh else bl_gapdh
    print(f"GAPDH: BLAST={len(bl_gapdh)} kw={len(kw_gapdh)} ∩={len(inter_gapdh)} → {len(family_contigs['GAPDH'])} contigs")

    # EF1a: HMM EF1a-class contigs
    fam_ef1a = set(r['contig'] for r in ef1a_hits)
    family_contigs['EF1a'] = list(fam_ef1a)
    print(f"EF1a (HMM EF1a-class): {len(family_contigs['EF1a'])} contigs")

    # H3: canonical H3 from HMM (already filtered)
    fam_h3 = set(r['contig'] for r in h3_hits)
    family_contigs['H3'] = list(fam_h3)
    print(f"H3 (HMM canonical): {len(family_contigs['H3'])} contigs")

    # eIF5A: PF01287 hits
    fam_eif = set(r['contig'] for r in eif_hits)
    family_contigs['eIF5A'] = list(fam_eif)
    print(f"eIF5A (HMM PF01287): {len(family_contigs['eIF5A'])} contigs")

    # polyUb: polyubiquitin by architecture + product filter (remove Fragment/-like/RUB1/NEDD8)
    fam_polyub_raw = [r for r in ubq_hits if r.get('kind') == 'polyubiquitin']
    fam_polyub_filtered = []
    removed = []
    for r in fam_polyub_raw:
        prod = (products.get(r['contig'], '') or '').lower()
        if 'fragment' in prod or '-like protein' in prod or 'nedd8' in prod or 'rub1' in prod:
            removed.append((r['contig'], products.get(r['contig'], '')))
        else:
            fam_polyub_filtered.append(r['contig'])
    family_contigs['polyUb'] = fam_polyub_filtered
    print(f"polyUb (HMM architecture): {len(fam_polyub_raw)} total → {len(fam_polyub_filtered)} after product filter")
    if removed:
        print(f"  Removed by product filter ({len(removed)}):")
        for c, p in removed:
            print(f"    {c}: {p}")

    # 4. Write output
    os.makedirs(OUTDIR, exist_ok=True)

    summary_rows = []
    for fam, contig_ids in family_contigs.items():
        if not contig_ids:
            print(f"\n  {fam}: NO CONTIGS, skipping")
            continue

        # Extract sequences
        seqs = extract_seqs(contig_ids, tx)
        missing = [c for c in contig_ids if c not in seqs and c[:] not in seqs]

        out_fasta = os.path.join(OUTDIR, f'{fam}_filtered.fasta')
        with open(out_fasta, 'w') as f:
            for c, s in seqs.items():
                f.write(f'>{c}\n')
                for i in range(0, len(s), 60):
                    f.write(s[i:i+60] + '\n')
        print(f"  {fam}: {len(seqs)} sequences → {out_fasta}")
        if missing:
            print(f"    WARNING: {len(missing)} contigs not found in transcriptome: {missing[:3]}")

        # Add to summary TSV
        for c in seqs:
            prod = products.get(c, '')
            summary_rows.append({'family': fam, 'contig': c, 'product': prod, 'source': FAMILIES[fam]['source']})

    # Write summary TSV
    tsv_path = os.path.join(OUTDIR, 'filtered_summary.tsv')
    with open(tsv_path, 'w') as f:
        f.write('Family\tContig\tProduct\tSource\n')
        for r in summary_rows:
            f.write(f"{r['family']}\t{r['contig']}\t{r['product']}\t{r['source']}\n")
    print(f"\nSummary written: {tsv_path}")
    print(f"Total: {len(summary_rows)} contigs in {len(family_contigs)} families")

if __name__ == '__main__':
    main()
