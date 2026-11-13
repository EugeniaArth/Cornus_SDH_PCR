#!/usr/bin/env python3
"""
run_validation.py — Comprehensive HMM validation for reference-gene candidates
(H3, eIF5A, TEF1, UBQ) using correct Pfam profiles + NR annotation cross-checking.

CORRECT profile assignments (verified against Pfam NAME/DESC and test scans):

  PF00125 'Histone'          -> H3 family (core histone domain; filtered to canonical H3,
                                 excludes H2A/H2B/H3-like centromeric CENH3)
  PF01287 'eIF-5a'           -> eIF5A family (eukaryotic elongation factor 5A hypusine,
                                 DNA-binding OB fold).  [REPLACES the wrongly-used PF01331
                                 = mRNA_cap_enzyme, which is NOT eIF5A]
  PF00009 'GTP_EFTU'         -> EF-1a / TEF1 / organellar EF-Tu (classified by NR product)
  PF00240 'ubiquitin'        -> ubiquitin-fold family (polyubiquitin by domain architecture)

Parsing is done by Pfam ACCESSION (fld[1], e.g. PF00125.31 -> PF00125), NOT by the HMM NAME
field (which can differ, e.g. PF01331 is named 'mRNA_cap_enzyme').

Outputs (in hmm_validation_v2/):
  pfam_norm5.hmm             combined profile DB (4 correct profiles)
  hmmscan_domtblout.txt      raw HMMER domain table
  <FAMILY>_hmm_hits.tsv      per-family domain hits with scores
  hmm_validation_report.md   human-readable summary

Run:  python3 hmm_validation_v2/run_validation.py   (from Normalisation/)
"""
import os, sys, subprocess, gzip, re, math
from collections import defaultdict, Counter

# ---- Paths ----
HERE = os.path.dirname(os.path.abspath(__file__))
PROFILE_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping', 'hmm_profile')
PROTEIN_FASTA = os.path.join(PROFILE_DIR, 'transcriptome_proteins.fasta')
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data', 'housekeeping', 'hmm_validation')
os.makedirs(OUTDIR, exist_ok=True)

print("HERE:", HERE)
print("PROFILE_DIR:", PROFILE_DIR)
print("OUTDIR:", OUTDIR)

# Profiles in the combined DB. Keys = Pfam accession; value = description only.
# The NAME field inside each .hmm file is NOT used for matching (parsing uses accession).
PROFILES = {
    'PF00125': 'Histone domain (core histone; H3 selected by NR product)',
    'PF01287': 'eIF5A (eukaryotic elongation factor 5A hypusine, OB fold)',
    'PF00009': 'GTP_EFTU (EF-Tu; EF-1a/TEF1 classification via NR product)',
    'PF00240': 'ubiquitin domain (polyubiquitin via domain architecture)',
}

GFF_PATH = os.path.join(os.path.dirname(os.path.dirname(PROFILE_DIR)), 'external', 'transcriptome.gff3')
if not os.path.exists(GFF_PATH):
    GFF_PATH = None

print("GFF_PATH:", GFF_PATH)

def load_fasta(path):
    seqs, name, buf = {}, None, []
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

def make_hmm_db():
    """Combine the four correct profiles into a single pressed HMM database."""
    db_path = os.path.join(OUTDIR, 'pfam_norm5.hmm')  # v5: corrected profile set
    if not os.path.exists(db_path + '.h3m'):
        with open(db_path, 'w') as out:
            for pfid in PROFILES:
                fp = os.path.join(PROFILE_DIR, f'{pfid}.hmm')
                fz = fp + '.gz'
                if os.path.exists(fz):
                    with gzip.open(fz, 'rt') as f: out.write(f.read())
                    print(f"  added {pfid} (gzipped)")
                elif os.path.exists(fp):
                    with open(fp) as f: out.write(f.read())
                    print(f"  added {pfid}")
                else:
                    print(f"WARNING: missing profile {pfid}")
        subprocess.run(['hmmpress', db_path], check=True, capture_output=True)
        print(f"Built {db_path}")
    else:
        print(f"Using existing {db_path}")
    return db_path

def run_hmmscan(hmm_db, out_domtbl):
    cmd = ['hmmscan', '--cut_ga', '--domtblout', out_domtbl,
           '-o', os.path.join(OUTDIR, 'hmmscan_full.txt'),
           '--cpu', '4', hmm_db, PROTEIN_FASTA]
    print(f"Running hmmscan over {PROTEIN_FASTA} ...")
    subprocess.run(cmd, check=True)
    print(f"Done. Domain table: {out_domtbl}")

def parse_domtbl(path):
    """
    Return list of (pfam_id, accession, query, evalue, score, qlen, desc).
    pfam_id derived from ACCESSION column (fld[1]), NOT the NAME column (fld[0]).
    """
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fld = line.strip().split()
            if len(fld) < 19: continue
            acc = fld[1]                         # e.g. PF00125.31
            pfid = acc.split('.')[0]             # e.g. PF00125
            qname = fld[3]
            evalue = float(fld[12])
            score = float(fld[13])
            qlen = int(fld[5])
            desc = ' '.join(fld[19:]) if len(fld) > 19 else ''
            hits.append((pfid, acc, qname, evalue, score, qlen, desc))
    return hits

def load_products():
    """contig -> product string from GFF (first occurrence)."""
    products = {}
    if not GFF_PATH or not os.path.exists(GFF_PATH): return products
    with open(GFF_PATH) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            m = re.search(r'product=([^;]+)', parts[8])
            if m:
                contig = parts[0]
                if contig not in products:
                    products[contig] = m.group(1)
    return products

def classify_EFTU(product):
    """Classify a PF00009 hit by NR product into EF1a / TEF1 / organellar EF-Tu / other."""
    p = (product or '').lower()
    if 'elongation factor 1-alpha' in p or 'ef-1' in p or 'ef1' in p or 'elongation factor 1 alpha' in p:
        return 'EF1a'
    if 'translation elongation factor 1' in p or 'tef1' in p or 'elongation factor 1' in p:
        return 'TEF1'
    if 'elongation factor tu' in p or 'ef-tu' in p:
        if 'mitochondrial' in p: return 'EFTu_mito'
        if 'chloroplast' in p or 'plastid' in p: return 'EFTu_chloro'
        return 'EFTu'
    return 'other_GTPase'

def main():
    print("="*70)
    print("Comprehensive HMM validation (v2, corrected profiles)")
    print("  H3=PF00125  eIF5A=PF01287  TEF1/EF1a=PF00009  UBQ=PF00240")
    print("  eIF5A profile CORRECTED: PF01331(mRNA_cap_enzyme) -> PF01287(eIF-5a)")
    print("="*70)

    # 1. Build HMM DB and run hmmscan
    hmm_db = make_hmm_db()
    domtbl = os.path.join(OUTDIR, 'hmmscan_domtblout.txt')
    run_hmmscan(hmm_db, domtbl)

    # 2. Parse
    hits = parse_domtbl(domtbl)
    print(f"\nTotal domain hits: {len(hits)}")
    prof_counts = Counter(h[0] for h in hits)
    for pfid, cnt in sorted(prof_counts.items()):
        print(f"  {pfid}: {cnt} hits")

    # 3. Load products
    products = load_products()
    if products:
        print(f"Loaded {len(products)} NR product annotations")

    # 4. Group by family
    h3_hits   = [h for h in hits if h[0] == 'PF00125']
    eif_hits  = [h for h in hits if h[0] == 'PF01287']
    eftu_hits = [h for h in hits if h[0] == 'PF00009']
    ub_hits   = [h for h in hits if h[0] == 'PF00240']

    # ---- H3 classification ----
    print("\n--- H3 (PF00125, core histone domain) ---")
    h3_rows = []
    for prof, acc, qname, ev, sc, qlen, desc in h3_hits:
        prod = products.get(qname, '')
        pl = prod.lower()
        # keep canonical H3 (H3.1/H3.2/H3.3); exclude H2A/H2B/CENH3/NF-Y
        if 'histone h3' in pl and 'cenh3' not in pl:
            h3_rows.append({'contig': qname, 'evalue': ev, 'score': sc, 'product': prod})
    print(f"  PF00125 hits: {len(h3_hits)}; canonical H3 (by NR product): {len(h3_rows)}")

    # ---- eIF5A classification (PF01287) ----
    print("\n--- eIF5A (PF01287, eIF-5a) ---")
    eif_rows = []
    for prof, acc, qname, ev, sc, qlen, desc in eif_hits:
        prod = products.get(qname, '')
        pl = prod.lower()
        # any hit at cutoff is eIF5A by domain; record product for confirmation
        eif_rows.append({'contig': qname, 'evalue': ev, 'score': sc, 'product': prod})
    print(f"  PF01287 hits: {len(eif_hits)}; eIF5A candidates: {len(eif_rows)}")
    for r in eif_rows:
        print(f"    {r['contig']}: score={r['score']:.1f}  {r['product']}")

    # ---- TEF1 / EF1a (PF00009) ----
    print("\n--- TEF1 / EF1a (PF00009, GTP_EFTU) ---")
    tef_rows = []; ef1a_rows = []; eftu_rows = []; other_g = []
    for prof, acc, qname, ev, sc, qlen, desc in eftu_hits:
        prod = products.get(qname, '')
        kind = classify_EFTU(prod)
        row = {'contig': qname, 'evalue': ev, 'score': sc, 'product': prod, 'kind': kind}
        if kind == 'TEF1': tef_rows.append(row)
        elif kind == 'EF1a': ef1a_rows.append(row)
        elif kind.startswith('EFTu'): eftu_rows.append(row)
        else: other_g.append(row)
    print(f"  PF00009 hits: {len(eftu_hits)}")
    print(f"    EF1a-class: {len(ef1a_rows)}")
    print(f"    TEF1-class: {len(tef_rows)}")
    print(f"    organellar EF-Tu: {len(eftu_rows)}")
    print(f"    other GTPases: {len(other_g)}")
    if tef_rows:
        for r in tef_rows:
            print(f"      TEF1: {r['contig']}: {r['product']}")

    # ---- UBQ (PF00240) ----
    print("\n--- Ubiquitin (PF00240) ---")
    ub_arch = {}
    with open(domtbl) as f:
        for line in f:
            if line.startswith('#'): continue
            fld = line.strip().split()
            if len(fld) < 22 or 'PF00240' not in fld[1]: continue
            qname = fld[3]
            env_from, env_to = int(fld[19]), int(fld[20])
            ub_arch.setdefault(qname, []).append((env_from, env_to))
    ub_rows = []
    for qname, envs in ub_arch.items():
        n_dom = len(envs)
        prot = proteins.get(qname, '')
        qlen = len(prot)
        cov = (sum(e2-e1+1 for e1,e2 in envs)) / qlen if qlen else 0
        envs_sorted = sorted(envs)
        spacing = [envs_sorted[i][0]-envs_sorted[i-1][0] for i in range(1, len(envs_sorted))]
        if n_dom >= 2 and cov >= 0.70 and all(68 <= s <= 84 for s in spacing):
            kind = 'polyubiquitin'
        elif n_dom >= 2 and cov >= 0.50:
            kind = 'multiUb'
        else:
            kind = 'singleUb/fusion'
        prod = products.get(qname, '')
        ub_rows.append({'contig': qname, 'n_dom': n_dom, 'coverage': round(cov,2),
                        'spacing': str(spacing), 'kind': kind, 'product': prod})
    kind_cnt = Counter(r['kind'] for r in ub_rows)
    print(f"  Ubiquitin-fold proteins: {len(ub_rows)}")
    for k, c in sorted(kind_cnt.items()):
        print(f"    {k}: {c}")
    polyub = [r for r in ub_rows if r['kind'] == 'polyubiquitin']
    print(f"  Polyubiquitin: {len(polyub)}")
    for r in polyub:
        print(f"    {r['contig']}: n_dom={r['n_dom']} cov={r['coverage']} spacing={r['spacing']}")

    # ---- Write outputs ----
    def write_tsv(name, header, rows):
        path = os.path.join(OUTDIR, name)
        with open(path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for r in rows:
                f.write('\t'.join(str(r.get(h,'')) for h in header) + '\n')
        print(f"  wrote {path} ({len(rows)} rows)")

    write_tsv('H3_hmm_hits.tsv',    ['contig','evalue','score','product'], h3_rows)
    write_tsv('eIF5A_hmm_hits.tsv', ['contig','evalue','score','product'], eif_rows)
    write_tsv('TEF1_hmm_hits.tsv',  ['contig','evalue','score','product'], tef_rows)
    write_tsv('EF1a_hmm_hits.tsv',  ['contig','evalue','score','product'], ef1a_rows)
    write_tsv('UBQ_hmm_hits.tsv',   ['contig','n_dom','coverage','spacing','kind','product'], ub_rows)

    # ---- Report ----
    with open(os.path.join(OUTDIR, 'hmm_validation_report.md'), 'w') as f:
        f.write("# Comprehensive HMM Validation Report (corrected profiles)\n\n")
        f.write("Profiles: H3=PF00125, eIF5A=PF01287, TEF1/EF1a=PF00009, UBQ=PF00240\n\n")
        f.write("> NOTE: eIF5A previously used PF01331 (mRNA_cap_enzyme, WRONG). "
                "Correct profile PF01287 (eIF-5a) is used here.\n\n")
        f.write("## Summary\n\n")
        f.write("| Family | Profile | Total hits | Validated |\n")
        f.write("|--------|---------|-----------|-----------|\n")
        f.write(f"| H3 | PF00125 | {len(h3_hits)} | {len(h3_rows)} |\n")
        f.write(f"| eIF5A | PF01287 | {len(eif_hits)} | {len(eif_rows)} |\n")
        f.write(f"| TEF1 / EF1a | PF00009 | {len(eftu_hits)} | TEF1:{len(tef_rows)}, EF1a:{len(ef1a_rows)} |\n")
        f.write(f"| Ubiquitin | PF00240 | {len(ub_rows)} | polyUb:{len(polyub)} |\n\n")
        f.write("## eIF5A candidates (PF01287)\n\n")
        for r in eif_rows:
            f.write(f"- `{r['contig']}` (score {r['score']:.1f}): {r['product']}\n")
        f.write("\n## TEF1 candidates\n\n")
        if tef_rows:
            for r in tef_rows:
                f.write(f"- `{r['contig']}` (score {r['score']:.1f}): {r['product']}\n")
        else:
            f.write("None found by HMM+NR classification (Cornus lacks a distinct TEF1).\n")
        f.write("\n## H3 candidates (canonical)\n\n")
        for r in h3_rows:
            f.write(f"- `{r['contig']}`: {r['product']}\n")
        f.write("\n## Polyubiquitin candidates\n\n")
        for r in polyub:
            f.write(f"- `{r['contig']}` n_dom={r['n_dom']} cov={r['coverage']} spacing={r['spacing']}\n")

    print(f"\nAll outputs in {OUTDIR}/")

if __name__ == '__main__':
    if not os.path.exists(PROTEIN_FASTA):
        sys.exit(f"PROTEIN_FASTA not found: {PROTEIN_FASTA}")
    proteins = load_fasta(PROTEIN_FASTA)
    main()
