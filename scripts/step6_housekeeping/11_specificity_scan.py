#!/usr/bin/env python3
"""
Whole‑transcriptome exact‑site scan + in‑silico PCR.
Fixes: proper contig matching via NODE ID, reports family and foreign hits.
"""
import os, sys, re
from collections import defaultdict


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P

def rc(s):
    comp={'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp[c] for c in reversed(s))

def load_fasta(path):
    seqs={}
    name,buf=None,[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if name: seqs[name]=''.join(buf)
                name,buf=line[1:].split()[0],[]
            else: buf.append(line.upper())
        if name: seqs[name]=''.join(buf)
    return seqs

def find_occ(primer, seq):
    occ=[]
    start=0
    while True:
        i=seq.find(primer,start)
        if i<0: break
        occ.append(i); start=i+1
    return occ

def contig2node(contig):
    m=re.search(r'NODE_\d+',contig)
    return m.group() if m else ''

def label2node(label):
    m=re.search(r'NODE_\d+',label)
    return m.group() if m else ''

FAMILIES=['ACT','BTUB','GAPDH','H3','eIF5A','polyUb']
FAMILY_MEMBERS=load_fasta(P.raw('family_transcripts.fasta'))
fam_labels=defaultdict(set)
for label in FAMILY_MEMBERS:
    for fam in FAMILIES:
        if label.startswith(fam):
            fam_labels[fam].add(label)

def main():
    primers={}
    tsv=P.results('FAMILY_primer_pairs.tsv')
    if not os.path.exists(tsv): sys.exit("Run design_family_pairs.py first")
    with open(tsv) as f:
        next(f)
        for line in f:
            p=line.strip().split('\t')
            if len(p)>=3: primers[p[0]]=(p[1].upper(),p[2].upper())

    print("Loading transcriptome...")
    tx=load_fasta(P.TRANSCRIPTOME)
    print(f"  {len(tx)} contigs")

    results_per_family = {}

    for fam in FAMILIES:
        if fam not in primers:
            print(f"\n{fam}: no pair"); results_per_family[fam] = []; continue
        fwd,rev=primers[fam]; rev_rc=rc(rev)
        print(f"\n--- {fam}: F={fwd} R={rev} ---")

        # Build node->label mapping for expected family members
        node2label={}
        for l in fam_labels[fam]:
            n=label2node(l)
            if n: node2label[n]=l

        family_hits={}
        foreign_hits={}
        fwd_sites_global=0; rev_sites_global=0

        for contig,seq in tx.items():
            fwd_sites=find_occ(fwd,seq)
            rev_sites=find_occ(rev_rc,seq)
            if fwd_sites: fwd_sites_global+=1
            if rev_sites: rev_sites_global+=1
            if not fwd_sites or not rev_sites: continue

            cn=contig2node(contig)
            in_family=cn in node2label

            for fi in fwd_sites:
                for ri in rev_sites:
                    if ri>fi and ri-fi<400:
                        prod=ri+20-fi
                        if in_family:
                            family_hits.setdefault(contig,set()).add(prod)
                        else:
                            foreign_hits.setdefault(contig,set()).add(prod)

        # Collect lines for this family
        lines = []
        def record(msg):
            lines.append(msg); print(msg)

        record(f"  Fwd exact sites: {fwd_sites_global} contigs, Rev: {rev_sites_global}")
        record(f"  Family amplicons: {len(family_hits)}")
        for c,prods in sorted(family_hits.items())[:10]:
            record(f"    {c}: {', '.join(f'{p}bp' for p in sorted(prods))}")
        if foreign_hits:
            record(f"  ⚠ FOREIGN amplicons: {len(foreign_hits)}")
            for c,prods in sorted(foreign_hits.items())[:5]:
                record(f"    {c}: {', '.join(f'{p}bp' for p in sorted(prods))}")
        else:
            record(f"  ✅ No foreign amplicons (exact match)")

        # Report unamplified family members
        nodes_amplified = set()
        for c in family_hits:
            n=contig2node(c)
            if n: nodes_amplified.add(n)
        not_amp = [(n,l) for n,l in node2label.items() if n not in nodes_amplified]
        if not_amp:
            record(f"  ⚠ Unamplified family members ({len(not_amp)}):")
            for n,l in not_amp:
                record(f"    {l} (node {n})")
        else:
            record(f"  ✅ All family members amplified")

        # Store for the report
        results_per_family[fam] = lines

    # Write report file
    report_path = P.spec('sensitivuty_report.md')
    with open(report_path, 'w') as f:
        f.write("# Family-level primer pairs — Whole-transcriptome specificity scan\n\n")
        for fam in FAMILIES:
            f.write(f"## {fam}\n\n")
            if fam in primers:
                f.write(f"Forward: `{primers[fam][0]}`  \nReverse: `{primers[fam][1]}`  \n\n")
            f.write('\n'.join(results_per_family.get(fam, ["No pair available."])))
            f.write("\n\n---\n\n")
    print(f"\nReport written: {report_path}")

if __name__=='__main__':
    main()
