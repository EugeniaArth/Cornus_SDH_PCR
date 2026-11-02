#!/usr/bin/env python3
"""
Whole-transcriptome (genome-wide) off-target check for the final Cornus DQD/SDH
primer set.  Runs fully LOCALLY with blastn; nothing is uploaded.

Template: Transcriptome/final.clust_transcripts_longest_iso.fasta (29,124 Trinity
          transcripts, ~48 MB) - the annotated, clustered transcript assembly
          whose IDs match Transcriptome/final.clust_annotation_longest_iso.gff3
          and ./annotated_final_nr.gff3.

Steps
  0. makeblastdb on the transcriptome (DB written inside primer_blast_check/).
  1. Map the six curated DQD/SDH transcripts (PZ518666-71) onto the assembly so
     that hits falling on a target contig are recognised as the intended gene.
  2. blastn-short each of the 12 primers against all 29,124 transcripts.
  3. Rebuild candidate amplicons with Primer-BLAST-style tolerance:
     full-length primer alignment, 80-600 bp product, <=3 total mismatches,
     perfect 5-nt 3'-terminal match, F/R facing each other on the same contig.
  4. Classify every candidate product:
        TARGET_CONTIG        - contig = one of the six genes (expected band)
        SDH_FAMILY_CONTIG    - annotated as SDH/DQD homolog but not one of the six
        UNRELATED_CONTIG     - true off-target candidate
     and pull the annotation product/name for each hit contig.
  5. Count strong single-primer priming sites transcriptome-wide.

Outputs (all inside primer_blast_check/):
  cornus_tx.n*              BLAST database
  tx_selfmap.tsv            six curated transcripts -> assembly contigs
  tx_blast_hits.tsv         raw blastn-short primer hits
  tx_priming_sites.tsv      strong single-primer sites (full-length + 3' perfect)
  tx_product_check.tsv      candidate amplicons + classification + annotation
  TRANSCRIPTOME_CHECK.md    human-readable verdict report
"""
import csv, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

TX_FASTA   = os.path.join(PARENT, 'Transcriptome', 'final.clust_transcripts_longest_iso.fasta')
GFF_NR     = os.path.join(PARENT, 'annotated_final_nr.gff3')          # rich annotation (nr blast)
GFF_TX     = os.path.join(PARENT, 'Transcriptome', 'final.clust_annotation_longest_iso.gff3')
TARGETS    = os.path.join(PARENT, 'DQD_SDH_possible_nt.fasta')
PRIMERS    = os.path.join(PARENT, 'FINAL_primers.fasta')
DB         = os.path.join(HERE, 'cornus_tx')

MIN_LEN, MAX_LEN = 80, 600
MAX_MM_TOTAL = 3
MM_3PRIME_NEED = 5

FAMILY_RE = re.compile(r'succinate dehydrogenase|dehydrogenase \[ubiquinone\]|'
                       r'SDH|fumarate reductase|flavoprotein subunit', re.I)

def rc(s):
    c = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(c.get(ch, 'N') for ch in reversed(s))

def run(cmd, inp=None):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write(' '.join(cmd) + '\n' + r.stderr + '\n')
        sys.exit(1)
    return r.stdout

def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if name: recs.append((name, ''.join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name: recs.append((name, ''.join(buf)))
    return recs

def blastn_short(queries, extra=()):
    cols = ('qseqid sseqid pident length mismatch gapopen qstart qend sstart send '
            'evalue bitscore qseq sseq')
    cmd = ['blastn', '-task', 'blastn-short', '-word_size', '7', '-evalue', '1000',
           '-perc_identity', '70', '-max_target_seqs', '5000', '-db', DB,
           '-outfmt', f'6 {cols}'] + list(extra)
    inp = '\n'.join(f'>{n}\n{s}' for n, s in queries)
    return run(cmd, inp)

def megablast(queries):
    cols = 'qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore'
    cmd = ['blastn', '-task', 'megablast', '-evalue', '1e-20',
           '-max_target_seqs', '20', '-db', DB, '-outfmt', f'6 {cols}']
    inp = '\n'.join(f'>{n}\n{s}' for n, s in queries)
    return run(cmd, inp)

def load_annotations():
    """contig -> (gene name, product text) from the nr-annotated GFF, fallback to tx GFF."""
    ann = {}
    for gff in (GFF_NR, GFF_TX):
        if not os.path.isfile(gff): continue
        with open(gff) as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 9 or p[2] != 'mRNA': continue
                cid = p[0]
                if cid in ann and ann[cid][1]: continue
                attrs = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
                name = attrs.get('Name', '')
                prod = attrs.get('product', '') or attrs.get('Note', '')
                if name.startswith('ORF type') : name = ''
                if cid not in ann or (prod and not ann[cid][1]):
                    ann[cid] = (name, prod)
    return ann

def main():
    # ---------------------------------------------------------------- 0. DB
    if not os.path.isfile(TX_FASTA):
        sys.exit(f"missing transcriptome: {TX_FASTA}")
    if not os.path.exists(DB + '.nin'):
        print("Building BLAST DB for the transcriptome (~48 MB, 29k transcripts)...")
        out = run(['makeblastdb', '-in', TX_FASTA, '-dbtype', 'nucl', '-out', DB])
        print('  ' + ' | '.join(out.splitlines()[:3]))

    primers = read_fastas(PRIMERS) if False else read_fasta(PRIMERS)
    plen = {n: len(s) for n, s in primers}
    print(f"{len(primers)} primers loaded")

    # ------------------------------------------------- 1. target-contig map
    tg = read_fasta(TARGETS)
    tg_len = {n: len(s) for n, s in tg}
    print("Mapping the six curated DQD/SDH transcripts onto the assembly...")
    selfmap = {}
    for line in megablast(tg).strip().splitlines():
        if not line: continue
        p = line.split('\t')
        q, s, pid, alen = p[0], p[1], float(p[2]), int(p[3])
        cov = alen / tg_len[q]
        if pid >= 95 and alen >= 300:
            selfmap.setdefault(q, []).append((s, pid, alen, cov))
    with open(os.path.join(HERE, 'tx_selfmap.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['curated_transcript','contig','pident','aln_len','qcov'])
        for q, hits in sorted(selfmap.items()):
            for h in sorted(hits, key=lambda x: -x[2])[:5]:
                w.writerow([q] + list(h))
    target_contigs = {h[0] for hits in selfmap.values() for h in hits}
    for q in tg:
        best = sorted(selfmap.get(q[0], []), key=lambda x: -x[2])[:1]
        print(f"  {q[0]:14s} -> {best[0][0] if best else 'NO >=95% CONTIG'}"
              + (f"  ({best[0][1]:.1f}% id, {best[0][2]} bp)" if best else ''))

    # ------------------------------------------- 2. primer hits (blastn-short)
    print("Running blastn-short of 12 primers against 29,124 transcripts...")
    raw = blastn_short(primers)
    aln, hits_rows = [], []
    for line in raw.strip().splitlines():
        if not line: continue
        p = line.split('\t')
        rec = {'primer': p[0], 'tx': p[1], 'pident': float(p[2]), 'alen': int(p[3]),
               'mm': int(p[4]), 'qstart': int(p[6]), 'qend': int(p[7]),
               'sstart': int(p[8]), 'send': int(p[9]), 'ev': float(p[10]),
               'qseq': p[12], 'sseq': p[13]}
        aln.append(rec)
        hits_rows.append(rec)
    with open(os.path.join(HERE, 'tx_blast_hits.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['primer','contig','pident','aln_len','mismatch','qstart','qend',
                    'sstart','send','evalue'])
        for r in hits_rows:
            w.writerow([r['primer'], r['tx'], r['pident'], r['alen'], r['mm'],
                        r['qstart'], r['qend'], r['sstart'], r['send'], r['ev']])
    print(f"  raw hits: {len(aln)}")

    # ------------------------------------- 3. strong single-primer priming sites
    def full_len(r, side):
        return (r['qstart'] == 1 and r['qend'] == plen[r['primer']]
                and r['alen'] == plen[r['primer']])

    def three_prime_perfect(r, need=MM_3PRIME_NEED):
        good = 0
        for qb, sb in reversed(list(zip(r['qseq'], r['sseq']))):
            if qb.upper() == sb.upper(): good += 1
            else: break
        return good >= need

    strong = [r for r in aln if full_len(r, None) and r['mm'] <= 3 and three_prime_perfect(r)]
    with open(os.path.join(HERE, 'tx_priming_sites.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['primer','contig','pident','mismatch','sstart','send','strand',
                    'is_target_contig'])
        for r in sorted(strong, key=lambda x: (x['primer'], x['tx'])):
            w.writerow([r['primer'], r['tx'], r['pident'], r['mm'], r['sstart'], r['send'],
                        'plus' if r['sstart'] < r['send'] else 'minus',
                        r['tx'] in target_contigs])
    print(f"  strong single-primer sites (full-length, <=3 mm, perfect 3' 5 nt): {len(strong)}")

    # ------------------------------------------------------- 4. candidate amplicons
    pairs = {}
    for n, s in primers:
        gene = n.rsplit('_', 1)[0]
        pairs.setdefault(gene, {})[n.rsplit('_', 1)[1]] = {'name': n, 'seq': s}

    ann = load_annotations()
    results = []
    for gene, d in sorted(pairs.items()):
        fset = [r for r in strong if r['primer'] == d['Fwd']['name'] and r['sstart'] < r['send']]
        rset = [r for r in strong if r['primer'] == d['Rev']['name'] and r['sstart'] > r['send']]
        for fh in fset:
            for rh in rset:
                if fh['tx'] != rh['tx']: continue
                size = rh['sstart'] - fh['sstart'] + 1
                if not (MIN_LEN <= size <= MAX_LEN): continue
                mm = fh['mm'] + rh['mm']
                if mm > MAX_MM_TOTAL: continue
        # keep only products; classify
                if fh['tx'] in target_contigs:
                    klass = 'TARGET_CONTIG'
                elif FAMILY_RE.search(' '.join(ann.get(fh['tx'], ('', '')))):
                    klass = 'SDH_FAMILY_CONTIG'
                else:
                    klass = 'UNRELATED_CONTIG'
                results.append({'pair_gene': gene, 'contig': fh['tx'], 'size': size,
                                'mm_fwd': fh['mm'], 'mm_rev': rh['mm'], 'mm_total': mm,
                                'class': klass,
                                'gene_name': ann.get(fh['tx'], ('', ''))[0],
                                'product': ann.get(fh['tx'], ('', ''))[1][:120]})
    with open(os.path.join(HERE, 'tx_product_check.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['pair_gene','contig','size','mm_fwd','mm_rev',
                                         'mm_total','class','gene_name','product'],
                           delimiter='\t')
        w.writeheader()
        for r in sorted(results, key=lambda x: (x['pair_gene'], x['class'], x['size'])):
            w.writerow(r)

    # ------------------------------------------------------------- 5. report
    own_contigs = {gene: {h[0] for h in selfmap.get(
        gene.replace('CserDQD_SDH', 'CserDQD/SDH'), [])} for gene in pairs}

    with open(os.path.join(HERE, 'TRANSCRIPTOME_CHECK.md'), 'w') as f:
        f.write("# Whole-transcriptome off-target check of the final primer set\n\n")
        f.write(f"Template: `Transcriptome/final.clust_transcripts_longest_iso.fasta` "
                f"(29,124 transcripts)  \n")
        f.write(f"Engine: local `blastn -task blastn-short -word_size 7 -perc_identity 70`, "
                f"Primer-BLAST-style product filter: full-length primer alignment, "
                f"{MIN_LEN}-{MAX_LEN} bp, ≤{MAX_MM_TOTAL} total mismatches, "
                f"perfect {MM_3PRIME_NEED}-nt 3' terminal match.\n")
        f.write("No network access; the six curated transcripts were mapped to contigs by "
                "megablast (≥95 % identity) to label expected bands.\n\n")

        f.write("## 1. Curated DQD/SDH transcripts → assembly contigs\n\n")
        f.write("| Curated transcript | Best matching contig | Identity | Aln length |\n")
        f.write("|---|---|---|---|\n")
        for q in tg:
            best = sorted(selfmap.get(q[0], []), key=lambda x: -x[2])[:1]
            if best:
                f.write(f"| {q[0]} | {best[0][0]} | {best[0][1]:.2f} % | {best[0][2]} bp |\n")
            else:
                f.write(f"| {q[0]} | *no ≥95 % contig* | – | – |\n")
        f.write("\n")

        f.write("## 2. Verdicts per primer pair\n\n")
        f.write("| Pair | Expected band on own contig | Other TARGET_CONTIG products | "
                "SDH-family contig products | Unrelated-contig products |\n")
        f.write("|---|---|---|---|---|\n")
        for gene in sorted(pairs):
            rs = [r for r in results if r['pair_gene'] == gene]
            tgt = [r for r in rs if r['class'] == 'TARGET_CONTIG']
            own = [r for r in tgt if r['contig'] in own_contigs[gene]]
            other_tgt = [r for r in tgt if r not in own]
            fam = [r for r in rs if r['class'] == 'SDH_FAMILY_CONTIG']
            unl = [r for r in rs if r['class'] == 'UNRELATED_CONTIG']
            fmt = lambda L: ', '.join(f"{r['contig']} {r['size']}bp" for r in L[:6]) or '0'
            f.write(f"| {gene} | {fmt(own)} | {fmt(other_tgt)} | {fmt(fam)} | {fmt(unl)} |\n")
        f.write("\n")

        f.write("## 3. All candidate products (detail)\n\n")
        if results:
            f.write("| Pair | Contig | Size | mm F | mm R | Class | Annotation product |\n")
            f.write("|---|---|---|---|---|---|---|\n")
            for r in sorted(results, key=lambda x: (x['class'], x['pair_gene'], x['size'])):
                f.write(f"| {r['pair_gene']} | {r['contig']} | {r['size']} | {r['mm_fwd']} | "
                        f"{r['mm_rev']} | {r['class']} | {r['product'] or r['gene_name'] or '-'} |\n")
        else:
            f.write("_No candidate products other than the six expected amplicons._\n")
        f.write("\n")

        f.write("## 4. Strong single-primer priming sites transcriptome-wide\n\n")
        f.write("A site = full-length 20/20 alignment, ≤3 mismatches, perfect 5-nt 3' end. "
                "Sites on the gene's own contig are expected; anything else is a partial-"
                "match site that needs a partner primer to matter.\n\n")
        f.write("| Primer | Strong sites total | on target contigs |\n|---|---|---|\n")
        for n, s in primers:
            st = [r for r in strong if r['primer'] == n]
            on_t = [r for r in st if r['tx'] in target_contigs]
            f.write(f"| {n} | {len(st)} | {len(on_t)} |\n")
        f.write("\n")

        f.write("## 5. Interpretation limits\n\n")
        f.write("- The search space is the **transcriptome** (assembled transcripts), which "
                "is the right template for RT-PCR / RT-qPCR from cDNA. Genomic DNA, "
                "un-assembled low-coverage transcripts and intronic regions are not covered.\n")
        f.write("- The assembly is clustered to longest isoforms, so allelic/paralogous "
                "transcripts collapsed out of it are under-represented; the six SDH paralogs "
                "were therefore mapped by megablast instead of assumed absent.\n")

    print("\n=== Candidate products (transcriptome-wide) ===")
    for gene in sorted(pairs):
        rs = [r for r in results if r['pair_gene'] == gene]
        own = [r for r in rs if r['class'] == 'TARGET_CONTIG']
        fam = [r for r in rs if r['class'] == 'SDH_FAMILY_CONTIG']
        unl = [r for r in rs if r['class'] == 'UNRELATED_CONTIG']
        print(f"  {gene}: target-contig {len(own)} | SDH-family {len(fam)} | unrelated {len(unl)}")
        for r in fam + unl:
            print(f"     -> {r['contig']} {r['size']} bp ({r['mm_total']} mm) : {r['product'][:70]}")
    print("\nWrote TRANSCRIPTOME_CHECK.md, tx_product_check.tsv, tx_priming_sites.tsv, tx_blast_hits.tsv")


if __name__ == '__main__':
    main()