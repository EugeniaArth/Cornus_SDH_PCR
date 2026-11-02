#!/usr/bin/env python3
"""
Primer-BLAST-style specificity check of the final Cornus DQD/SDH primer set,
run LOCALLY with blastn (no data leaves the project folder).

Checks performed per the NCBI Primer-BLAST logic (simplified):
  1. blastn-short every primer against the six Cornus DQD/SDH transcripts.
  2. Count priming sites per transcript with relaxed identity (allows mismatches),
     like Primer-BLAST's "Select primers inside transcripts" + specificity screen.
  3. Pair forward/reverse hits that face each other within 80-600 bp and score
     potential mis-priming products (total 3'-end mismatches <= 3, similar to
     Primer-BLAST's default 'In a product, max mismatches' = 3, 3' 5 nt must match).
  4. Report per-pair verdict: OK / AMBIGUOUS / NON-SPECIFIC.

Outputs (all inside primer_blast_check/, existing results untouched):
  - blast_hits.tsv        : raw blastn hits
  - product_check.tsv     : candidate amplicons per pair x transcript
  - PRIMER_BLAST_CHECK.md : human-readable verdict table
"""
import csv, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
DB = os.path.join(HERE, 'cornus_dqdsdh_6tx')
PRIMER_FASTA = os.path.join(PARENT, 'FINAL_primers.fasta')

MIN_LEN, MAX_LEN = 80, 600      # product window (like Primer-BLAST default range)
MAX_MM_TOTAL = 3                # total mismatches allowed in a product (both primers)
MM_3PRIME_LAST5_MAX = 0         # last 5 nt of the 3' end must match perfectly

def load_primers(path):
    primers, name, buf = {}, None, []
    seq_all = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    seq_all.append((name, ''.join(buf)))
                name = line[1:]
                buf = []
            else:
                buf.append(line)
        if name:
            seq_all.append((name, ''.join(buf)))
    return seq_all

def blast(queries, db, outfmt_cols, extra=()):
    cmd = ['blastn', '-task', 'blastn-short', '-word_size', '7',
           '-evalue', '1000', '-perc_identity', '75',
           '-max_target_seqs', '50', '-db', db,
           '-outfmt', f'6 {outfmt_cols}'] + list(extra)
    inp = '\n'.join(f'>{n}\n{s}' for n, s in queries)
    res = subprocess.run(cmd, input=inp, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr + '\n')
        sys.exit(1)
    return res.stdout

def main():
    primers = load_primers(PRIMER_FASTA)
    print(f"Loaded {len(primers)} primers from FINAL_primers.fasta")

    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    raw = blast(primers, DB, cols)

    hits = []
    with open(os.path.join(HERE, 'blast_hits.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['primer','transcript','pident','alen','mismatch','qstart','qend',
                    'sstart','send','strand','bitscore'])
        for line in raw.strip().splitlines():
            if not line: continue
            p = line.split('\t')
            rec = {
                'primer': p[0], 'transcript': p[1], 'pident': float(p[2]),
                'alen': int(p[3]), 'mismatch': int(p[4]), 'qstart': int(p[6]),
                'qend': int(p[7]), 'sstart': int(p[8]), 'send': int(p[9]),
                'sstrand': 'plus' if int(p[8]) < int(p[9]) else 'minus',
                'bits': float(p[11]),
            }
            hits.append(rec)
            w.writerow([rec['primer'], rec['transcript'], rec['pident'], rec['alen'],
                        rec['mismatch'], rec['qstart'], rec['qend'],
                        rec['sstart'], rec['send'], rec['sstrand'], rec['bits']])

    print(f"blastn-short hits: {len(hits)}")

    # ---- group by gene pair ----
    # primer names look like: CserDQD_SDH1.1_Fwd
    pairs = {}
    for name, seq in primers:
        gene = name.rsplit('_', 1)[0]
        side = name.rsplit('_', 1)[1]
        pairs.setdefault(gene, {})[side] = {'name': name, 'seq': seq}

    # ---- product search allowing mismatches ----
    def hits_for(primer_name):
        return [h for h in hits if h['primer'] == primer_name]

    def three_prime_clean(h, qseq):
        """Require perfect match on the last 5 bases of the primer at the template 3' end."""
        # qstart/qend are primer coordinates in query orientation; determine which
        # end of the primer contacts template at the extension side.
        # For a plus-strand genomic hit (sstart<send) and primer in forward orientation,
        # the 3' end of the primer aligns to qend / send.
        if h['sstrand'] == 'plus':
            tail_q = qseq[-5:]
        else:
            tail_q = qseq[-5:]
        # We need the aligned subject segment corresponding to the primer 3' end.
        # Reconstruct from blast alignment via qstart/qend on stored qseq/sseq rows.
        return True  # placeholder; replaced below by direct alignment scan

    # For robustness, rebuild alignments from blast output with qseq/sseq fields.
    aln = []
    for line in raw.strip().splitlines():
        if not line: continue
        p = line.split('\t')
        aln.append({
            'primer': p[0], 'transcript': p[1], 'pident': float(p[2]),
            'length': int(p[3]), 'mismatch': int(p[4]),
            'qstart': int(p[6]), 'qend': int(p[7]),
            'sstart': int(p[8]), 'send': int(p[9]),
            'qseq': p[12], 'sseq': p[13],
        })

    def count_mismatches(rec):
        return rec['mismatch']

    def three_prime_matches(rec):
        """Number of perfectly matched bases counting back from the primer 3' end.
        Stop at the first mismatch."""
        q, s = rec['qseq'], rec['sseq']
        # 3' end of the primer is at qend (if primer is + on template) else qstart.
        # For extension, the polymerase extends from the primer 3'-OH, so the critical
        # bases are at the genomic side of the 3' end.
        # Determine alignment direction: if sstart<send, template plus-strand.
        # Query (primer) is always written 5'->3'. If template hit is minus strand,
        # then blast still aligned primer against minus strand; qstart<qend but
        # sstart>send. In that case the 3' end of the primer aligns at send.
        plus = rec['sstart'] < rec['send']
        # Walk the alignment from the appropriate end:
        pairs = list(zip(q, s))
        if plus:
            walk = reversed(pairs)     # 3' end of primer is at the end of alignment
        else:
            walk = pairs               # minus-strand template: alignment already reversed
            walk = reversed(pairs)     # blast writes qseq 5'->3', so 3' end still at the end
        good = 0
        for qb, sb in walk:
            if qb.upper() == sb.upper():
                good += 1
            else:
                break
        return good

    results = []
    for gene, d in sorted(pairs.items()):
        fwd, rev = d['Fwd'], d['Rev']
        plen_f = len(fwd['seq'])
        plen_r = len(rev['seq'])
        # Forward-strand products: Fwd hits (plus) followed by Rev hits (minus)
        for tx in sorted({a['transcript'] for a in aln}):
            # require FULL-LENGTH alignment of the primer (qstart=1, qend=primer len)
            f_hits = [a for a in aln if a['primer'] == fwd['name'] and a['transcript'] == tx
                      and a['sstart'] < a['send']
                      and a['qstart'] == 1 and a['qend'] == plen_f
                      and a['length'] == plen_f]
            r_hits = [a for a in aln if a['primer'] == rev['name'] and a['transcript'] == tx
                      and a['sstart'] > a['send']
                      and a['qstart'] == 1 and a['qend'] == plen_r
                      and a['length'] == plen_r]
            for fh in f_hits:
                for rh in r_hits:
                    start = fh['sstart']
                    end = rh['sstart']          # minus-strand hit: sstart is the right end
                    size = end - start + 1
                    if not (MIN_LEN <= size <= MAX_LEN):
                        continue
                    mm_f = count_mismatches(fh)
                    mm_r = count_mismatches(rh)
                    ok3_f = three_prime_matches(fh)
                    ok3_r = three_prime_matches(rh)
                    total_mm = mm_f + mm_r
                    if total_mm > MAX_MM_TOTAL:
                        continue
                    if ok3_f < 5 or ok3_r < 5:
                        continue
                    is_self = (tx == gene.replace('CserDQD_SDH', 'CserDQD/SDH'))
                    results.append({
                        'pair_gene': gene, 'transcript': tx,
                        'size': size, 'mm_fwd': mm_f, 'mm_rev': mm_r,
                        'mm_total': total_mm,
                        'match': 'TARGET' if is_self else 'OFF-TARGET',
                    })

    with open(os.path.join(HERE, 'product_check.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['pair_gene','transcript','size','mm_fwd','mm_rev','mm_total','match'],
                           delimiter='\t')
        w.writeheader()
        for r in sorted(results, key=lambda x: (x['pair_gene'], x['match'], x['size'])):
            w.writerow(r)

    # ---- per-primer priming-site summary ----
    print("\nPrimng sites per primer (any identity >=75%, e<1000):")
    per_primer = {}
    for a in aln:
        per_primer.setdefault(a['primer'], {}).setdefault(a['transcript'], 0)
        per_primer[a['primer']][a['transcript']] += 1

    verdicts = {}
    for gene, d in sorted(pairs.items()):
        offs = [r for r in results if r['pair_gene'] == gene and r['match'] == 'OFF-TARGET']
        tgts = [r for r in results if r['pair_gene'] == gene and r['match'] == 'TARGET']
        if not tgts:
            verdict = 'NO TARGET PRODUCT?'
        elif offs:
            verdict = 'AMBIGUOUS (off-target products possible with mismatches)'
        else:
            verdict = 'SPECIFIC'
        verdicts[gene] = (verdict, tgts, offs)

    # ---- write markdown report ----
    with open(os.path.join(HERE, 'PRIMER_BLAST_CHECK.md'), 'w') as f:
        f.write("# Primer-BLAST-style local check (blastn-short vs the six Cornus DQD/SDH transcripts)\n\n")
        f.write("Parameters: `-task blastn-short -word_size 7 -evalue 1000 -perc_identity 75`.  \n")
        f.write(f"Product window: {MIN_LEN}-{MAX_LEN} bp, total mismatches ≤ {MAX_MM_TOTAL}, "
                "perfect 5-nt 3'-terminal match required for each primer.\n")
        f.write("Database: `DQD_SDH_possible_nt.fasta` (6 transcripts, 12,888 bp).  \n")
        f.write("Note: this is a LOCAL approximation of Primer-BLAST specificity. "
                "It does not scan the whole transcriptome yet (see note at bottom).\n\n")
        f.write("## Verdicts\n\n")
        f.write("| Pair | Verdict | Target product(s) | Off-target product(s) |\n")
        f.write("|---|---|---|---|\n")
        for gene, (v, tgts, offs) in verdicts.items():
            tg = '; '.join(f"{r['size']} bp ({r['mm_total']} mm)" for r in tgts) or '-'
            of = '; '.join(f"{r['transcript']} {r['size']} bp ({r['mm_total']} mm)" for r in offs) or 'none'
            f.write(f"| {gene} | {v} | {tg} | {of} |\n")
        f.write("\n## Priming sites per primer (single-primer hits, any identity ≥75 %)\n\n")
        f.write("| Primer | Transcript | Hits |\n|---|---|---|\n")
        for pr, per_tx in sorted(per_primer.items()):
            for tx, n in sorted(per_tx.items()):
                f.write(f"| {pr} | {tx} | {n} |\n")
        f.write("\n## Full candidate products (allowing mismatches)\n\n")
        f.write("| Pair | Transcript | Size | mm(F) | mm(R) | mm total | Class |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in sorted(results, key=lambda x: (x['pair_gene'], x['size'])):
            f.write(f"| {r['pair_gene']} | {r['transcript']} | {r['size']} | "
                    f"{r['mm_fwd']} | {r['mm_rev']} | {r['mm_total']} | {r['match']} |\n")
        f.write("\n## Still missing for a real Primer-BLAST\n\n")
        f.write("- Whole-transcriptome database (the Trinity `transcripts.fasta` in the "
                "sibling `Cornus_transcriptome` folder) is outside this project directory — "
                "it needs explicit read permission.\n")
        f.write("- The *Cornus* genome/transcriptome is not in NCBI, so the web Primer-BLAST "
                "would require pasting the six transcripts as `primer-orientation` template "
                "or uploading a custom database.\n")

    # console summary
    print("\n=== Verdicts (local blastn Primer-BLAST-style) ===")
    for gene, (v, tgts, offs) in verdicts.items():
        print(f"  {gene}: {v}  (target products: {len(tgts)}, off-target: {len(offs)})")
        for r in offs:
            print(f"     ⚠ {r['transcript']} {r['size']} bp, mm_total={r['mm_total']}")

if __name__ == '__main__':
    main()