#!/usr/bin/env python3
import os, sys, subprocess, glob, re, json
from collections import defaultdict

# Paths
TRANSCRIPTOME = 'PCR_Cornus_SDH/data/external/transcriptome.fasta'
AGROSTEMMA_DIR = 'PCR_Cornus_SDH/data/housekeeping/agrostemma_seqs'
BLAST_DIR = 'PCR_Cornus_SDH/data/housekeeping/agrostemma_blast'
KEYWORD_DIR = 'PCR_Cornus_SDH/data/housekeeping/keyword_search'

os.makedirs(BLAST_DIR, exist_ok=True)

def build_blast_db():
    """Build nucleotide BLAST database from transcriptome if not exists"""
    db_name = os.path.join(BLAST_DIR, 'cornus_transcriptome')
    if not os.path.exists(db_name + '.nsq'):
        print(f"Building BLAST database from {TRANSCRIPTOME}...")
        cmd = ['makeblastdb', '-in', TRANSCRIPTOME, '-dbtype', 'nucl', '-out', db_name]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error building BLAST database: {result.stderr}")
            sys.exit(1)
        print("BLAST database created.")
    else:
        print("BLAST database already exists.")
    return db_name

def run_blast_queries(db_path):
    """Run blastn for each Agrostemma FASTA file"""
    fasta_files = glob.glob(os.path.join(AGROSTEMMA_DIR, '*.fasta'))
    fasta_files = [f for f in fasta_files if not f.endswith('ALL_AGROSTEMMA.fasta')]

    all_results = {}

    for fasta in fasta_files:
        family = os.path.basename(fasta).replace('.fasta', '')
        print(f"BLASTing {family}...")
        outfile = os.path.join(BLAST_DIR, f'{family}_blast.tsv')

        # Run blastn
        cmd = [
            'blastn',
            '-query', fasta,
            '-db', db_path,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-evalue', '1e-5',
            '-num_threads', '4',
            '-max_target_seqs', '5',  # get top 5 hits per query
            '-out', outfile
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Warning: BLAST for {family} failed: {result.stderr}")
            continue

        # Parse results
        hits = parse_blast_output(outfile)
        all_results[family] = hits

    return all_results

def parse_blast_output(outfile):
    """Parse BLAST tabular output, return list of hits"""
    hits = []
    if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:
        return hits

    with open(outfile, 'r') as f:
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    hit = {
                        'query': fields[0],
                        'subject': fields[1],
                        'pident': float(fields[2]),
                        'length': int(fields[3]),
                        'evalue': float(fields[10]),
                        'bitscore': float(fields[11])
                    }
                    hits.append(hit)
    return hits

def group_hits_by_family(blast_results):
    """Group hits by family and select best hit per subject contig"""
    family_to_subjects = defaultdict(list)

    for family, hits in blast_results.items():
        # Sort hits by bitscore (descending)
        hits_sorted = sorted(hits, key=lambda x: x['bitscore'], reverse=True)
        # Take unique subject contigs (first occurrence = best bitscore)
        seen_subjects = set()
        for hit in hits_sorted:
            subj = hit['subject']
            if subj not in seen_subjects:
                seen_subjects.add(subj)
                family_to_subjects[family].append({
                    'contig': subj,
                    'pident': hit['pident'],
                    'bitscore': hit['bitscore'],
                    'evalue': hit['evalue']
                })
        print(f"  {family}: {len(family_to_subjects[family])} unique Cornus contigs")

    return family_to_subjects

def load_keyword_contigs():
    """Load contigs from keyword search results"""
    keyword_files = glob.glob(os.path.join(KEYWORD_DIR, '*_contigs.txt'))
    keyword_map = {}
    for f in keyword_files:
        family = os.path.basename(f).replace('_contigs.txt', '')
        with open(f, 'r') as fp:
            contigs = [line.strip() for line in fp if line.strip()]
            keyword_map[family] = set(contigs)
    return keyword_map

def compare_methods(blast_contigs, keyword_contigs):
    """Compare BLAST hits with keyword search results"""
    comparison = {}

    all_families = set(list(blast_contigs.keys()) + list(keyword_contigs.keys()))

    for fam in all_families:
        blast_set = set([c['contig'] for c in blast_contigs.get(fam, [])])
        keyword_set = keyword_contigs.get(fam, set())

        union = blast_set | keyword_set
        intersection = blast_set & keyword_set
        blast_only = blast_set - keyword_set
        keyword_only = keyword_set - blast_set

        comparison[fam] = {
            'blast_count': len(blast_set),
            'keyword_count': len(keyword_set),
            'intersection_count': len(intersection),
            'union_count': len(union),
            'blast_only': sorted(blast_only),
            'keyword_only': sorted(keyword_only),
            'intersection': sorted(intersection)
        }

    return comparison

def write_results(blast_contigs, comparison):
    """Write results to files"""
    # Write BLAST hits per family
    with open(os.path.join(BLAST_DIR, 'blast_hits_summary.tsv'), 'w') as f:
        f.write("Family\tContig\tPident\tBitscore\tEvalue\n")
        for fam, hits in blast_contigs.items():
            for h in hits:
                f.write(f"{fam}\t{h['contig']}\t{h['pident']:.2f}\t{h['bitscore']:.1f}\t{h['evalue']:.2e}\n")

    # Write comparison
    with open(os.path.join(BLAST_DIR, 'method_comparison.tsv'), 'w') as f:
        f.write("Family\tBlastHits\tKeywordHits\tIntersection\tUnion\tBlastOnly\tKeywordOnly\n")
        for fam in sorted(comparison.keys()):
            c = comparison[fam]
            f.write(f"{fam}\t{c['blast_count']}\t{c['keyword_count']}\t{c['intersection_count']}\t{c['union_count']}\t")
            f.write(f"{len(c['blast_only'])}\t{len(c['keyword_only'])}\n")

    # Write detailed comparison per family
    os.makedirs(os.path.join(BLAST_DIR, 'comparison_details'), exist_ok=True)
    for fam in comparison.keys():
        c = comparison[fam]
        with open(os.path.join(BLAST_DIR, 'comparison_details', f'{fam}_comparison.txt'), 'w') as f:
            f.write(f"Family: {fam}\n")
            f.write(f"BLAST hits: {c['blast_count']}\n")
            f.write(f"Keyword hits: {c['keyword_count']}\n")
            f.write(f"Intersection: {c['intersection_count']}\n")
            f.write(f"Union: {c['union_count']}\n")
            f.write("\n=== Contigs in intersection ===\n")
            for contig in c['intersection']:
                f.write(f"{contig}\n")
            f.write("\n=== Contigs only in BLAST ===\n")
            for contig in c['blast_only']:
                f.write(f"{contig}\n")
            f.write("\n=== Contigs only in keyword search ===\n")
            for contig in c['keyword_only']:
                f.write(f"{contig}\n")

    print(f"\nResults written to {BLAST_DIR}/")

def main():
    print("=== Step 2: BLAST Agrostemma sequences against Cornus transcriptome ===")

    # 1. Build BLAST database
    db_path = build_blast_db()

    # 2. Run BLAST queries
    blast_results = run_blast_queries(db_path)

    # 3. Group hits by family
    blast_contigs = group_hits_by_family(blast_results)

    # 4. Load keyword search results
    keyword_contigs = load_keyword_contigs()

    # 5. Compare methods
    comparison = compare_methods(blast_contigs, keyword_contigs)

    # 6. Write results
    write_results(blast_contigs, comparison)

    # Print summary
    print("\n=== Comparison Summary ===")
    for fam in sorted(comparison.keys()):
        c = comparison[fam]
        print(f"{fam}: BLAST={c['blast_count']}, Keyword={c['keyword_count']}, Intersection={c['intersection_count']}")

if __name__ == '__main__':
    main()
