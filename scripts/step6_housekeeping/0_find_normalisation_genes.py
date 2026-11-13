#!/usr/bin/env python3
"""Find C. sericea normalisation (reference) gene candidates in the annotated transcriptome.

Targets: GAPDH (glyceraldehyde-3-phosphate dehydrogenase), EF-1alpha (eukaryotic
translation elongation factor 1 alpha) and polyubiquitin (UBQ).

The script scans ../PCR_Cornus_SDH/data/external/transcriptome.gff3 for gene features whose product matches
each target, classifies the match, then extracts every candidate contig's nucleotide
sequence from ../PCR_Cornus_SDH/data/external/transcriptome.fasta in one pass.

Inputs are read-only. Outputs are written only inside this Normalisation/ folder:
  Normalisation_candidates.tsv   one row per candidate, grouped by family
  GAPDH_candidates.fasta         sequences (GAPDH family)
  EF1a_candidates.fasta          sequences (EF-1alpha)
  UBQ_candidates.fasta           sequences (polyubiquitin)

Nothing outside this folder is touched. Re-run safe (overwrites its own outputs).
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # PCR_SDH/
GFF3 = os.path.join(ROOT, "transcriptome.gff3")
TX   = os.path.join(ROOT, "transcriptome.fasta")

# Minimum a candidate must look real to be considered the "primary" reference gene.
MIN_LEN_FOR_PRIMARY = 1000
MIN_COV_FOR_PRIMARY = 50.0


def attr(attribs, key):
    m = re.search(rf"{key}=([^;]+)", attribs)
    return m.group(1) if m else ""


def has(seq, *subs):
    return any(s in seq for s in subs)


# ---------------------------------------------------------------- target matchers
def match_gapdh(product):
    if not re.search(r"glyceraldehyde[- ]3[- ]phosphate[- ]dehydrogenase", product, re.I):
        return None
    p = product.lower()
    if "nadp" in p:
        return "NADP-GAPDH (atypical, EC 1.1.1.354) - NOT reference"
    if "chloroplastic" in p:
        return "GAPDH, chloroplastic (EC 1.2.1.12)"
    return "GAPDH, cytosolic (EC 1.2.1.12) - reference"


def match_ef1a(product):
    p = product.lower()
    is_ef1a = ("elongation factor 1-alpha" in p or "elongation factor 1 alpha" in p
               or "elongation factor 1-alpha" in p)
    if not is_ef1a:
        return None
    if "beta" in p or "gamma" in p:            # EF-1b / EF-1g are different subunits
        return None
    return "EF-1alpha (eEF1A) - reference"


def match_ubq(product):
    if re.fullmatch(r"(poly)?ubiquitin", product, re.I):   # exact Ubiquitin / Polyubiquitin
        return "Polyubiquitin (UBQ) - reference"
    return None


def match_ubq_ext(product):
    p = product.lower()
    if re.search(r"ubiquitin[- ](ribosomal|extension|fusion)", p) or "extension protein" in p:
        return "Ubiquitin extension/fusion (eS31/eL40) - UBQ7-type"
    return None


TARGETS = [
    ("GAPDH",  match_gapdh,    "GAPDH_candidates.fasta"),
    ("EF1a",   match_ef1a,     "EF1a_candidates.fasta"),
    ("UBQ",    match_ubq,      "UBQ_candidates.fasta"),
    ("UBQext", match_ubq_ext,  "UBQext_candidates.fasta"),
]

# Rank order within GAPDH so the cytosolic reference comes first.
GAPDH_ORDER = {
    "GAPDH, cytosolic (EC 1.2.1.12) - reference": 0,
    "GAPDH, chloroplastic (EC 1.2.1.12)": 1,
    "NADP-GAPDH (atypical, EC 1.1.1.354) - NOT reference": 2,
}


def scan_gff3():
    """family -> list of candidate dicts (one per gene locus)."""
    buckets = {fam: {} for fam, _, _ in TARGETS}
    with open(GFF3) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            seqid, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
            a = p[8]
            product = attr(a, "product")
            for fam, matcher, _ in TARGETS:
                cls = matcher(product)
                if cls is None:
                    continue
                if seqid in buckets[fam]:
                    continue
                buckets[fam][seqid] = {
                    "contig": seqid,
                    "name": attr(a, "Name") or "-",
                    "product": product,
                    "uniprot_hit": attr(a, "UniProt_best_hit").lstrip("="),
                    "class": cls,
                    "strand": strand,
                    "tx_len": int(re.search(r"_length_(\d+)_", seqid).group(1)),
                    "cov": float(re.search(r"_cov_([\d.]+)_", seqid).group(1)),
                }
    return {fam: list(v.values()) for fam, v in buckets.items()}


def extract_sequences(wanted):
    seqs, keep, buf, cid = {}, False, [], None
    with open(TX) as f:
        for line in f:
            if line[0] == ">":
                if keep and cid:
                    seqs[cid] = "".join(buf)
                cid = line[1:].split()[0]
                keep = cid in wanted
                buf = []
            elif keep:
                buf.append(line.strip())
    if keep and cid:
        seqs[cid] = "".join(buf)
    return seqs


def sort_family(fam, rows):
    if fam == "GAPDH":
        return sorted(rows, key=lambda r: (GAPDH_ORDER.get(r["class"], 9), -r["cov"]))
    return sorted(rows, key=lambda r: -r["cov"])


def primary(rows):
    ok = [r for r in rows if "NOT reference" not in r["class"]]
    pool = [r for r in ok if r["tx_len"] >= MIN_LEN_FOR_PRIMARY and r["cov"] >= MIN_COV_FOR_PRIMARY]
    if pool:
        return max(pool, key=lambda r: r["cov"])
    # fallback for short families (e.g. ubiquitin-fusion ORFs are ~0.8 kb): best expressed
    pool = [r for r in ok if r["cov"] >= MIN_COV_FOR_PRIMARY]
    return max(pool, key=lambda r: r["cov"]) if pool else None


def write_fasta(path, rows, seqs):
    with open(path, "w") as f:
        for r in rows:
            s = seqs.get(r["contig"], "")
            f.write(f">{r['contig']} {r['name']} | {r['class']} | cov={r['cov']:.1f} len={len(s)}\n")
            for j in range(0, len(s), 60):
                f.write(s[j:j+60] + "\n")


def write_family_tsv(path, fam, rows):
    with open(path, "w") as f:
        f.write("family\trank\tcontig\tname\tclass\ttx_len\tcov\tstrand\tuniprot_hit\tproduct\n")
        for i, r in enumerate(rows, 1):
            f.write(f"{fam}\t{i}\t{r['contig']}\t{r['name']}\t{r['class']}\t"
                    f"{r['tx_len']}\t{r['cov']:.1f}\t{r['strand']}\t{r['uniprot_hit']}\t{r['product']}\n")


def main():
    buckets = scan_gff3()
    all_ids = {r["contig"] for rows in buckets.values() for r in rows}
    seqs = extract_sequences(all_ids)

    order = ["GAPDH", "EF1a", "UBQ", "UBQext"]
    combined = os.path.join(HERE, "Normalisation_candidates.tsv")
    with open(combined, "w") as f:
        f.write("family\trank\tcontig\tname\tclass\ttx_len\tcov\tstrand\tuniprot_hit\tproduct\n")
        for fam in order:
            rows = sort_family(fam, buckets[fam])
            for i, r in enumerate(rows, 1):
                f.write(f"{fam}\t{i}\t{r['contig']}\t{r['name']}\t{r['class']}\t"
                        f"{r['tx_len']}\t{r['cov']:.1f}\t{r['strand']}\t{r['uniprot_hit']}\t{r['product']}\n")

    for fam in order:
        rows = sort_family(fam, buckets[fam])
        fasta_name = dict((t[0], t[2]) for t in TARGETS)[fam]
        write_fasta(os.path.join(HERE, fasta_name), rows, seqs)
        write_family_tsv(os.path.join(HERE, fasta_name.replace(".fasta", ".tsv")), fam, rows)

    # summary
    for fam in order:
        rows = sort_family(fam, buckets[fam])
        p = primary(rows)
        print(f"\n== {fam}: {len(rows)} locus/loci ==")
        for i, r in enumerate(rows, 1):
            star = "  <= PRIMARY" if p and r["contig"] == p["contig"] else ""
            print(f"  {i:2d}. {r['name']:10s} {r['contig'][:40]:40s} "
                  f"len={r['tx_len']:5d} cov={r['cov']:9.1f} {r['class']}{star}")
    print("\nWrote Normalisation_candidates.tsv + per-family *_candidates.tsv/.fasta "
          "(GAPDH, EF1a, UBQ, UBQext)")


if __name__ == "__main__":
    main()
