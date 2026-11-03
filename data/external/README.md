# data/external/ — large inputs, not committed

Everything in this folder is gitignored. The pipeline reads these files through
`scripts/paths.py` (`P.ext(...)`), which looks for exactly the names below.

| Expected name | What it is | Where it comes from |
|---|---|---|
| `transcriptome.fasta` | Assembly of 29 081 transcripts | Trinity output of the *Cornus sericea* RNA-seq |
| `transcriptome_annotation.gff3` | Matching annotation | same assembly run |
| `blastdb/` | BLAST nucleotide databases | created on first run by the scripts (`makeblastdb`) |

Assembly and annotation are available at https://drive.google.com/drive/folders/1aDozkY-kDjnAsJzh3PdtklgKbpfJFE5D?usp=sharing.
