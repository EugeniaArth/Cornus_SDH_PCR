# data/external/ — large inputs, not committed

Everything in this folder is gitignored. The pipeline reads these files through
`scripts/paths.py` (`P.ext(...)`), which looks for exactly the names below.

| Expected name | What it is | Where it comes from |
|---|---|---|
| `transcriptome.fasta` | Assembly of 29 081 transcripts, 48.1 Mb (`final.clust_transcripts_longest_iso.fasta`) | Trinity output of the *Cornus sericea* RNA-seq |
| `transcriptome_annotation.gff3` | Matching annotation (`final.clust_annotation_longest_iso.gff3`) | same assembly run |
| `annotated_final_nr.gff3` | BLASTx-vs-NR annotation, used to label contigs as SDH-family | annotation pipeline |
| `blastdb/` | BLAST nucleotide databases | created on first run by the scripts (`makeblastdb`) |

In the working copy these are symlinks to the files that stay on the lab disk, one
directory up (`PCR_SDH/Transcriptome/`, `PCR_SDH/annotated_final_nr.gff3`). They are
0-byte-looking in `ls` because they are links, not copies.

To use this repository elsewhere, either copy the real files in under these names, or
re-create the links:

```sh
ln -s ../../../Transcriptome/final.clust_transcripts_longest_iso.fasta  data/external/transcriptome.fasta
ln -s ../../../Transcriptome/final.clust_annotation_longest_iso.gff3   data/external/transcriptome_annotation.gff3
ln -s ../../../annotated_final_nr.gff3                                 data/external/annotated_final_nr.gff3
```

Checks 1 and 4, and all of stages 1-3, run without anything in this folder. Checks 2, 3 and 3b need
it, as does stage 5: `verify_sdh3_backup.py` rescans the assembly unconditionally, while
`design_sdh3_backup.py` skips its transcriptome verification when the file is absent. Those scripts
stop with a clear message rather than producing a partial result.
