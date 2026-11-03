# Whole-transcriptome off-target check of the final primer set

Template: `Transcriptome/final.clust_transcripts_longest_iso.fasta` (29,124 transcripts)  
Engine: local `blastn -task blastn-short -word_size 7 -perc_identity 70`, Primer-BLAST-style product filter: full-length primer alignment, 80-600 bp, ≤3 total mismatches, perfect 5-nt 3' terminal match.
No network access; the six curated transcripts were mapped to contigs by megablast (≥95 % identity) to label expected bands.

## 1. Curated DQD/SDH transcripts → assembly contigs

| Curated transcript | Best matching contig | Identity | Aln length |
|---|---|---|---|
| CserDQD/SDH1.1 | NODE_10980_length_2779_cov_121.715385_g1243_i5 | 100.00 % | 1585 bp |
| CserDQD/SDH1.2 | NODE_10980_length_2779_cov_121.715385_g1243_i5 | 100.00 % | 2779 bp |
| CserDQD/SDH2 | NODE_18517_length_2241_cov_255.196168_g8755_i0 | 100.00 % | 2241 bp |
| CserDQD/SDH3 | NODE_20612_length_2133_cov_1535.833013_g8101_i1 | 100.00 % | 2133 bp |
| CserDQD/SDH4 | NODE_21648_length_2083_cov_1398.299410_g8244_i1 | 100.00 % | 2083 bp |
| CserDQD/SDH5 | NODE_29040_length_1767_cov_22.706636_g12409_i1 | 100.00 % | 1767 bp |

## 2. Verdicts per primer pair

| Pair | Expected band on own contig | Other TARGET_CONTIG products | SDH-family contig products | Unrelated-contig products |
|---|---|---|---|---|
| CserDQD_SDH1.1 | 0 | 0 | 0 | 0 |
| CserDQD_SDH1.2 | NODE_10980_length_2779_cov_121.715385_g1243_i5 180bp | 0 | 0 | 0 |
| CserDQD_SDH2 | NODE_18517_length_2241_cov_255.196168_g8755_i0 240bp | 0 | 0 | 0 |
| CserDQD_SDH3 | NODE_20612_length_2133_cov_1535.833013_g8101_i1 300bp | 0 | 0 | 0 |
| CserDQD_SDH4 | NODE_21648_length_2083_cov_1398.299410_g8244_i1 360bp | 0 | 0 | 0 |
| CserDQD_SDH5 | NODE_29040_length_1767_cov_22.706636_g12409_i1 420bp | 0 | 0 | 0 |

## 3. All candidate products (detail)

| Pair | Contig | Size | mm F | mm R | Class | Annotation product |
|---|---|---|---|---|---|---|
| CserDQD_SDH1.2 | NODE_10980_length_2779_cov_121.715385_g1243_i5 | 180 | 0 | 0 | TARGET_CONTIG | XP_059625762.1 bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase, chloroplastic-like isoform X3 [Cornus  |
| CserDQD_SDH2 | NODE_18517_length_2241_cov_255.196168_g8755_i0 | 240 | 0 | 0 | TARGET_CONTIG | XP_059632657.1 bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase, chloroplastic isoform X1 [Cornus flori |
| CserDQD_SDH3 | NODE_20612_length_2133_cov_1535.833013_g8101_i1 | 300 | 0 | 0 | TARGET_CONTIG | XP_059652600.1 bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase, chloroplastic-like [Cornus florida] |
| CserDQD_SDH4 | NODE_21648_length_2083_cov_1398.299410_g8244_i1 | 360 | 0 | 0 | TARGET_CONTIG | XP_059652599.1 bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase, chloroplastic-like [Cornus florida] |
| CserDQD_SDH5 | NODE_29040_length_1767_cov_22.706636_g12409_i1 | 420 | 0 | 0 | TARGET_CONTIG | XP_059652601.1 bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase, chloroplastic-like [Cornus florida] |

## 4. Strong single-primer priming sites transcriptome-wide

A site = full-length 20/20 alignment, ≤3 mismatches, perfect 5-nt 3' end. Sites on the gene's own contig are expected; anything else is a partial-match site that needs a partner primer to matter.

| Primer | Strong sites total | on target contigs |
|---|---|---|
| CserDQD_SDH1.1_Fwd | 0 | 0 |
| CserDQD_SDH1.1_Rev | 1 | 1 |
| CserDQD_SDH1.2_Fwd | 1 | 1 |
| CserDQD_SDH1.2_Rev | 1 | 1 |
| CserDQD_SDH2_Fwd | 1 | 1 |
| CserDQD_SDH2_Rev | 1 | 1 |
| CserDQD_SDH3_Fwd | 1 | 1 |
| CserDQD_SDH3_Rev | 1 | 1 |
| CserDQD_SDH4_Fwd | 1 | 1 |
| CserDQD_SDH4_Rev | 1 | 1 |
| CserDQD_SDH5_Fwd | 1 | 1 |
| CserDQD_SDH5_Rev | 1 | 1 |

## 5. Interpretation limits

- The search space is the **transcriptome** (assembled transcripts), which is the right template for RT-PCR / RT-qPCR from cDNA. Genomic DNA, un-assembled low-coverage transcripts and intronic regions are not covered.
- The assembly is clustered to longest isoforms, so allelic/paralogous transcripts collapsed out of it are under-represented; the six SDH paralogs were therefore mapped by megablast instead of assumed absent.
