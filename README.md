# transcript_counts_Salmon
Create PCA of read counts using Salmon and DESeq2  
  Salmon: https://github.com/COMBINE-lab/salmon  DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html


# 1. Index the transcriptomes with Salmon
```bash
sbatch -n 14 --mem 96G -t 4:00:00 -p Lewis,BioCompute --wrap="salmon index -t pan_trans_cds.fa -i pan_trans_cds_index -k 31 --keepDuplicates -p 14"
```
