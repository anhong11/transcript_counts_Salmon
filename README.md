# transcript_counts_Salmon
Create PCA of read counts using Salmon and DESeq2  
	Salmon: https://github.com/COMBINE-lab/salmon  
	DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html


## 1. Index the transcriptomes with Salmon
```bash
sbatch -n 14 --mem 96G -t 4:00:00 -p Lewis,BioCompute --wrap="salmon index -t pan_trans_cds.fa -i pan_trans_cds_index -k 31 --keepDuplicates -p 14"
```
## 2. Run salmon quant to get read counts  
run salmon_count.sh as:
```bash
sh salmon_count.sh sample.txt
```
	if the PE fq file named as 2109.pair_1.fq.gz, 2109.pair_2.fq.gz; 2110.pair_1.fq.gz, 2110.pair_2.fq.gz..., the sample.txt will be:
```bash
2109
2110
```
	and the salmon_count.sh is:
```bash
while read line
do
	sbatch -n 14 --mem 96G -t 2:00:00 -p Lewis,BioCompute --wrap="salmon quant -l A -i ../pan_trans_cds_index -1 '$line'.pair_1.fq.gz -2 '$line'.pair_2.fq.gz -o quants/'$line'_quant -p 14 --seqBias --validateMappings"
done <$1 
```
