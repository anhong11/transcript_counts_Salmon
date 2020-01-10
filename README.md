# transcript_counts_Salmon
Create PCA of read counts using Salmon and DESeq2  
>Salmon: https://github.com/COMBINE-lab/salmon  
>DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html


## 1. Index the transcriptomes with Salmon
```bash
sbatch -n 14 --mem 96G -t 4:00:00 -p Lewis,BioCompute --wrap="salmon index -t pan_trans_cds.fa -i pan_trans_cds_index -k 31 --keepDuplicates -p 14"
```
## 2. Run salmon quant to get read counts  
run salmon_count.sh as:
```bash
sh salmon_count.sh sample.txt
```
>if the PE fq file named as 2109.pair_1.fq.gz, 2109.pair_2.fq.gz; 2110.pair_1.fq.gz, 2110.pair_2.fq.gz..., the `sample.txt` will be:
```bash
2109
2110
```
>and the `salmon_count.sh` is:
```bash
while read line
do
	sbatch -n 14 --mem 96G -t 2:00:00 -p Lewis,BioCompute --wrap="salmon quant -l A -i ../pan_trans_cds_index -1 '$line'.pair_1.fq.gz -2 '$line'.pair_2.fq.gz -o quants/'$line'_quant -p 14 --seqBias --validateMappings"
done <$1 
```
split *B. napus* quant.sf into A- and C- subgenomes, run split_AC_quant.sh as:
```bash
sh split_AC_quant.sh sample.txt
```
>and the `split_AC_quant.sh` is:
```bash
while read line
do
	cd "$line"_quant
	head -n1 quant.sf > head.txt
	awk '{if($1 ~ /^Cab/ || $1 ~ /^BnaA/)print}' quant.sf | cat head.txt - > ../chrA_"$line"_quant.sf
	awk '{if($1 ~ /^Bo/ || $1 ~ /^BnaC/)print}' quant.sf | cat head.txt - > ../chrC_"$line"_quant.sf
	cd ..
done<$1
```
## 3. Use tximport and DESeq2 to get PCA
Prepare packages for tximport and DESeq2
```bash
module load r/r-3.6.0-python-2.7.14-tk
export R_LIBS=~/R/:${R_LIBS} #~/R/: your own R library path
srun --pty -p Interactive --mem 64G R
```
Install R packages **only at the first time**, and then load them for R (in R cmd model)
```bash
install.packages("BiocManager")
install.packages("jsonlite")
install.packages("plyr", lib="~/R/")
install.packages("ggplot2")
install.packages("Hmisc")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
install.packages("RColorBrewer")
BiocManager::install("apeglm")

library(ggplot2)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(apeglm)
```

