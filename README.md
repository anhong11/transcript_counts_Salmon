# transcript_counts_Salmon
Create PCA of read counts using Salmon and DESeq2  

## 1. Index the transcriptomes with [Salmon](https://github.com/COMBINE-lab/salmon)
```bash
sbatch -n 14 --mem 96G -t 4:00:00 -p Lewis,BioCompute --wrap="salmon index -t pan_trans_cds.fa -i pan_trans_cds_index -k 31 --keepDuplicates -p 14"

For B. rapa:
sbatch -n 14 --mem 96G -t 4:00:00 -p Lewis,BioCompute --wrap="salmon index -t pan_A_cds.fa -i pan_A_cds_index -k 31 --keepDuplicates -p 14"
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
	###for B. rapa change "pan_trans_cds_index" to "pan_A_cds_index"
done <$1 
```
split *B. napus* quant.sf into A- and C- subgenomes, run split_AC_quant.sh as:
> !!only *B. napus* need this step
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
## 3. Use tximport and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to get PCA
Prepare packages for tximport and DESeq2
```bash
module load r/r-3.6.0-python-2.7.14-tk
export R_LIBS=~/R/:${R_LIBS} #~/R/: your own R library path
srun --pty -p Interactive -c8 --mem 64G R

###Install R packages only at the first time, and then load them for R (in R cmd model)###
install.packages("BiocManager")
install.packages("jsonlite")
install.packages("plyr", lib="~/R/")
install.packages("ggplot2")
install.packages("Hmisc")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
install.packages("RColorBrewer")
BiocManager::install("apeglm")
##########################################################################################

#Load library
library(ggplot2)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(apeglm)
```

B. napus only R script
```bash
#Calculate
samples <- read.table("phenotype_Bn_noAdmix.txt", header = TRUE)
files <- file.path("/storage/htc/pireslab/salmon/b_napus/quants", samples$sample, "quant.sf")
trans2gene <- read.csv("trans2gene_all")
txi <- tximport(files, type = "salmon", tx2gene = trans2gene)
write.csv(txi$counts, "txi.csv") #use this if you just want a table of counts per sample
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~1)
keep <- rowSums(counts(dds)) >= 2
dds <- dds[keep,]
colSums(counts(dds)) #this just checks how many reads per sample
dds <- estimateSizeFactors(dds)  #I think this is correcting for library size, which is important for all Brassica dataset
dds_vst <- vst(dds) #normalized with respect to library size or other normalization factors.
pcaData <- plotPCA(dds_vst, intgroup=c("group"), ntop = 500, returnData=TRUE)  #PCA plot; MMbary's cmd: pcaData <- plotPCA(dds_vst, intgroup=c("phenotype", "group"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Plot
cols <-c("WEAm"='yellow', "WEsA"='orange', "S"='red', "R"='blue', "SK"='purple', "WeA"='green')
PCA_Expression <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
  						geom_point(size=3) +
 						xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 						ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
 						scale_color_manual(values=cols) +
 						ggtitle("Expression Variation between different B. napus") +
 	 					coord_fixed()
ggsave("PCA_Expression_Bnapus.pdf", PCA_Expression, width = 8.5, height = 9)
```

A subgenome(B. rapa and B. napus) only R script\
```bash
#Calculate
samples <- read.table("phyenotype_A_noAdmix.txt", header = TRUE)
files <- file.path("/storage/htc/pireslab/salmon/sub_A_quant", samples$sample, "quant.sf")
trans2gene <- read.csv("trans2gene_A")
txi <- tximport(files, type = "salmon", tx2gene = trans2gene)
write.csv(txi$counts, "txi.csv") #use this if you just want a table of counts per sample
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~1)
keep <- rowSums(counts(dds)) >= 2
dds <- dds[keep,]
colSums(counts(dds)) #this just checks how many reads per sample
dds <- estimateSizeFactors(dds)  #I think this is correcting for library size, which is important for all Brassica dataset
dds_vst <- vst(dds) #normalized with respect to library size or other normalization factors.
pcaData <- plotPCA(dds_vst, intgroup=c("phenotype", "group"), ntop = 500, returnData=TRUE)  #PCA plot
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Plot
cols <-c("WEAm"='yellow', "WEsA"='orange', "S"='red', "R"='blue', "SK"='purple', "WeA"='green')
PCA_Expression <- ggplot(pcaData, aes(PC1, PC2, color=group, shape=phenotype.1)) +
  						geom_point(size=3) +
 						xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 						ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
 						scale_color_manual(values=cols) +
 						ggtitle("Expression Variation among A subgenome") +
 	 					coord_fixed()
ggsave("PCA_Expression_Asub.pdf", PCA_Expression, width = 8.5, height = 9)

```

