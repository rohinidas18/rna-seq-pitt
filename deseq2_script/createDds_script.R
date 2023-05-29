## Differential expression analysis using DESeq2

## Setup: Import libraries as required

library(GenomicFeatures)
library(readr)
library(tximport)
library(DESeq2)

## download the annotated gtf file, and write the transcipt name v/s geneid to .csv

gtf_file <- "Mus_musculus.GRCm39.106.chr.gtf.gz"
#file.exists(gtf_file)

txdb <- makeTxDbFromGFF(gtf_file)
#keytypes(txdb)
#columns(txdb)
 
k <- keys(txdb, keytype="TXNAME")
tx_map <- AnnotationDbi::select(txdb, k, columns=c("TXNAME","GENEID"), "TXNAME")
#head(tx_map)

write.csv(tx_map,file="txgene.csv",row.names=F)

## read the csv to a vector

tx2gene <- read_csv("txgene.csv")
#head(tx2gene)

## List all directories containing data and to the quant result files

dirs <- list.files("data/")
quant_files <- list.files(path = "data/", full.names = TRUE, recursive = TRUE, pattern="quant.sf")
names(quant_files) <- dirs
#quant_files

## use tximport to display the reads per gene

txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
#head(txi$counts)

condition = as.factor(c("KO","WT"))
samples = data.frame(row.names = c("F2_KO_S23","F7_WT_S22"),condition)
#samples

## contruct a DESeq2 dataset

dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
dds@assays@data@listData$counts
#dds <- DESeq(dds)
#res <- results(dds)
