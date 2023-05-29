
library(GenomicFeatures)
library(readr)
library(tximport)

# read the txdb to map the transcripts to the genes
txdb = loadDb(file = 'txdb_mmus.sqlite')
k <- keys(txdb, keytype="TXNAME")
tx_map <- AnnotationDbi::select(txdb, k, columns=c("TXNAME","GENEID"), "TXNAME")
write.csv(tx_map,file="txgene.csv",row.names=F)

tx2gene <- read_csv("txgene.csv")

#read the quant.sf directories
dirs <- list.files("data/")
quant_files <- list.files(path = "data/", full.names = TRUE, recursive = TRUE, pattern="quant.sf")
names(quant_files) <- dirs

#import the gene-level mappings using tximport
txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

write.csv(txi$abundance,file="txabnd.csv",row.names=T)
og_df <- data.frame(read_csv("txabnd.csv"))
colnames(og_df)[1] <- "ensembl_gene_id"
colnames(og_df)[2] <- "KO_expression"
colnames(og_df)[3] <- "WT_expression"

# map the gene ids to the gene-symbols
g_sym <- read_csv("g_sym.csv")
og_df <- data.frame(merge(og_df, g_sym, by = "ensembl_gene_id"))

# specify the conditions to build the dds
condition = as.factor(c("KO","WT"))
samples = data.frame(row.names = c("F2_KO_S23","F7_WT_S22"),condition)

#build the dds
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
vsd <- vst(dds)
lfc_df <- data.frame(vsd@assays@data@listData)
library(dplyr)
lfc_df <- dplyr::as_tibble(lfc_df, rownames = "ensembl_gene_id")

# calculate lfc
lfc_df$lfc <- lfc_df$F7_WT_S22 - lfc_df$F2_KO_S23
lfc_df$F2_KO_S23 <- NULL
lfc_df$F7_WT_S22 <- NULL

og_df <- data.frame(merge(og_df, lfc_df, by = "ensembl_gene_id"))
og_df$ensembl_gene_id <- NULL
og_df <- og_df[, c("gene_symbols", "lfc", "WT_expression", "KO_expression")]

# extract profilin1
pfn <- og_df[og_df$gene_symbols=="Pfn1", ]

# filter the tpms by >= 1
og_df <- og_df[og_df$WT_expression>=1 & og_df$KO_expression>=1,]

# find out the top 10 upregulated and downregulated genes
up <- og_df %>% slice_max(lfc, n = 10)

# add profilin to the up dataframe
up[nrow(up) + 1, ] <- pfn
#reset the index
rownames(up) <- NULL
write.csv(up,file="up_regulated_table.csv",row.names=F)

down <- og_df %>% slice_min(lfc, n = 10)
write.csv(down,file="down_regulated_table.csv",row.names=F)

# stack up the dataframe tpms for dot plot
ggup <- data.frame(up[1:2], stack(up[3:4]))
colnames(ggup)[3] <- "tpm_expression"
colnames(ggup)[4] <- "sample_type"

ggdown <- data.frame(down[1:2], stack(down[3:4]))
colnames(ggdown)[3] <- "tpm_expression"
colnames(ggdown)[4] <- "sample_type"

# plot the dot plot
library(ggplot2)
ggplot(ggup, aes(x=gene_symbols, y=tpm_expression, fill=sample_type)) + geom_dotplot(binaxis='y', stackdir='center') + facet_wrap(~reorder(gene_symbols,lfc,decreasing=T), scales="free", ncol=3) + xlab("genes") + ggtitle("Up-Regulated Genes")
ggplot(ggdown, aes(x=gene_symbols, y=tpm_expression, fill=sample_type)) + geom_dotplot(binaxis='y', stackdir='center') + facet_wrap(~reorder(gene_symbols,lfc,), scales="free", ncol=3) + xlab("genes") + ggtitle("Down-Regulated Genes")

