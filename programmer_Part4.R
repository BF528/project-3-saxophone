setwd("C:/Users/jason/Desktop/BU/BF528/project_3")
library(DESeq2)
library(apeglm)
library(ggplot2)
library(dplyr)
library(tidyverse)
#Slibrary(EnhancedVolcano)

CCounts<-read.csv(file = "control_counts.csv", header = TRUE, row.names=1)
sCCounts = subset(CCounts, select = c('SRR1178004', 'SRR1178006', 'SRR1178013', 'SRR1178064', 'SRR1178074', 'SRR1178075'))

conc_counts<-read.csv(file = "conc_counts.csv")
conc_control = cbind(conc_counts, sCCounts)

conc_control <- subset(conc_control,rowSums(conc_control==0)==0)

info <- read.csv('group_4_rna_info.csv', header = TRUE)

#Subsetting combined data to match info by sample and order
DNA_group <- subset(conc_control, select = c(2, 3, 4, 11, 12, 13, 14, 15, 16))
DNA_info <- info[-c(4, 5, 6, 7, 8, 9), ]

ER_group <- subset(conc_control, select = c(5, 6, 7, 11, 12, 13, 14, 15, 16))
ER_info <- info[-c(1, 2, 3, 7, 8, 9), ]

PPARA_group <- subset(conc_control, select = c(8, 9, 10, 11, 12, 13, 14, 15, 16))
PPARA_info <- info[-c(1, 2, 3, 4, 5, 6), ]

#DESeq objects
dds_DNA_group <- DESeqDataSetFromMatrix(
  countData = DNA_group,
  colData = DNA_info,
  design= ~ mode_of_action
)

dds_ER_group <- DESeqDataSetFromMatrix(
  countData = ER_group,
  colData = ER_info,
  design= ~ mode_of_action
)

dds_PPARA_group <- DESeqDataSetFromMatrix(
  countData = PPARA_group,
  colData = PPARA_info,
  design= ~ mode_of_action
)

#Relevel mode_of_action as factor
dds_DNA_group$mode_of_action <- relevel(dds_DNA_group$mode_of_action, ref='Control')
dds_ER_group$mode_of_action <- relevel(dds_ER_group$mode_of_action, ref='Control')
dds_PPARA_group$mode_of_action <- relevel(dds_PPARA_group$mode_of_action, ref='Control')

#Running DESeq
dds_DNA <- DESeq(dds_DNA_group)
dds_ER <- DESeq(dds_ER_group)
dds_PPARA <- DESeq(dds_PPARA_group)

#DESeq results
res_DNA <- results(dds_DNA, contrast=c('mode_of_action', 'DNA_Damage', 'Control'))
res_ER <- results(dds_ER, contrast=c('mode_of_action', 'ER', 'Control'))
res_PPARA <- results(dds_PPARA, contrast=c('mode_of_action', 'PPARA', 'Control'))

#Shrink
res_DNA <- lfcShrink(dds_DNA, coef=2)
res_ER <- lfcShrink(dds_ER, coef=2)
res_PPARA <- lfcShrink(dds_PPARA, coef=2)

#Writing results
write.csv(res_DNA,'4__DNA_deseq_results.csv')
write.csv(res_ER,'4__ER_deseq_results.csv')
write.csv(res_PPARA,'4__PPARA_deseq_results.csv')

#Norm counts
write.csv(counts(dds_DNA, normalized=TRUE),'DNA_deseq_norm_counts.csv')
write.csv(counts(dds_ER, normalized=TRUE),'ER_deseq_norm_counts.csv')
write.csv(counts(dds_PPARA, normalized=TRUE),'PPARA_deseq_norm_counts.csv')

#Adjust p-value of DESeq results

res_DNA_ordered <- res_DNA[order(res_DNA$pvalue),]
res_ER_ordered <- res_ER[order(res_ER$pvalue),]
res_PPARA_ordered <- res_PPARA[order(res_PPARA$pvalue),]

write.csv(res_DNA_ordered,'padj_DNA_summary.csv')
write.csv(res_ER_ordered,'padj_ER_summary.csv')
write.csv(res_PPARA_ordered,'padj_PPARA_summary.csv')

sum(res_DNA_ordered$padj < 0.05, na.rm=TRUE)
sum(res_ER_ordered$padj < 0.05, na.rm=TRUE)
sum(res_PPARA_ordered$padj < 0.05, na.rm=TRUE)

#Seperate sig DE to be placed in histogram against log2fold change

padj_DNA <- as.data.frame(res_DNA_ordered) %>% filter(padj < 0.05)
padj_ER <- as.data.frame(res_ER_ordered) %>% filter(padj < 0.05)
padj_PPARA <- as.data.frame(res_PPARA_ordered) %>% filter(padj < 0.05)

padj_DNA <- as.data.frame(padj_DNA)
padj_ER <- as.data.frame(padj_ER)
padj_PPARA <- as.data.frame(padj_PPARA)

ggplot(padj_DNA, aes(log2FoldChange)) +
  geom_histogram() +
  xlim(-5,5) +
  ggtitle("Histogram MoA: DNA")

ggplot(padj_ER, aes(log2FoldChange)) +
  geom_histogram() +
  xlim(-5,5)+
  ggtitle("Histogram MoA: ER")

ggplot(padj_PPARA, aes(log2FoldChange)) +
  geom_histogram() +
  xlim(-5,5)+
  ggtitle("Histogram MoA: PPARA")

#Scatter plots

ggplot(data=padj_DNA, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(col = "blue") +
  ggtitle("Volcano Plot MoA: DNA")

ggplot(data=padj_ER, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(col = "blue") +
  ggtitle("Volcano Plot MoA: ER")

ggplot(data=padj_PPARA, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(col = "blue") +
  ggtitle("Volcano Plot MoA: PPARA")


#Unused plot code
#EnhancedVolcano(padj_DNA,
#                lab = '',
#                x = 'log2FoldChange',
#                y = 'padj',
#                title = "Volcano Plot of MoA: DNA")


#EnhancedVolcano(padj_ER,
#                lab = rownames(padj_ER),
 #               x = 'log2FoldChange',
  #              y = 'padj',
   #             title = "Volcano Plot of MoA: ER")


#EnhancedVolcano(padj_PPARA,
 #               lab = rownames(padj_PPARA),
  #              x = 'log2FoldChange',
   #             y = 'padj',
    #            title = "Volcano Plot of MoA: PPARA")
