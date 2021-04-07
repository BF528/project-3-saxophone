#4.1 
# subset the data to get the significant DE genes (adjusted p-value <0.05) 
ER_DEG<- read.csv("4__ER_deseq_results.csv", header=TRUE)
PPARA_DEG<-read.csv("4__PPARA_deseq_results.csv", header=TRUE)
DNA_DEG<-read.csv("4__DNA_deseq_results.csv", header=TRUE)
sig_ER <- subset(ER_DEG, pvalue < 0.05)
sig_PPARA <- subset(PPARA_DEG, pvalue < 0.05)
sig_DNA <- subset(DNA_DEG, pvalue < 0.05)
head(sig_ER)

#cut off at log2foldchange > 1.5
sig_ER1 <- subset(sig_ER, abs(log2FoldChange) > 1.5)
sig_PPARA1 <- subset(sig_PPARA, abs(log2FoldChange) > 1.5)
sig_DNA1 <- subset(sig_DNA, abs(log2FoldChange) > 1.5)
dim(sig_ER1)
dim(sig_PPARA1)
dim(sig_DNA1)

#write to csv files
sig_ER1<- cbind(sig_ER1)
sig_PPARA1 <- cbind(sig_PPARA1)
sig_DNA1<- cbind(sig_DNA1)
write.csv(sig_ER1, "sig_ER1.csv")
write.csv(sig_PPARA1, "sig_PPARA1.csv")
write.csv(sig_DNA1, "sig_DNA1.csv")


#4.2
ER_norm_counts<-read.csv("ER_deseq_norm_counts.csv", header=TRUE)
PPARA_norm_counts<-read.csv("PPARA_deseq_norm_counts.csv", header=TRUE)
DNA_norm_counts<-read.csv("DNA_deseq_norm_counts.csv", header=TRUE)

combined<- as.data.frame(c(ER_norm_counts, PPARA_norm_counts, DNA_norm_counts))

#renaming first column to subset them together 
colnames(ER_norm_counts)[1] <- "ID"
colnames(PPARA_norm_counts)[1] <- "ID"
colnames(DNA_norm_counts)[1] <- "ID"

#merged all three normalized counts together 
merged_MOA<- merge(ER_norm_counts[,1:10],PPARA_norm_counts[,1:10], by= "ID") 
merged_MOA<- merge(merged_MOA, DNA_norm_counts[, 1:10], by="ID")

#gets rid of first column
merged_MOA<- merged_MOA[, -1, drop=FALSE]
merged_MOA_new <- as.matrix(merged_MOA)

#make final heatmap 
labels<- c("Control_BEZ_60", "Control_BEZ_23", "Control_BEZ_49", "Control_BEZ_67", "PPARA_68", "PPARA_71", "PPARA_04", "PPARA_06",  "DNA_Damage_13", "DNA_Damage_64", "DNA_Damage_74", "DNA_Damage_75" )
heatmap(merged_MOA_new,labCol=labels , ylab = "Genes")