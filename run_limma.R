#Set working directory to project directory
setwd('/projectnb/bf528/users/saxophone/project3/analyst')

#Import library
library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_4_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# subset the samples table for each sample
samples.BE <- samples[(samples$chemical == 'BETA-ESTRADIOL') | (samples$chemical == 'Control'),]
samples.BF <- samples[(samples$chemical == 'BEZAFIBRATE') | (samples$chemical == 'Control'),]
samples.N <- samples[(samples$chemical == 'N-NITROSODIETHYLAMINE') | (samples$chemical == 'Control'),]

# subset the full expression matrix to just those in this comparison
rma.subset.BE <- rma[,paste0('X',samples.BE$array_id)]
rma.subset.BF <- rma[,paste0('X',samples.BF$array_id)]
rma.subset.N <- rma[,paste0('X',samples.N$array_id)]

# design matrix for Beta-Estradiol
design.BE <- model.matrix(
  ~factor(
    samples.BE$chemical,
    levels=c('Control','BETA-ESTRADIOL')
  )
)

# design matrix for Bezafibrate
design.BF <- model.matrix(
  ~factor(
    samples.BF$chemical,
    levels=c('Control','BEZAFIBRATE')
  )
)

# design matrix for N-Nitrosodiethylamine
design.N <- model.matrix(
  ~factor(
    samples.N$chemical,
    levels=c('Control','N-NITROSODIETHYLAMINE')
  )
)

# column names
colnames(design.BE) <- c('Intercept','BETA-ESTRADIOL')
colnames(design.BF) <- c('Intercept','BEZAFIBRATE')
colnames(design.N) <- c('Intercept','N-NITROSODIETHYLAMINE')

# run limma
# limma -> BE
fit.BE <- lmFit(rma.subset.BE, design.BE)
fit.BE <- eBayes(fit.BE)
t.BE <- topTable(fit.BE, coef=2, n=nrow(rma.subset.BE), adjust='BH')
t.de.BE <- t.BE[t.BE$adj.P.Val<0.05,]
t.de.BE <- t.de.BE[abs(t.de.BE$logFC)>1.5,]

# limma -> BF
fit.BF <- lmFit(rma.subset.BF, design.BF)
fit.BF <- eBayes(fit.BF)
t.BF <- topTable(fit.BF, coef=2, n=nrow(rma.subset.BF), adjust='BH')
t.de.BF <- t.BF[t.BF$adj.P.Val<0.05,]
t.de.BF <- t.de.BF[abs(t.de.BF$logFC)>1.5,]

# limma -> N
fit.N <- lmFit(rma.subset.N, design.N)
fit.N <- eBayes(fit.N)
t.N <- topTable(fit.N, coef=2, n=nrow(rma.subset.N), adjust='BH')
t.de.N <- t.N[t.N$adj.P.Val<0.05,]
t.de.N <- t.de.N[t.de.N$logFC>1.5,]

# write out the results to file
write.csv(t.de.BE,'BE_limma_results.csv')
write.csv(t.de.BF,'BF_limma_results.csv')
write.csv(t.de.N,'N_limma_results.csv')

##Histogram of LFC vs. Number of Genes
hist(t.de.BE$logFC, freq = TRUE, breaks = seq(from = -3, to = 4, by = 0.5),xlab = 'Log Fold Change', ylab = 'Frequency', ylim = c(0,25), main = NULL)
hist(t.de.BF$logFC, freq = TRUE, breaks = seq(from = -3, to = 9, by = 0.5),xlab = 'Log Fold Change', ylab = 'Frequency', ylim = c(0,70), main = NULL)
hist(t.de.N$logFC, freq = TRUE, breaks = seq(from = 0, to = 6, by = 0.5),xlab = 'Log Fold Change', xlim = c(0,6),ylab = 'Frequency', ylim = c(0,25),main = NULL)
box(which="plot",lty="solid")

##Scatterplot of LFC vs. Nominal P-value
par(mar=c(6,5,4,2)) #Adjust dimensions of the plot
options(scipen=999) #Standard Notation

plot(x = t.de.BE$logFC, y = -log10(t.de.BE$P.Value), xlab = "Log Fold Change", ylab = "-Log10(P-value)")
plot(x = t.de.BF$logFC, y = -log10(t.de.BF$P.Value), xlab = "Log Fold Change", ylab = "-Log10(P-value)")
plot(x = t.de.N$logFC, y = -log10(t.de.N$P.Value), xlab = "Log Fold Change", ylab = "-Log10(P-value)")
