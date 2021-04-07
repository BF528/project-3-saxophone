#Set working directory to project directory
setwd('/projectnb/bf528/users/saxophone/project3/analyst')

#Import tables from DESeq and Limma
#Beta-Estradiol (ER)
deseq <- read.csv('/projectnb/bf528/users/saxophone/project3/programmer/deliverables_part4/4__ER_deseq_results.csv')
#Bezafibrate (PPARA)
#deseq <- read.csv('/projectnb/bf528/users/saxophone/project3/programmer/deliverables_part4/4__PPARA_deseq_results.csv')
#N-nitrosodiethylamine (DNA)
#deseq <- read.csv('/projectnb/bf528/users/saxophone/project3/programmer/deliverables_part4/4__DNA_deseq_results.csv')
deseq <- deseq[deseq$padj<0.05,]
deseq <- deseq[abs(deseq$log2FoldChange)>1.5,]
deseq <- na.omit(deseq)

#Beta-Estradiol (ER)
limma <- read.csv('/projectnb/bf528/users/saxophone/project3/analyst/BE_limma_results.csv')
#Bezafibrate (PPARA)
#limma <- read.csv('/projectnb/bf528/users/saxophone/project3/analyst/BF_limma_results.csv')
#N-nitrosodiethylamine (DNA)
#limma <- read.csv('/projectnb/bf528/users/saxophone/project3/analyst/N_limma_results.csv')
limma <- na.omit(limma)

#Import table for ID mapping
ID.map <- read.csv('/project/bf528/project_3/refseq_affy_map.csv')

### -------- MATCH THE REFSEQ IDs AND PROBESET IDs -------- ###

#Change the ID column in DESeq, Limma from X to REFSEQ, PROBEID
names(deseq)[names(deseq) == 'X'] <- 'REFSEQ'
names(limma)[names(limma) == 'X'] <- 'PROBEID'

#Merge ID columns in DESeq and Limma with the ID mapping table
deseq.ID.merge <- merge(deseq, ID.map, c('REFSEQ'))
final.merge <- merge(limma, deseq.ID.merge, c('PROBEID'))

### -------- CALCULATE THE CONCORDANCE OF RNASEQ AND MICROARRAY DATA -------- ###

##Determine agreement in directionality of fold change
#Compare the sign of fold change between DESeq and Limma - matching sign == in agreement  
final.merge$deseqLFCsign <- sign(final.merge$log2FoldChange)
final.merge$limmaLFCsign <- sign(final.merge$logFC)
final.merge <- final.merge[final.merge$deseqLFCsign == final.merge$limmaLFCsign,]

n1 <- nrow(deseq)
n2 <- nrow(limma)
N <- 22229                #Total number of genes in Norway Rat, according to NCBI  
n0 <- nrow(final.merge)   #Observed intersection -> overlapping genes between both experiments

#Background corrected intersection
x <- (n0 * N - n1 * n2)/(n0 + N - n1 - n2)

#Concordance calculation
concordance <- (2 * x)/(n1 + n2)

#Concordance Plot vs. Number of Differentially Expressed Genes
df1 <- data.frame(DEGs = c(427,472,929), conc = c(0.106804370875599,0.249569241843679,0.0183976061659833)*100, labels = c("ER","PPARA","DNA"))
plot(x = df1$DEGs, y = df1$conc, xlab = "Number of DEGs from RNASeq", ylab = "Concordance of DEGs (%)")

df2 <- data.frame(DEGs = c(37,115,36), conc = c(0.106804370875599,0.249569241843679,0.0183976061659833)*100, labels = c("ER","PPARA","DNA"))
plot(x = df2$DEGs, y = df2$conc, xlab = "Number of DEGs from Microarray", ylab = "Concordance of DEGs (%)", xlim = c(30,120))

### -------- CALCULATE ABOVE-MEDIAN AND BELOW-MEDIAN EXPRESSION FOR GENES IN RNASEQ AND MICROARRAY DATA  -------- ###

##Median Expression Values with Above-Median and Below-Median Tables
#RNASeq data
deseq.median <- median(deseq$baseMean)
deseq.abovemed <- deseq[deseq$baseMean > deseq.median,]
deseq.belowmed <- deseq[deseq$baseMean < deseq.median,]

#Microarray data
limma.median <- median(limma$AveExpr)
limma.abovemed <- limma[limma$AveExpr > limma.median,]
limma.belowmed <- limma[limma$AveExpr < limma.median,]

##Merge Above-Median and Below-Median Tables
#Above-Median
abovemed.merge <- merge(deseq.abovemed, ID.map, c('REFSEQ'))
final.abovemed.merge <- merge(limma.abovemed, abovemed.merge, c('PROBEID'))

#Below-Median
belowmed.merge <- merge(deseq.belowmed, ID.map, c('REFSEQ'))
final.belowmed.merge <- merge(limma.belowmed, belowmed.merge, c('PROBEID'))

### -------- CALCULATE THE CONCORDANCE OF ABOVE-MEDIAN AND BELOW-MEDIAN EXPRESSION -------- ###

##Agreement in Directionality of Fold Change for Median Expression Analysis
#Above-Median
final.abovemed.merge$deseqLFCsign <- sign(final.abovemed.merge$log2FoldChange)
final.abovemed.merge$limmaLFCsign <- sign(final.abovemed.merge$logFC)
final.abovemed <- final.abovemed.merge[final.abovemed.merge$deseqLFCsign == final.abovemed.merge$limmaLFCsign,]

#Below-Median
final.belowmed.merge$deseqLFCsign <- sign(final.belowmed.merge$log2FoldChange)
final.belowmed.merge$limmaLFCsign <- sign(final.belowmed.merge$logFC)
final.belowmed <- final.belowmed.merge[final.belowmed.merge$deseqLFCsign == final.belowmed.merge$limmaLFCsign,]

##Calculate background-corrected intersection and concordance
#Above-Median
n1.a <- nrow(deseq.abovemed)
n2.a <- nrow(limma.abovemed)
n0.a <- nrow(final.abovemed.merge)
N.a <- 22229

x.a <- (n0.a * N.a - n1.a * n2.a)/(n0.a + N.a - n1.a - n2.a)
conc.a <- (2 * x.a)/(n1.a + n2.a)

#Below-Median
n1.b <- nrow(deseq.belowmed)
n2.b <- nrow(limma.belowmed)
n0.b <- nrow(final.belowmed.merge)
N.b <- 22229 

x.b <- (n0.b * N.b - n1.b * n2.b)/(n0.b + N.b - n1.b - n2.b)
conc.b <- (2 * x.b)/(n1.b + n2.b)

##Bar plot of overall, above-median, and below-median concordance
par(mar=c(5,5,2,2))
barplot(c(0.106804370875599,0.0859412004450172,0.00723960494859172, 0.249569241843679, 0.320232248965575, 0.0234782419319211, 0.0183976061659833, 0.0196086323535471, 0.0068884742942228)*100, ylab="Concordance (%)", ylim = c(0,40))
box(which="plot",lty="solid")
