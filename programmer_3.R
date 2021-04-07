setwd("C:/Users/jason/Desktop/BU/BF528/project_3")
library(ggplot2)
library(dplyr)
library(tidyverse)

#Reading in feature counts, subsetting data and writing files for each sample
n_a<-read.table(file = "SRR1177960", header = T)
names(n_a)[7]<- "SRR1177960"
n_a = subset(n_a, select = -c(2, 3, 4, 5, 6))
write.csv(n_a, "SRR1177960_counts.csv", row.names = F)

n_b<-read.table(file = "SRR1177967", header = T)
names(n_b)[7]<- "SRR1177967"
n_b = subset(n_b, select = -c(2, 3, 4, 5, 6))
write.csv(n_b, "SRR1177967_counts.csv", row.names = F)

n_c<-read.table(file = "SRR1177968", header = T)
names(n_c)[7]<- "SRR1177968"
n_c = subset(n_c, select = -c(2, 3, 4, 5, 6))
write.csv(n_c, "SRR1177968_counts.csv", row.names = F)

n_d<-read.table(file = "SRR1177971", header = T)
names(n_d)[7]<- "SRR1177971"
n_d = subset(n_d, select = -c(2, 3, 4, 5, 6))
write.csv(n_d, "SRR1177971_counts.csv", row.names = F)

n_e<-read.table(file = "SRR1177984", header = T)
names(n_e)[7]<- "SRR1177984"
n_e = subset(n_e, select = -c(2, 3, 4, 5, 6))
write.csv(n_e, "SRR1177984_counts.csv", row.names = F)

n_f<-read.table(file = "SRR1177985", header = T)
names(n_f)[7]<- "SRR1177985"
n_f = subset(n_f, select = -c(2, 3, 4, 5, 6))
write.csv(n_f, "SRR1177985_counts.csv", row.names = F)

n_g<-read.table(file = "SRR1177986", header = T)
names(n_g)[7]<- "SRR1177986"
n_g = subset(n_g, select = -c(2, 3, 4, 5, 6))
write.csv(n_g, "SRR1177986_counts.csv", row.names = F)

n_h<-read.table(file = "SRR1178049", header = T)
names(n_h)[7]<- "SRR1178049"
n_h = subset(n_h, select = -c(2, 3, 4, 5, 6))
write.csv(n_h, "SRR1178049_counts.csv", row.names = F)

n_i<-read.table(file = "SRR1178023", header = T)
names(n_i)[7]<- "SRR1178023"
n_i = subset(n_i, select = -c(2, 3, 4, 5, 6))
write.csv(n_i, "SRR1178023_counts.csv", row.names = F)


#Concatenating files
total1 = cbind(n_e, n_f, n_g, n_a, n_i, n_h, n_b, n_c, n_d)
total2 = subset(total1, select = -c(3,5,7,9,11,13,15,17))
 
write.csv(total2, "conc_counts.csv", row.names = F, row.names = 1)

#Boxplot

#boxplt <- total2[,-1]
#rownames(boxplt) <- total2[,1]

#ggplot(data = boxplt, aes(x = samples, y = counts)) +
#  geom_boxplot()



