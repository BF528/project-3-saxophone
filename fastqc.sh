#!/bin/bash -l

#$ -N fastqc
#$ -j y
#$ -pe omp 6

module load fastqc

fastqc -t 6 SRR1177984_1.fastq.gz SRR1177984_2.fastq.gz
fastqc -t 6 SRR1177985_1.fastq.gz SRR1177985_2.fastq.gz
fastqc -t 6 SRR1177986_1.fastq.gz SRR1177986_2.fastq.gz
fastqc -t 6 SRR1177960_1.fastq.gz SRR1177960_2.fastq.gz
fastqc -t 6 SRR1178023_1.fastq.gz SRR1178023_2.fastq.gz
fastqc -t 6 SRR1178049_1.fastq.gz SRR1178049_2.fastq.gz
fastqc -t 6 SRR1177967_1.fastq.gz SRR1177967_2.fastq.gz
fastqc -t 6 SRR1177968_1.fastq.gz SRR1177968_2.fastq.gz
fastqc -t 6 SRR1177971_1.fastq.gz SRR1177971_2.fastq.gz
