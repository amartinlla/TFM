#!/bin/bash
#$ -P prod
#$ -N trimmomatic
#$ -A TFM_Alba
#$ -l thread=44
#$ -l h_vmem=124G
#$ -o _log/trimmomatic.$TASK_ID.stdout
#$ -e _log/trimmomatic.$TASK_ID.stderr

#
#
#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/FASTQFilename.lst)

#
# Create dirs
mkdir -p _out
mkdir -p _tmp/javaio
mkdir -p _log

/usr/lib/jvm/java-1.8.0/jre/bin/java -Xmx90G -jar /home/jdelabarrera/Apps/Trimmomatic/trimmomatic.jar PE -threads 44  \
_in/${item}_1.fastq.gz _in/${item}_2.fastq.gz \
_out/${item}.R1.fastq.gz _out/${item}_unpaired.fastq.gz _out/${item}.R2.fastq.gz _out/${item}_unpaired.fastq.gz \
ILLUMINACLIP:/home/jdelabarrera/Apps/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
