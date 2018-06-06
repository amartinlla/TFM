#!/bin/bash
#$ -P prod
#$ -N mappingBwa
#$ -A TFM_Alba
#$ -l thread=45
#$ -l h_vmem=123G
#$ -o _log/mappingBwa.$TASK_ID.stdout
#$ -e _log/mappingBwa.$TASK_ID.stderr

#
#
#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/FASTQFilename.lst)

#
# Create dirs
mkdir -p _out/
mkdir -p _tmp/javaio
mkdir -p _log

#
# Generate ReadGroups
rgID='NA12878'
rgPU='Unknown'
rgSM='NA12878'
rgDT=$fecha

rgLB='SRR1517848'
rgPL='ILLUMINA'
rgPM='Unknown'
rgDS='Alba Martin'
rgCN='CNIC'
rgPI='Alba Martin'

rgString=$(printf ""@RG\\\tID:%s\\\tLB:%s\\\tPL:%s\\\tPM:%s\\\tPU:%s\\\tSM:%s\\\tCN:%s\\\tDT:%s\\\tDS:%s"" "$rgID" "$rgLB" "$rgPL" "$rgPM" "$rgPU" "$rgSM" "$rgCN" "$rgDT" "$rgDS")

bwa mem -t 45 -M -R "$rgString" ../../REFERENCES/bwa/human_g1k_v37.fasta _in/${item}.R1.fastq.gz _in/${item}.R2.fastq.gz > _tmp/${item}.sam

/usr/lib/jvm/java-1.8.0/jre/bin/java -jar -Xmx16G /programs/prod/apps/picard-tools/MergeSamFiles.jar \
INPUT=_tmp/${item}.sam \
OUTPUT=_out/${item}.bam  \
SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

samtools index _out/${item}.bam
