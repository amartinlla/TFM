#!/bin/bash
#$ -P prod
#$ -N MarkDuplicates
#$ -A TFM_Alba
#$ -l thread=3
#$ -l h_vmem=32G
#$ -o _log/MarkDuplicates.$TASK_ID.stdout
#$ -e _log/MarkDuplicates.$TASK_ID.stderr

#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/Samples.lst)

mkdir -p _out
mkdir -p _tmp/javaio

#Markduplicates
/usr/lib/jvm/java-1.8.0/jre/bin/java -jar -Xmx16G /programs/prod/apps/picard-tools/MarkDuplicates.jar \
I=_in/${item}.bam \
O=_out/${item}.bam  \
VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=true \
M=_log/${item}.output.metrics
