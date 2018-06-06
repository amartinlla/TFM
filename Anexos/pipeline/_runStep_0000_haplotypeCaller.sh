#!/bin/bash
#$ -P prod
#$ -N haplotypeCaller
#$ -A TFM_Alba
#$ -l thread=45
#$ -l h_vmem=124G
#$ -o _log/haplotypeCaller.$TASK_ID.stdout
#$ -e _log/haplotypeCaller.$TASK_ID.stderr

#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/Samples.lst)

mkdir -p _out
mkdir -p _tmp/javaio

/etc/alternatives/java_sdk_1.8.0/bin/java -jar -Xmx100G /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 44 \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
--dbsnp ../../REFERENCES/GATK_37/dbsnp_137.b37.vcf.gz \
-L ../../REFERENCES/Regions.interval_list \
-I _in/${item}.bam -o _out/${item}.vcf \
-G Standard -G AS_Standard
