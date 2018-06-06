#!/bin/bash
#$ -P prod
#$ -N BQSR
#$ -A TFM_Alba
#$ -l thread=41
#$ -l h_vmem=123G
#$ -o _log/BQSR..$TASK_ID.stdout
#$ -e _log/BQSR..$TASK_ID.stderr

#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/Samples.lst)

mkdir -p _out
mkdir -p _tmp/javaio


/usr/lib/jvm/java-1.8.0/jre/bin/java -Xmx100G -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 40 \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-I _in/${item}.bam \
-o _tmp/${item}.recal_data.table \
-knownSites ../../REFERENCES/GATK_37/1000G_phase1.snps.high_confidence.b37.vcf.gz \
-knownSites ../../REFERENCES/GATK_37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
-knownSites ../../REFERENCES/GATK_37/dbsnp_137.b37.vcf.gz

/usr/lib/jvm/java-1.8.0/jre/bin/java -Xmx100G -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 40 \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-I _in/${item}.bam -BQSR _tmp/${item}.recal_data.table \
-o _tmp/${item}.post_recal_data.table \
-knownSites ../../REFERENCES/GATK_37/1000G_phase1.snps.high_confidence.b37.vcf.gz \
-knownSites ../../REFERENCES/GATK_37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
-knownSites ../../REFERENCES/GATK_37/dbsnp_137.b37.vcf.gz

/usr/bin/time --verbose -a -o _log/runtime.${item}.log \
/usr/lib/jvm/java-1.8.0/jre/bin/java -Xmx100G -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T PrintReads -nct 40 \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-I _in/${item}.bam -BQSR _tmp/${item}.post_recal_data.table \
-o _out/${item}.bam

