#!/bin/bash
#$ -P prod
#$ -N VariantFiltration
#$ -A TFM_Alba
#$ -l thread=45
#$ -l h_vmem=124G
#$ -o _log/VariantFiltration.$TASK_ID.stdout
#$ -e _log/VariantFiltration.$TASK_ID.stderr

#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/Samples.lst)

mkdir -p _out
mkdir -p _tmp/javaio


#
#Selecciona los SNPs
/etc/alternatives/java_sdk_1.8.0/bin/java -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T SelectVariants \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-V _in/${item}.vcf  -selectType SNP \
-o _tmp/${item}.vcf.snp.vcf
#Aplica los filtros recomendados
/etc/alternatives/java_sdk_1.8.0/bin/java -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-V _in/${item}.vcf  --filterName "standard_SNP_hard_filter" \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
-o _tmp/${item}.snp.vcf

#
# Selecciona los INDEL
/etc/alternatives/java_sdk_1.8.0/bin/java -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T SelectVariants \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-V _in/${item}.vcf  -selectType INDEL \
-o _tmp/${item}.vcf.indel.vcf
#Aplica los filtros recomendados
/etc/alternatives/java_sdk_1.8.0/bin/java -jar /home/jdelabarrera/Apps/GATK/GenomeAnalysisTK.jar -T VariantFiltration \
-R ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-V _in/${item}.vcf --filterName "standard INDEL hard filter" \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
-o _tmp/${item}.indel.vcf


#Combina los VCF correspondiente a SNPs e INDELs.
bgzip _in/${item}.snp.vcf
tabix -p vcf _tmp/${item}.snp.vcf.gz 

bgzip _in/${item}.indel.vcf
tabix -p vcf _tmp/${item}.indel.vcf

bcftools concat _tmp/${item}.snp.vcf _tmp/${item}.indel.vcf -o _tmp/${item}.vcf 

#Nos quedamos solo con los PASS
bcftools view -f "PASS"  -O z -o _out/${item}.vcf.gz _tmp/${item}.vcf

#indexamos VCF
tabix -p vcf _out/${item}.vcf.gz

