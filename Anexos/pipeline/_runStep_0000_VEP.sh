#!/bin/bash
#$ -P prod
#$ -N VEP
#$ -A TFM_Alba
#$ -l thread=2
#$ -l h_vmem=16G
#$ -o _log/VEP.$TASK_ID.stdout
#$ -e _log/VEP.$TASK_ID.stderr


mkdir -p _out/
mkdir -p _log/
mkdir -p _tmp/javaio

#
# Get item
#
item=$(awk "NR==$SGE_TASK_ID" ../0000_dataset/Samples.lst)

/programs/prod/apps/ensembl-tools/scripts/variant_effect_predictor/vep \
--offline --cache -species homo_sapiens --cache_version 90 --dir ../../REFERENCES/VEP \
--fasta ../../REFERENCES/GATK_37/human_g1k_v37.fasta \
-i _in/SRR1517848.vcf.gz \
-o _out/SRR1517848.vcf \
--vcf \
--variant_class \
--sift b \
--polyphen b \
-gene_phenotype ClinVar dbGaP DDG2P GIANT GOA HGMD-PUBLIC HGMD-PUBLIC NHGRI-EBI GWAS catalog OMIM Teslovich \
--domains \
--hgvs \
--protein \
--symbol \
--ccds \
--uniprot \
--biotype \
--pubmed \
--plugin LoFtool \
--plugin CADD,../../REFERENCES/CADD/whole_genome_SNVs.tsv.gz
