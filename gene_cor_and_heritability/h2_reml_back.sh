#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -N h2.back
#$ -j y
#$ -l h_vmem=10G
#$ -pe smp 10

## SNP-based h2 for chronic pain
## for MV-pain project

DIRbolt=/mnt/lustre/users/DTR/k1471250/tools/BOLT-LMM_v2.3.2/
DIRbgen=/mnt/lustre/datasets/ukbiobank/June2017/Imputed/

$DIRbolt/bolt \
--bed=/mnt/lustre/datasets/ukbiobank/June2017/Genotypes/ukb_binary_v2.bed \
--bim=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/ukb18219.bim \
--fam=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/ukb18219.fam \
--remove=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/phen/chron/mv_mar18_cron.disc.EA_excl.txt \
--phenoFile=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/phen/chron/mv_mar18_chron.phen.txt \
--phenoCol=back \
--covarFile=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/phen/chron/mv_mar18_chron.phen.txt \
--covarCol=batch \
--covarCol=Sex \
--qCovarCol=Age \
--qCovarCol=pc{1:10} \
--LDscoresFile=$DIRbolt/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=$DIRbolt/tables/genetic_map_hg19_withX.txt.gz \
--reml \
--noMapCheck \
--numThreads=10 \
--statsFile=h2_reml_back.gz \
--verboseStats 

## 2>&1 | tee h2_reml_mv.log

##--bgenMinMAF=0.01 \
##--bgenMinINFO=0.8 \
##--maxMissingPerSnp=0.05 \
##--maxMissingPerIndiv=0.05 \
##--bgenFile=$DIRbgen/ukb_imp_chr{1:22}_v3.bgen \
##--sampleFile=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/v3/ukb18219_imp_chr1_v3_s487395.sample \
##--statsFileBgenSnps=h2_reml_back.gz \

