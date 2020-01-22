#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -N Knee_MV_disc_bolt
#$ -j y
#$ -l h_vmem=19G
#$ -pe smp 10

## Discovery for Knee pain
## Chronic pain
## White British only, N = 265000

DIRbolt=/mnt/lustre/users/DTR/k1471250/tools/BOLT-LMM_v2.3.2/
DIRphen=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/mar18/phen/chron
DIRbgen=/mnt/lustre/datasets/ukbiobank/June2017/Imputed/

$DIRbolt/bolt \
--bed=/mnt/lustre/datasets/ukbiobank/June2017/Genotypes/ukb_binary_v2.bed \
--bim=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/mar18/ukb18219.bim \
--fam=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/mar18/ukb18219.fam \
--remove=$DIRphen/mv_mar18_cron.disc.EA_excl.txt \
--remove=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/mar18/phen/bolt.in_plink_but_not_imputed.FID_IID.968.txt \
--phenoFile=$DIRphen/mv_mar18_chron.phen.txt \
--phenoCol=knee \
--covarFile=$DIRphen/mv_mar18_chron.phen.txt \
--covarCol=Sex \
--covarCol=batch \
--qCovarCol=Age \
--qCovarCol=pc{1:10} \
--LDscoresFile=$DIRbolt/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=$DIRbolt/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--bgenMinMAF=0.0002 \
--bgenMinINFO=0.7 \
--maxMissingPerSnp=0.02 \
--maxMissingPerIndiv=0.02 \
--numThreads=10 \
--statsFile=mv_knee_disc_gwas.stats.txt.gz \
--bgenFile=$DIRbgen/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=/mnt/lustre/users/DTR/k1471250/UKBB/ukbb500/mar18/ukb18219_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=MV_Knee_Disc_gwas.BGEN.stats.txt.gz \
--noMapCheck \
--verboseStats
