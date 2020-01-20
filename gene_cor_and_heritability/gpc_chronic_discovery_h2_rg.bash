#!/bin/bash

for NPC in {1..4}

do

/home/common/projects/ldsc_tool/munge_sumstats.py \
--sumstat /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC".txt \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC" \
--merge-alleles /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--N-col N \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p p \
--signed-sumstats b,0 \
--frq eaf \
--ignore Z,pos,se,chr 

/home/common/projects/ldsc_tool/ldsc.py \
--h2 /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC".sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC"_h2

/home/common/projects/ldsc_tool/ldsc.py \
--rg /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC".sumstats.gz,/mnt/polyomica/projects/mv_gwas/MA/all_eur/06_GPC/20181004_DISC_GPC"$NPC".sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC"_GPC"$NPC"_rg

done
