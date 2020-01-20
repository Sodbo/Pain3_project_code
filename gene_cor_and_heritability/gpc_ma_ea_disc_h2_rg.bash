#!/bin/bash

for NPC in {1..4}

do

/home/common/projects/ldsc_tool/munge_sumstats.py \
--sumstat /mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/GPC"$NPC"/gpc"$NPC"_cau1.tbl \
--out /mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/gene_cor/GPC"$NPC"/20181012_GPC"$NPC"_MA \ 
--merge-alleles /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--N-col n_total \
--snp MarkerName \
--a1 Allele1 \
--a2 Allele2 \
--p P-value \
--signed-sumstats Effect,0 \
--frq Freq1 \
--ignore FreqSE,StdErr,Direction,HetISq,HetChiSq,HetDf,HetPVal 

/home/common/projects/ldsc_tool/ldsc.py \
--h2 /mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/gene_cor/GPC"$NPC"/20181012_GPC"$NPC"_MA.sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/gene_cor/GPC"$NPC"/20181012_GPC"$NPC"_MA_h2

/home/common/projects/ldsc_tool/ldsc.py \
--rg /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC".sumstats.gz,/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/gene_cor/GPC"$NPC"/20181012_GPC"$NPC"_MA.sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/gene_cor/GPC"$NPC"/20181012_GPC"$NPC"_MA_GPC"$NPC"_DISC_rg

done
