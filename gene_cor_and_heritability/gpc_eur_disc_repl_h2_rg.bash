#!/bin/bash

for NPC in {2..4}

do

/home/common/projects/ldsc_tool/munge_sumstats.py \
--sumstat /mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_GPC"$NPC".txt \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_GPC"$NPC" \
--merge-alleles /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--N-col N \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p p \
--signed-sumstats b,0 \
--frq eaf \
--ignore Z,se,chr,pos 

/home/common/projects/ldsc_tool/ldsc.py \
--h2 /mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_GPC"$NPC".sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_REPL_GPC"$NPC"_h2

/home/common/projects/ldsc_tool/ldsc.py \
--rg /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC"$NPC".sumstats.gz,/mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_GPC"$NPC".sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_DISC_GPC"$NPC"_REPL_GPC"$NPC"_rg

done
