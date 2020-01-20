#!/bin/bash

GPC=4

mkdir /storage/projects/PainOmics/mv_gwas_results/GPC_disc_COJO/cojo_out_100k/gpc"$GPC"/ -m 777

for i in {1..22}
do
   /storage/projects/PainOmics/BP_crude_GWAS/gcta_1.90.0beta/gcta64 \
   --bfile /storage/projects/PainOmics/MV_GWAS/100k_bed_filtered/100k_chr"$i" \
   --maf 0.0002 \
   --cojo-slct \
   --cojo-p 6e-8 \
   --chr $i \
   --cojo-wind 5000 \
   --cojo-file /storage/projects/PainOmics/mv_gwas_results/GPC_disc_COJO/gpc_sumstats/20181011_DISC_GPC${GPC}.for_cojo \
   --out /storage/projects/PainOmics/mv_gwas_results/GPC_disc_COJO/cojo_out_100k/gpc"$GPC"/GPC${GPC}_265K_chr"$i".cojo
#   --out /storage/projects/PainOmics/mv_gwas_results/GPC_disc_COJO/tmp/GPC${GPC}_265K_chr"$i".cojo
done
