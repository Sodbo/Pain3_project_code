#!/bin/bash
# This script is to run COJO on GPC1 MA Eur data with p = 1e-05 threshold

path='/storage/projects/PainOmics/mv_gwas_results/GPC_MA_eur_COJO/data'

for i in {1..22}
do
   /storage/projects/PainOmics/BP_crude_GWAS/gcta_1.90.0beta/gcta64 \
   --bfile /storage/projects/PainOmics/MV_GWAS/100k_bed_filtered/100k_chr"$i" \
   --maf 0.05 \
   --cojo-slct \
   --cojo-p 1e-5 \
   --chr $i \
   --cojo-wind 5000 \
   --cojo-file ${path}/input_data/gpc1_ma_eur.for_cojo \
   --out ${path}/output_data/gpc1_ma_eur_chr"$i"_1e-5.cojo
done
