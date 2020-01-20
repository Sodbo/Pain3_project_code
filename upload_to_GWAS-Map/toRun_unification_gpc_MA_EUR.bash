#!/bin/bash
for i in {1..4}
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/GPC${i}/gpc${i}_Pval \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/upload_to_db/mapping_gpc_ma_chron_pain_eur.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/upload_to_db/gpc${i}_chron_ma_eur_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC${i} \
  --output-file=gpc${i}_output
done
