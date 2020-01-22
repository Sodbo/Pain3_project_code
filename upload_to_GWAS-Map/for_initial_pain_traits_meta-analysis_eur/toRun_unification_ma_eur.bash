#!/bin/bash
for pain in 'Back' 'Neck' 'Hip' 'Knee' 'Head' 'Face' 'Stom'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/01_METAL_out/${pain}_chronic_EUR_MA_440K1_Pval.tbl \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/upload_to_db/mapping_ma_chron_pain_eur.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/upload_to_db/chronic_${pain}_ma_eur_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/unification_results/${pain}/ \
  --output-file=${pain}_output
done
