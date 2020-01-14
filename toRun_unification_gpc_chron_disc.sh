#!/bin/bash
for i in {1..4}
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC${i}.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/upload_to_db/mapping_gpc_chron_pain_discovery.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/upload_to_db/gpc${i}_chron_disc_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC${i} \
  --output-file=gpc${i}_output
done
