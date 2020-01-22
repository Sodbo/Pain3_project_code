#!/bin/bash
for pain in 'Back' 'Neck' 'Knee' 'Head' 'Face' 'Stom'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/scaled_filtered/MV_${pain}_Disc_gwas.BGEN.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/upload_to_db/mapping_chronic_pain_discovery.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/upload_to_db/chronic_${pain}_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/${pain}/ \
  --output-file=${pain}_output
done
