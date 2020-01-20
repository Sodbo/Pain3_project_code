#!/bin/bash
for pain in 'Back' 'Neck' 'Face' 'Head' 'Hip' 'Knee' 'Stomach'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/discovery/filtered_scaled/MV_${pain}_Discovery_gwas.bgen.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/discovery/mapping_MV_pain.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/discovery/descriptors/${pain}_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/discovery/unification_results/ \
  --output-file=${pain}_output
done

