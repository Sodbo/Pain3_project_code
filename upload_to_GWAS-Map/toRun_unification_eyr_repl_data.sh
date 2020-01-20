#!/bin/bash
for pain in 'Back' 'Neck' 'Face' 'Head' 'Hip' 'Knee' 'Stom'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/replication/filtered_scaled/MV_${pain}_Repl_CAU_gwas.bgen.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/CAU_repl_07092018/mapping_MV_pain.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/CAU_repl_07092018/descriptors/${pain}_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/CAU_repl_07092018/unification_results_prod/$pain/ \
  --output-file=${pain}_output
done

