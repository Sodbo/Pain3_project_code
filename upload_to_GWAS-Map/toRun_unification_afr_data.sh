#!/bin/bash
for pain in 'Back' 'Neck' 'Face' 'Head' 'Knee' 'Stomach'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/filtered_scaled/MV_${pain}_Repl_AFR_gwas.bgen.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/for_mapping/mapping_MV_pain_replication.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/descriptors/${pain}_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/unification_results/ \
  --output-file=${pain}_output
done
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/filtered_scaled/MV_Hip_Repl_AFR_gwas.bgen.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/for_mapping/mapping_Hip_pain_replication.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/descriptors/Hip_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/unification_results/ \
  --output-file=Hip_output

