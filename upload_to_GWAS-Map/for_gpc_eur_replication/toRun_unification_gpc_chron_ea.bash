#!/bin/bash
for i in {1..4}
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_GPC${i}.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/upload_gpc/mapping_gpc_chron_pain_ea.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/upload_gpc/gpc${i}_chron_repl_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/unification_gpc_results/GPC${i} \
  --output-file=gpc${i}_output
done
