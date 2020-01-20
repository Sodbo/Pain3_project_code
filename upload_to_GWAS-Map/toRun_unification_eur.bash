#!/bin/bash
for pain in 'Back' 'Neck' 'Hip' 'Knee' 'Head' 'Face' 'Stom'
do
python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
  --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_${pain}_Repl.EA_gwas.BGEN.stats.txt \
  --mapping-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/upload_to_db/mapping_chronic_pain_replication.json \
  --descriptors-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/upload_to_db/chronic_${pain}_eur_descriptor.json \
  --output-dir=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/unification_results/ \
  --output-file=${pain}_output
done

