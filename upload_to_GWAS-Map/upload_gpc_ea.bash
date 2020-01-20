#!/bin/bash
for i in {1..4}
do
run_upload --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/unification_gpc_results/GPC${i}/gpc${i}_output_done.csv
done
