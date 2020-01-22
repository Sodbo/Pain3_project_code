#!/bin/bash
for i in {1..4}
do
run_upload --gwas-path=/home/ubuntu/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC${i}/gpc${i}_output_done.csv
done
