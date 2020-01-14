#!/bin/bash
source ~/.bashrc
for pain in 'Back' 'Neck' 'Face' 'Head' 'Hip' 'Knee' 'Stom'
do
python3 ~/code_folder/gwas_v2/gwas_v2/unifier/run_upload.py --gwas-path=~/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/${pain}/${pain}_output_done.csv
done
