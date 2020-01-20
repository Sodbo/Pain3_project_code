#!/bin/bash

for GWAS_ID in 4461 4462 4466 4465 4469 4470 4473 4474;
do
python3 ~/code_folder/gwas_v2/gwas_v2/analysis/ldscore/run_ldscore.py --h2 --gwas-id ${GWAS_ID};
done
