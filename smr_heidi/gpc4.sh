export PROD=TRUE

g2s=$(cat ~/polyomica/projects/mv_gwas/results/20190212/20190213_gwasids_for_SMR_HEIDI.txt)

for g2 in $g2s
do
echo $g2

run_smrheidi \
--gwas-id-1 4533 \
--gwas-id-2 $g2 \
--snp-list snp_list_gpc4.txt \
--version 1 \
--overwrite

done

