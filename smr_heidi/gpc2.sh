export PROD=TRUE

g2s=$(cat 20190217_gwasid_for_smrheidi.txt)

for g2 in $g2s
do
echo $g2

run_smrheidi \
--gwas-id-1 4530 \
--gwas-id-2 $g2 \
--snp-list snp_list_gpc2.txt \
--version 3 \
--overwrite

done

