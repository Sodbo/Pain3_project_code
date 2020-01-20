export PROD=TRUE

g2s=$(cat ./20190212_gwasids_for_ldsc.txt)

for g2 in $g2s
do
echo $g2

run_ldscore \
--rg \
--gwas-id-1 4529 \
--gwas-id-2 $g2

done
