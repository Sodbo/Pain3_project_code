for s in {1150,4249,790,1307,3432,1076} 
do
echo $s
run_smrheidi --gwas-id-1 4532 --gwas-id-2="$s" --snp-id='rs6857'
done

run_smrheidi_report --gwas-id-1=4532 --output-dir='./' --output-file='cBP' --only-table

