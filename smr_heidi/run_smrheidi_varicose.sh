for s in {4053,1544,2065,1282656,3432} 
do
echo $s
run_smrheidi --gwas-id-1 4041 --gwas-id-2="$s" --snp-id='rs505922' --overwrite --version 11
done

run_smrheidi_report --gwas-id-1=4041 --output-dir='./' --output-file='srVar' --only-table --version 11




