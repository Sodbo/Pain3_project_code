# This code is for reformating of GWAS summary statistics of GPC1 MA Eur
# for further COJO analythis

library(data.table)

gpc_in <- '/storage/projects/PainOmics/mv_gwas_results/GPC_MA_eur_COJO/data/input_data/gpc1_Pval'

tmp <- fread(gpc_in)

tmp <- tmp[,c('MarkerName', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P', 'n_total')]

colnames(tmp) <- c('SNP','A1','A2','freq','b','se','p','N')

tmp$A1 <- toupper(tmp$A1)

tmp$A2 <- toupper(tmp$A2)

gpc_out <- '/storage/projects/PainOmics/mv_gwas_results/GPC_MA_eur_COJO/data/input_data/gpc1_ma_eur.for_cojo'

fwrite(tmp, file = gpc_out, sep = '\t', dec = '.')

