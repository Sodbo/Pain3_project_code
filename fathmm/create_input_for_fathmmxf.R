# Collect SNPs with P inside unit order of magnitude

library(data.table)

setwd('/mnt/polyomica/projects/mv_gwas/funk_an/fathmm/input/raw')

tab_1 <- NULL

for (f in list.files()){
	tmp <- fread(f)
	tab_1 <- rbind(tab_1, tmp)
}
