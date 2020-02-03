# run on mga-n4
# srun -w mga-n4 --pty bash -i
# run using R (Rscript name-of-script.R)

setwd("/mnt/polyomica/projects/mv_gwas/results/20200120/03_data_twostep_app/")
path_input="chr/chr"
path_output="./"

library(data.table)
chr=1
x=fread(paste0(path_input,chr,".raw"),data.table=F)
y=x

for(chr in c(2:20,22)){
		x=fread(paste0(path_input,chr,".raw"),data.table=F)
	y=cbind(y,x)
}

full=y[,grepl("rs",colnames(y))]
full=cbind(y[,c(2,5)],full)
fwrite(x=full,file=paste0(path_output,"01_joint_all_snps.txt"),sep="\t",col.names=T,row.names=F,quote=F)


