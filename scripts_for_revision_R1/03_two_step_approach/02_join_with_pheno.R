# run on mga-n4
# srun -w mga-n4 --pty bash -i
# run using R (Rscript name-of-script.R)

setwd("/mnt/polyomica/projects/mv_gwas/results/20200120/03_data_twostep_app/")

x=fread("mv_phen_GIPs.txt",data.table=F)
x=na.omit(x)
y=fread("01_joint_all_snps.txt",data.table=F)

ind=which(x[,1]%in%y[,1])
length(ind)
x=x[ind,]

ind=which(y[,1]%in%x[,1])
length(ind)
y=y[ind,]

ind=match(x[,1],y[,1])
table(y[ind,1]==x[,1])
y=y[ind,]
table(y[,1]==x[,1])

g=y
save(file="02_pheno_and_geno.RData",list=c("g","x"))

