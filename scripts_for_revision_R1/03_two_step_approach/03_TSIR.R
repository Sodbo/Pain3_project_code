#,run,on,mga-n4
#,srun,-w,mga-n4,--pty,bash,-i
#,run,using,R,(Rscript,name-of-script.R)

setwd("/mnt/polyomica/projects/mv_gwas/results/20200120/03_data_twostep_app/")

load("02_pheno_and_geno.RData")

N=nrow(x)
ind_rs=3:64
N_snps=length(ind_rs)

covs=c("Sex","pc1","pc2","pc3","pc4",
       "pc5","pc6","pc7","pc8","pc9","pc10","batch",
       "bmi","Age")

pip=0.8
alfa1=1e-5/N_snps
alfa2=0.05/N_snps
Nloops=100
TSIR_thr=0.2

Ndisc=round(pip*N)
Nrepl=N-Ndisc

#GIP1
i=1
out_d=array(NA,c(Nloops,N_snps))
out_r=array(NA,c(Nloops,N_snps))
for(i in 1:Nloops){
	print(i)
	ind_d=sample(1:N,Ndisc)
	ind_r=(1:N)[-ind_d]
	G_d=g[ind_d,ind_rs]
	G_r=g[ind_r,ind_rs]

	j=1
	for (j in 1:ncol(G_d)){
		y=x[ind_d,"GIP1"]
		G=cbind(SNP=G_d[,j],x[ind_d,covs])
		l=summary(lm(y~as.matrix(G)))
		pvald=l$coeff[2,4]
		betad=l$coeff[2,1]

		y=x[ind_r,"GIP1"]
		G=cbind(SNP=G_r[,j],x[ind_r,covs])
		l=summary(lm(y~as.matrix(G)))
		pvalr=l$coeff[2,4]
		betar=l$coeff[2,1]

		if (sign(betar)==sign(betad)){
			out_d[i,j]=pvald
			out_r[i,j]=pvalr
		} else{
			out_d[i,j]=out_r[i,j]=1
		}
	}
}

i=1
out=NULL
for (i in 1:N_snps){
	if ((sum(out_d[,i]<=alfa1 & out_r[,i]<=alfa2)/Nloops)>=TSIR_thr){
		out=c(out,i)
	}
}

l=colnames(G_d)[out]
l=unlist(strsplit(l,split="_"))[seq(1,by=2,length(l)*2)]
#[1] "rs3737240_T"   "rs9436127_A"   "rs67627084_T"  "rs13135092_G"
#[5] "rs140753832_T" "rs274623_G"    "rs967823_G"    "rs62098013_A"

ref=fread("../../../data/chronic_discovery/unification_gpc_results/GPC1/gpc1_output_done.csv",data.table=F)
table(l%in%ref$rs_id)

ref_c=ref[ref$rs_id%in%l,]





