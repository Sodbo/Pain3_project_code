library(data.table)
library(dplyr)

back <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_Back_Repl.EA_gwas.BGEN.stats.txt',data.table=F)
neck <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_Neck_Repl.EA_gwas.BGEN.stats.txt',data.table=F)
knee <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_Knee_Repl.EA_gwas.BGEN.stats.txt',data.table=F)
hip <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_Hip_Repl.EA_gwas.BGEN.stats.txt',data.table=F)

load('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/20181010_GPCs.RData')

aa <- as.matrix(four_pains$eigens)
covm <- fread('/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/raw/chron_corrs_repl.EA.txt',
  data.table = FALSE)

row.names(covm) <- covm$V1

covm <- covm[c('back','neck','knee','hip'),c('back','neck','knee','hip')]

covm <- as.matrix(covm)

snps=back$SNP
snps=intersect(neck$SNP,snps)
snps=intersect(hip$SNP,snps)
snps=intersect(knee$SNP,snps)

ind=match(snps,back$SNP)
back=back[ind,]
table(back$SNP==snps)

ind=match(snps,neck$SNP)
neck=neck[ind,]
table(neck$SNP==snps)

ind=match(snps,knee$SNP)
knee=knee[ind,]
table(knee$SNP==snps)

ind=match(snps,hip$SNP)
hip=hip[ind,]
table(hip$SNP==snps)

betas=select(back,back_b=BETA)	
betas=mutate(betas, neck_b=neck$BETA, knee_b=knee$BETA, hip_b=hip$BETA)	
betas <- as.matrix(betas)

se1 <- back$SE

GWAS_linear_combination=function(a,betaa,se1,vary1=1,covm,N){
    vary2=sum(covm*(a%o%a))
    b=(betaa%*%a)
    
    vary1_varg_n=se1^2+(betaa[,1]^2)/N
    
    se2=vary1_varg_n*(vary2/vary1)
    
    se2=se2-b^2/N
    
    se=sqrt(se2)
    
    b=b/sqrt(vary2)
    se=se/sqrt(vary2)
    
    out=cbind(b=b,se=se)
    colnames(out)=c("b","se")
    out=as.data.frame(out)
    return(out)
}

for (NPC in 1:4){
  x=GWAS_linear_combination(a=as.vector(aa[,NPC]),betaa=as.matrix(betas),se1=as.vector(se1),vary1=1,covm=as.matrix(covm),N=174831)
  x=mutate(x,Z=b/se,p=pchisq(Z^2,1,low=F))
  x=mutate(x,SNP=back$SNP)
  x=mutate(x,A1=back$ALLELE1,A2=back$ALLELE0,N=174831,chr=back$CHR,pos=back$BP,
             eaf=back$A1FREQ)
  
  head(x,n=2)
  
  fnme=paste0("/mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/20181011_REPL_RAW_GPC",NPC,".txt")
    
  data.table::fwrite(
    x, 
    row.names=F,
    file = fnme,
    sep = '\t')
  }
  
