library(data.table)
library(dplyr)

back <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Back/Back_output_done.csv',data.table=F)
neck <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Neck/Neck_output_done.csv',data.table=F)
knee <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Knee/Knee_output_done.csv',data.table=F)
hip <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Hip/Hip_output_done.csv',data.table=F)

load('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/20181010_GPCs.RData')
aa <- as.matrix(four_pains$eigens)
covm <- as.matrix(four_pains$cov_y)
covm <- covm[1:4, 1:4]

snps=back$rs_id
snps=intersect(neck$rs_id,snps)
snps=intersect(hip$rs_id,snps)
snps=intersect(knee$rs_id,snps)

ind=match(snps,back$rs_id)
back=back[ind,]
table(back$rs_id==snps)

ind=match(snps,neck$rs_id)
neck=neck[ind,]
table(neck$rs_id==snps)

ind=match(snps,knee$rs_id)
knee=knee[ind,]
table(knee$rs_id==snps)

ind=match(snps,hip$rs_id)
hip=hip[ind,]
table(hip$rs_id==snps)

betas=select(back,back_b=beta)	
betas=mutate(betas,neck_b=neck$beta,knee_b=knee$beta,hip_b=hip$beta)	
betas <- as.matrix(betas)

se1 <- back$se

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
  x=GWAS_linear_combination(a=as.vector(aa[,NPC]),betaa=as.matrix(betas),se1=as.vector(se1),vary1=1,covm=as.matrix(covm),N=265000)
  x=mutate(x,Z=b/se,p=pchisq(Z^2,1,low=F))
  x=mutate(x,SNP=back$rs_id)
  x=mutate(x,A1=back$ea,A2=back$ra,N=265000,chr=back$chr,pos=back$bp,
             eaf=back$eaf)
  
  head(x,n=2)
  
  fnme=paste0("/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC",NPC,".txt")
    
  data.table::fwrite(
    x, 
    row.names=F,
    file = fnme,
    sep = '\t')
  }
  
