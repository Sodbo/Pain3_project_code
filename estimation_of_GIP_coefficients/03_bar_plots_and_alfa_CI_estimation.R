#### Alfa PLOTs and CI

load("data/four_pains_GIPs.RData")
load("data/gcov_phe_matrices.RData")
source("00_core_functions.R")
out=four_pains

library(dplyr)
library(ggplot2)
library(tidyr)

################## BAR PLOTS
dta=out$eigens
colnames(dta)=paste0("GPC",1:4)
rownames(dta)=colnames(out$cov_g)[1:4]

dta=cbind(dta,Trait=rownames(dta))
dta=as.data.frame(dta)
dta[,1:4]=apply(dta[,1:4],MAR=2,as.numeric)
dta=dta %>% gather(GPC, a,-Trait)

dta$GPC=as.factor(dta$GPC)
dta$GPC=factor(dta$GPC,levels(dta$GPC)[4:1])


###########
### PERMUTATIONS
Chi2=qchisq(as.matrix(pvls),1,low=F)
Z=sqrt(Chi2)
se=rgs_cov/Z
diag(se)=h2_long[,"h2_se"]

phem=phe
gcovm=rgs_cov

noise_sd=function(se){
  n=nrow(se)
  out=array(NA,c(n,n))
  for (i in 1:(n)){
    j=i
    for (j in (i):n){
      out[i,j]=rnorm(1,mean=0,sd=se[i,j])
    }
  }
  #diag(out)=0
  out[lower.tri(out)]=0
  out_1=t(out)
  out=out+out_1
  out
}

comp=1
N_permut=1000
lll=0
for (comp in 1:4){
  
  out=array(NA,c(N_permut,nrow(gcovm)))
  dim(out)
  dim(rgs_cov)
  i=1
  
  for (i in 1:N_permut){
    NS=noise_sd(se)
    Z_p=gcovm+NS
    
    
    lambda=rep((1-lll),nrow(Z_p))
    lambda=sqrt(lambda)%o%sqrt(lambda)
    diag(lambda)=1
    Z_p=Z_p*lambda
    
    
    
    eigens=eigen(Z_p)$vectors
    eigens=apply(eigens,MAR=2,FUN=function(x){if (x[1]<0) {x=-x}; x})
    
    vars=apply(eigens,MAR=2,phen,phem=phem)
    sds=sqrt(vars)
    j=1
    for (j in 1:ncol(eigens)){
      eigens[,j]=eigens[,j]/sds[j]
    }
    
    
    
    l=eigens[,comp]
    
    out[i,]=l
  }
  
  j=1
  for (j in 1:4){
    qq=quantile(x=out[,j],probs = c(0.025,0.975))
    dta[(comp-1)*4+j,"CI"]=abs(qq[1]-qq[2])/2
  }
}



### Brplots ploting
library(ggplot2)

dodge <- position_dodge(width = 1)

l=ggplot(dta, aes(x=GPC, y=a,fill=Trait))+
  geom_bar(stat="identity",position = position_dodge()) +
  xlab("GPC") +
  ylab("Input of each pain trait into GPC") +
  coord_flip()+
  scale_fill_brewer(palette="Spectral")+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymax = a + CI, ymin = a-CI), position = dodge, width = 0.2)

#pdf("tmp3.pdf",width = 7,height = 7)
l
#dev.off()

