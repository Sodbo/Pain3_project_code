######### LOAD DATA

#setwd("/mnt/Disk_D/DB_ux/Dropbox/MV_PAIN/results/20180928/data/")
#setwd("C:/Users/ND/Dropbox/MV_PAIN/results/20181010/data")
#setwd("Dropbox/MV_PAIN/results/20181010/data")
setwd("~/Dropbox/Dropbox/MV_PAIN/results/20200115/02_indep_traits_data/")
source("~/Dropbox/Dropbox/MV_PAIN/results/20200115/20200115_function.R")


covm<-read.table('pheno_corr_matrix.txt', row.names=1,  check.names=F)
gcovm<-read.table('gene_cov_matrix.txt', row.names=1,  check.names=F)
#gcovm_se<-read.table('gene_cov_se_matrix.txt', row.names=1,  check.names=F)


nms=c("HBMD","VVs","LH","TL_VLDL")
names(nms)=c(4128,4041,1112,377)
phe=as.matrix(covm)
rgs_cov=as.matrix(gcovm)
colnames(rgs_cov)=rownames(rgs_cov)=colnames(phe)=rownames(phe)=nms

#rgs_cov=diag(nrow = 4)*0.5
#phe=diag(nrow=4)

#rgs_cov=matrix(1,nrow=4,ncol=4)*sqrt(rep(0.5,4)%o%rep(0.5,4))
#phe=matrix(1,nrow=4,ncol=4)
#colnames(rgs_cov)=rownames(rgs_cov)=colnames(phe)=rownames(phe)=paste0("Trait",1:4)

four_pains=add_gpc(gcovm=as.matrix(rgs_cov),phem=as.matrix(phe),l=0)

out=four_pains
out_rg=out$cor_g
diag(out_rg)=out$H2
library(corrplot)
toPlot=out_rg

colnames(toPlot)[5:8]=rownames(toPlot)[5:8]=paste0("GIP",1:4)
corrplot(corr = toPlot,method = "square",
         tl.col="black",p.mat = toPlot,insig = "p-value", sig.level = -1,tl.srt = 75)

toPlot=four_pains$cov_y
colnames(toPlot)[5:8]=rownames(toPlot)[5:8]=paste0("GIP",1:4)
corrplot(corr = toPlot,method = "square",
         tl.col="black",p.mat = toPlot,insig = "p-value", sig.level = -1,tl.srt = 75)



library(dplyr)
library(ggplot2)
library(tidyr)
dta=out$cor_g[(nrow(out$cor_g)/2+1):nrow(out$cor_g),1:(nrow(out$cor_g)/2)]
dta=dta^2
dta=cbind(dta,GPC=rownames(dta))
dta=as.data.frame(dta)
dta[,1:(nrow(out$cor_g)/2)]=apply(dta[,1:(nrow(out$cor_g)/2)],MAR=2,as.numeric)
dta=dta %>% gather(Trait, r2,-GPC)
#dta=cbind(dta,COL=dta$GPC)
#dta$COL=as.factor(dta$COL)
#levels(dta$COL)=c("red","blue","darkgreen","yellow")

#pdf("tmp2.pdf",width = 7,height = 7)
dta %>%
  ggplot(aes(x = Trait, y = r2, fill = GPC)) + 
  labs(title="Genetic variance explained by GPC in pain traits",
       x ="Pain type", y = "Genetic variance explained") +
  scale_fill_brewer(palette="Paired") +
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()


dta=out$eigens
colnames(dta)=paste0("GPC",1:4)
rownames(dta)=colnames(out$cov_g)[1:4]

dta=cbind(dta,Trait=rownames(dta))
dta=as.data.frame(dta)
dta[,1:4]=apply(dta[,1:4],MAR=2,as.numeric)
dta=dta %>% gather(GPC, a,-Trait)

dta$GPC=as.factor(dta$GPC)
dta$GPC=factor(dta$GPC,levels(dta$GPC)[4:1])




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
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

#pdf("tmp3.pdf",width = 7,height = 7)
l
#dev.off()



