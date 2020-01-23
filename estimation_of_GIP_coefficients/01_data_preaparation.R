######### LOAD DATA

#phenotypic correalations
phe=read.table("data/chron_corrs_disc.txt",header = T,row.names = 1)

#genotypic correalations
rgs_long=read.table("data/20181010_gene_cor_14_traits.csv",header=T,stringsAsFactors = F,sep=",")
trts=unique(c(rgs_long[,"trait_abbreviation_1"],rgs_long[,"trait_abbreviation_2"]))

#H2 values
h2_long=read.table("data/hereredity.csv",header=T,stringsAsFactors = F,sep=",")
tmp=rgs_long[rgs_long$gwas_id_1%in%h2_long$gwas_id,]
tmp=tmp[!duplicated(tmp$gwas_id_1),]
ind=match(h2_long$gwas_id,tmp$gwas_id_1)
tmp=tmp[ind,]
table(h2_long$gwas_id==tmp$gwas_id_1)

h2_long=cbind(h2_long,trait_abbreviation_1=tmp$trait_abbreviation_1)


#cleaning data
library(tidyr)


tmp=rgs_long[,c("trait_abbreviation_1","trait_abbreviation_2","rg")] %>% spread(trait_abbreviation_2,rg,fill=0)
rownames(tmp)=tmp[,1]
tmp=tmp[,-1]

tmp=cbind(t1=0,tmp)
colnames(tmp)[1]=trts[!trts%in%colnames(tmp)]

tmp=rbind(t2=0,tmp)
rownames(tmp)[1]=trts[!trts%in%rownames(tmp)]

tmp=tmp[trts,trts]

tmp_t=t(tmp)
tmp=tmp_t+tmp
diag(tmp)=1
rgs=tmp

tmp=rgs_long[,c("trait_abbreviation_1","trait_abbreviation_2","pval")] %>% spread(trait_abbreviation_2,pval,fill=0)
rownames(tmp)=tmp[,1]
tmp=tmp[,-1]

tmp=cbind(t1=0,tmp)
colnames(tmp)[1]=trts[!trts%in%colnames(tmp)]

tmp=rbind(t2=0,tmp)
rownames(tmp)[1]=trts[!trts%in%rownames(tmp)]

tmp=tmp[trts,trts]

tmp_t=t(tmp)
tmp=tmp_t+tmp
diag(tmp)=0
pvls=tmp

ind=c(8,9,10,11)
rgs=rgs[ind,ind]
pvls=pvls[ind,ind]
library(corrplot)
to_plot=as.matrix(rgs)
corrplot(corr = to_plot,method = "square",tl.col="black",p.mat = to_plot,insig = "p-value", sig.level = -1,tl.srt = 75)

#rgs[rgs>1]=1
#rgs[pvls>1e-3]=0
#rgs=as.matrix(rgs)


#l=hclust(d = as.dist(1-rgs))
#plot(l)
#to_plot=rgs[l$order[-(1:2)],l$order[-(1:2)]]

#library(corrplot)
#corrplot(corr = to_plot,method = "square",tl.col="black",p.mat = to_plot,insig = "p-value", sig.level = -1,tl.srt = 75)

phe=phe[c(1,2,6,3),c(1,2,6,3)]

colnames(phe)
h2_long$trait_abbreviation_1

colnames(rgs)
ind=c(2,3,4,1)
rgs=rgs[ind,ind]

traits=colnames(phe)

colnames(rgs)=rownames(rgs)=colnames(pvls)=rownames(pvls)=colnames(phe)


h2=h2_long[,"h2"]
names(h2)=colnames(rgs)


# caluclating covariances
rgs_cov=rgs*(sqrt(h2)%o%sqrt(h2))

env_cov=phe-rgs_cov
vare=diag(as.matrix(env_cov))
env_cor=env_cov*(sqrt(1/vare)%o%sqrt(1/vare))


#CMD SCALE
#library(plotly)


l=cmdscale(d=1-abs(rgs),k=2)
plot(x=l[,1],y=l[,2],xlab = "PC1",ylab ="PC2",main = "MDS using genetic correlation 
     matrix for seven pain types",asp=1)
text(x=l[,1],y=l[,2],labels = rownames(l),adj = c(-0.1,-0.3))



#### SAVING
rgs_cov=as.matrix(rgs_cov)
phe=as.matrix(phe)
save(file="data/gcov_phe_matrices.RData",list=c("phe","rgs_cov","pvls","h2_long"))

