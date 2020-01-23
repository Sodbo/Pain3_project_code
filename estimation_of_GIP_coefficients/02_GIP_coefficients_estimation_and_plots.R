#GIP coeffisients estiamtion and plots

source("00_core_functions.R")

load("data/gcov_phe_matrices.RData")

four_pains=add_gpc(gcovm=as.matrix(rgs_cov),phem=as.matrix(phe),l=0)

save(file="data/four_pains_GIPs.RData",list="four_pains")
#eigen(rgs_cov)$values/sum(eigen(rgs_cov)$values)
#[1] 0.78443550 0.14670590 0.04449125 0.02436735

########
out=four_pains
out_rg=out$cor_g
diag(out_rg)=out$H2
library(corrplot)
toPlot=out_rg


#pdf("tmp.pdf",width = 7,height = 7)
  corrplot(corr = toPlot,method = "square",
         tl.col="black",p.mat = toPlot,insig = "p-value", sig.level = -1,tl.srt = 75)
#dev.off()

toPlot=four_pains$cov_y
#pdf("tmp4.pdf",width = 7,height = 7)
corrplot(corr = toPlot,method = "square",
         tl.col="black",p.mat = toPlot,insig = "p-value", sig.level = -1,tl.srt = 75)
#dev.off()


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



