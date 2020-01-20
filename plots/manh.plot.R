setwd('/mnt/polyomica/projects/mv_gwas/data/downloaded_from_GWASMAP')

library(data.table)
library(dplyr)

gpc_1 <- fread('unpigz -c GPC1_disc.tsv.gz')
gpc_2 <- fread('unpigz -c GPC2_disc.tsv.gz')
gpc_3 <- fread('unpigz -c GPC3_disc.tsv.gz')
gpc_4 <- fread('unpigz -c GPC4_disc.tsv.gz')

gpc_1 <- gpc_1[eaf<1-2e-4 & eaf > 2e-4]
gpc_2 <- gpc_2[eaf<1-2e-4 & eaf > 2e-4]
gpc_3 <- gpc_3[eaf<1-2e-4 & eaf > 2e-4]
gpc_4 <- gpc_4[eaf<1-2e-4 & eaf > 2e-4]

ldscr_intercept <- c(1.01601733091182,
  1.00091809801273,
  1.01301076500069,
  1.02066516433101)

gpc_1_cl <- gpc_1 %>% select(rs_id,chr,bp,z)
gpc_2_cl <- gpc_2 %>% select(rs_id,chr,bp,z)
gpc_3_cl <- gpc_3 %>% select(rs_id,chr,bp,z)
gpc_4_cl <- gpc_4 %>% select(rs_id,chr,bp,z)

gpc_1_cl$p_gc <- gpc_1_cl$z^2 / ldscr_intercept[1]
gpc_2_cl$p_gc <- gpc_2_cl$z^2 / ldscr_intercept[1]
gpc_3_cl$p_gc <- gpc_3_cl$z^2 / ldscr_intercept[1]
gpc_4_cl$p_gc <- gpc_4_cl$z^2 / ldscr_intercept[1]

gpc_1_cl$p_gc <- pchisq(gpc_1_cl$p_gc,df=1,lower.tail=FALSE)
gpc_2_cl$p_gc <- pchisq(gpc_2_cl$p_gc,df=1,lower.tail=FALSE)
gpc_3_cl$p_gc <- pchisq(gpc_3_cl$p_gc,df=1,lower.tail=FALSE)
gpc_4_cl$p_gc <- pchisq(gpc_4_cl$p_gc,df=1,lower.tail=FALSE)

library(qqman)

setwd('/mnt/polyomica/projects/mv_gwas/manhattan_plot')

tiff('manh_gpc_1.tiff',width=1024, height=720)

gpc_1_cl[p_gc<1e-2] %>% manhattan(.,main='GPC1',p='p_gc',chr = 'chr', bp='bp', snp = 'rs_id',
                      suggestiveline = FALSE,genomewideline = -log10(5e-8/4))

dev.off()

tiff('manh_gpc_2.tiff',width=1024, height=720)

gpc_2_cl[p_gc<1e-3] %>% manhattan(.,main='GPC2',p='p_gc',chr = 'chr', bp='bp', snp = 'rs_id',
                      suggestiveline = FALSE,genomewideline = -log10(5e-8/4))

dev.off()


tiff('manh_gpc_3.tiff',width=1024, height=720)

gpc_3_cl[p_gc<1e-2] %>% manhattan(.,main='GPC3',p='p_gc',chr = 'chr', bp='bp', snp = 'rs_id',
                      suggestiveline = FALSE,genomewideline = -log10(5e-8/4))

dev.off()


tiff('manh_gpc_4.tiff',width=1024, height=720)

gpc_4_cl[p_gc<1e-2] %>% manhattan(.,main='GPC4',p='p_gc',chr = 'chr', bp='bp', snp = 'rs_id',
                      suggestiveline = FALSE,genomewideline = -log10(5e-8/4))

dev.off()

tiff('qq_gpc_1.tiff')
qq(gpc_1_cl$p_gc, main='GPC1')
dev.off()

tiff('qq_gpc_2.tiff')
qq(gpc_2_cl$p_gc, main='GPC2')
dev.off()

tiff('qq_gpc_3.tiff')
qq(gpc_3_cl$p_gc, main='GPC3')
dev.off()

tiff('qq_gpc_4.tiff')
qq(gpc_4_cl$p_gc, main='GPC4')
dev.off()

