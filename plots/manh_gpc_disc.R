setwd('/mnt/polyomica/projects/mv_gwas/manhattan_plot/')

library(data.table)
library(dplyr)
library(R.utils)
library(qqman)

ldscr_intercept <- c(1.01702829832718,
  1.01480503061896,
  1.01773786921597,
  1.0302928308865)

for(gpc in 1){
  
  gpc_file <- "/mnt/polyomica/projects/mv_gwas/data/downloaded_from_GWASMAP/GPC"
  
  gpc_file <- gpc_file %>% paste0(., gpc, "_disc.tsv.gz")
  
  gpc_tab <- fread(gpc_file)
  
  gpc_sm <- gpc_tab %>% select(rs_id,chr,bp,p)
  
  p_gc <- qchisq(gpc_sm$p,df=1,lower.tail=FALSE) / ldscr_intercept[gpc]
  
  p_gc <- pchisq(p_gc, df=1,lower.tail=FALSE)
  
  gpc_sm <- gpc_sm %>% mutate (p_gc = p_gc)
  
  gpc_sm %>% manhattan(.,chr = 'chr', bp='bp', snp = 'rs_id',
                      suggestiveline = FALSE,genomewideline = -log10(5e-8/4))
  
}