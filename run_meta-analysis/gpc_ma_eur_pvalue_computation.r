for (i in c(1,2,3,4)){
  ma <- data.table::fread(input = paste0('~/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/GPC', i, '/gpc', i, '_cau1.tbl'))
  ma$P <- pchisq((ma$Effect / ma$StdErr)^2, df = 1, lower.tail = FALSE)
  data.table::fwrite(ma, file = paste0('~/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/GPC', i, '/gpc', i, '_Pval'))
  }
  
  
  