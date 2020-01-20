for (pain in c('Back', 'Neck', 'Knee', 'Hip', 'Head', 'Face', 'Stom')){
  ma <- data.table::fread(input = paste0('~/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/01_METAL_out/', pain, '_chronic_EUR_MA_440K1.tbl'))
  ma$P <- pchisq((ma$Effect / ma$StdErr)^2, df = 1, lower.tail = FALSE)
  data.table::fwrite(ma, file = paste0('~/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/01_METAL_out/', pain, '_chronic_EUR_MA_440K1_Pval.tbl'))
  }
  