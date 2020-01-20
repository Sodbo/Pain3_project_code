library(data.table)
for (i in c(1:4)) {
  gpc.disc <- fread(input = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC', i, '/gpc', i, '_output_done.csv'))
  gpc.disc.short <- gpc.disc[p < 1e-5]
  summary(gpc.disc.short$p)
  fwrite(gpc.disc.short, file = paste0('/mnt/polyomica/projects/mv_gwas/funk_an/shorter_files_p1e-5/gpc', i, '_disc_short.csv'))
  gpc.eur.ma <- fread(input = paste0('/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC', i, '/gpc', i, '_output_done.csv'))
  gpc.eur.ma.short <- gpc.eur.ma[p < 1e-5]
  summary(gpc.eur.ma.short$p)
  fwrite(gpc.eur.ma.short, file = paste0('/mnt/polyomica/projects/mv_gwas/funk_an/shorter_files_p1e-5/gpc', i, '_eur_ma_short.csv'))
}

  