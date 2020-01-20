library(qqman)
for (i in c(1:4)) {
  gwasResults <- data.table::fread(input = paste0("/mnt/polyomica/projects/mv_gwas/funk_an/shorter_files_p1e-5/gpc", i, "_disc_short.csv"), data.table = F)
  gwasResults <- gwasResults[which(gwasResults$eaf > 0.0002 & gwasResults$eaf < 0.9998), ]
  png(file = paste0("/mnt/polyomica/projects/mv_gwas/funk_an/plots/GPC", i, "_disc_manh.png"), height=720, width=1080)
  manhattan(gwasResults, chr="chr", bp="bp", snp="rs_id", p="p", main = paste0("Manhattan Plot GPC", i, " discovery"), ylim = c(5, max(-log10(gwasResults$p)) + 2), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = -log10(0.00000005))
  dev.off()
  gwasFull <- data.table::fread(input = paste0("/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC", i, ".txt"), data.table = F)
  gwasFull <- gwasFull[which(gwasFull$eaf > 0.0002 & gwasFull$eaf < 0.9998), ]
  png(file = paste0("/mnt/polyomica/projects/mv_gwas/funk_an/plots/GPC", i, "_disc_qq.png"), height=720, width=1080)
  qq(gwasFull$p, main = paste0("Q-Q plot for GPC", i, " discovery"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
  dev.off()
  }
