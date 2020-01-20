library(dplyr)

sd.table <- read.csv(
  file = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/SA_10102018/raw/cases.controls_repl.SA.csv',
  header = TRUE,
  sep = ';',
  skip = 1,
  stringsAsFactors=FALSE,
  dec=','
)

ptpse=c('Back', 'Neck', 'Hip', 'Face', 'Stom','Knee', 'Head')
for (col.number in 2:ncol(sd.table)) {
    
    pain.type=ptpse[col.number-1]
    sd.trait <- sd.table[5, col.number]
    
    cat(pain.type,"; ","SD: ",sd.trait,"\n")	
    # Reading raw GWAS-file

    gwas.name <- paste('MV', pain.type, 'Repl.SA_gwas.BGEN.stats.txt', sep = "_")
    raw.gwas <- data.table::fread(
      input = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_replication/SA_10102018/raw/',gwas.name),
        data.table=F,
          header=T,
          stringsAsFactors=F)
    
    #
    gwas.info.filtered <- filter(raw.gwas,INFO >= 0.7) # Filtering by info
    gwas.info.filtered=mutate(gwas.info.filtered,MAF=pmin(A1FREQ,1-A1FREQ))
    gwas.MAF.filtered <- filter(gwas.info.filtered,MAF >= 5e-3) # Filtering by maf
    
    cat("Nsnps after filtering:", nrow(gwas.MAF.filtered),"\n")
    gwas.standart <- gwas.MAF.filtered
    gwas.standart$BETA <- gwas.standart$BETA / sd.trait
    gwas.standart$SE <- gwas.standart$SE / sd.trait
    
    # Writing an output file

    data.table::fwrite(gwas.standart, 
        file = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_replication/SA_10102018/scaled_filtered/', gwas.name))

}


