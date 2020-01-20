library(dplyr)

sd.table <- read.csv(
  file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/cases_controls_chronic.csv',
  header = TRUE,
  sep = ';',
  skip = 0,
  stringsAsFactors=FALSE,
  dec=','
)

ptpse=c('Back', 'Neck', 'Hip', 'Face', 'Stom','Knee', 'Head')
for (col.number in 2:ncol(sd.table)) {
    
    pain.type=ptpse[col.number-1]
    sd.trait <- sd.table[5, col.number]
    
    cat(pain.type,"; ","SD: ",sd.trait,"\n")	
    # Reading raw GWAS-file

    gwas.name <- paste('MV', pain.type, 'Disc_gwas.BGEN.stats.txt', sep = "_")
    raw.gwas <- data.table::fread(
      input = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/raw/',gwas.name),
        data.table=F,
          header=T,
          stringsAsFactors=F)
    
    #
    gwas.info.filtered <- filter(raw.gwas,INFO >= 0.7) # Filtering by info
    gwas.info.filtered=mutate(gwas.info.filtered,MAF=pmin(A1FREQ,1-A1FREQ))
    gwas.MAF.filtered <- filter(gwas.info.filtered,MAF >= 1e-5) # Filtering by info
    
    cat("Nsnps after filtering:", nrow(gwas.MAF.filtered),"\n")
    gwas.standart <- gwas.MAF.filtered
    gwas.standart$BETA <- gwas.standart$BETA / sd.trait
    gwas.standart$SE <- gwas.standart$SE / sd.trait
    
    # Writing an output file

    data.table::fwrite(gwas.standart, 
        file = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/scaled_filtered/', gwas.name))

}


