# Reading a table with sd for each type of pain

sd.table <- read.csv(
  file = '/mnt/polyomica/projects/mv_gwas/data/replication/Afr_MV_replication.csv', 
  header = TRUE, 
  sep = ';',  
  skip = 1,
  stringsAsFactors=FALSE,
  dec=','
) 

for (row.number in 1:nrow(sd.table)) {
  
  pain.type <- sd.table[row.number, 1]
  
  sd.trait <- sd.table[row.number, 4] 
  
  # Reading raw GWAS-file
  
  gwas.name <- paste('MV', pain.type, 'Repl_AFR_gwas.bgen.stats.txt', sep = "_")
  
  raw.gwas <- data.table::fread(
    input = paste0(
      '/mnt/polyomica/projects/mv_gwas/data/replication/afr_06092018/', 
       gwas.name
    )
  )
  
  # MAF computation & Filtering
  
  gwas.info.filtered <- raw.gwas[INFO  >= 0.7] # Filtering by info
  
  gwas.info.filtered$MAF <- pmin(gwas.info.filtered[, A1FREQ], 1 - gwas.info.filtered[, A1FREQ]) # Adding MAF column
  
  gwas.MAF.filtered <- gwas.info.filtered[MAF >= 1e-5] # Filtering by MAF
  
  # Standartization

  gwas.standart <- gwas.MAF.filtered
  
  gwas.standart$BETA <- gwas.standart$BETA / sd.trait
  
  gwas.standart$SE <- gwas.standart$SE / sd.trait
  
  # Writing an output file

  data.table::fwrite(gwas.standart, file = paste0('/mnt/polyomica/projects/mv_gwas/data/upload_2_gwasmap/AFR_repl_06092018/', gwas.name))
  
}

