# Reading a table with sd for each type of pain

sd.table <- read.csv(
  file = '/mnt/polyomica/projects/mv_gwas/data/replication/As_MV_replication.csv', 
  header = TRUE, 
  sep = ';',  
  skip = 1,
  stringsAsFactors=FALSE,
  dec=','
) 

for (row.number in 1:nrow(sd.table)) {
  
  pain.type <- sd.table[row.number, 'Pain.type']
  
  sd.trait <- sd.table[row.number, 'sd_as'] 
  
  # Reading raw GWAS-file
  
  gwas.name <- paste('MV', pain.type, 'Repl_asi_gwas.bgen.stats.txt', sep = "_")
  
  raw.gwas <- data.table::fread(
    input = paste0(
      '/mnt/polyomica/projects/mv_gwas/data/replication/asi_08092018/', 
       gwas.name
    )
  )
  
  # MAF computation & Filtering
  
  gwas.filtered <- raw.gwas[INFO >= 0.7] # Filtering by info
  

  #Adding MAF column

  gwas.filtered$MAF <- pmin(gwas.filtered[,A1FREQ], 1 - gwas.filtered[,A1FREQ]) #
  
  gwas.filtered <- gwas.filtered[MAF >= 1e-5] # Filtering by MAF
  
  # Standartization
  
  gwas.filtered$BETA <- gwas.filtered$BETA / sd.trait
  
  gwas.filtered$SE <- gwas.filtered$SE / sd.trait
  
  # Writing an output file

  data.table::fwrite(gwas.filtered, file = paste0('/mnt/polyomica/projects/mv_gwas/data/upload_2_gwasmap/ASI_repl_08092018/filtered_scaled/', gwas.name))
  
}
