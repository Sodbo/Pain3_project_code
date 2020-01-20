library(dplyr)

sdtable <- read.table("/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/raw/cases.controls_repl.EA.txt",
                     header=T,row.names=1)
prev=sdtable[2,]/(sdtable[1,]+sdtable[2,])
SDs=sqrt(prev*(1-prev))
#SDs
traits=names(SDs)
#traits
traits=unlist(lapply(traits,FUN = function(name){paste(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)), sep="")}))
#traits

for (col.number in 1:length(traits)) {

    pain.type=traits[col.number]
    sd.trait <- as.numeric(SDs[col.number])

    cat(paste0(pain.type,";"),"SD:",sd.trait,"\n")
    # Reading raw GWAS-file

    gwas.name <- paste('MV', pain.type, 'Repl.EA_gwas.BGEN.stats.txt', sep = "_")
    raw.gwas <- data.table::fread(
        input = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/raw/',gwas.name),
        data.table=F,header=T,stringsAsFactors=F)

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
        file = paste0('/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/', gwas.name))

}

