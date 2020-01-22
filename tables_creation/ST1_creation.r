library(data.table)
library(dplyr)

# Creating an SNP-information-part of a table

backGWAS <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Back/Back_output_done.csv', data.table = F)
cojoFile <- fread(input = '/mnt/polyomica/projects/mv_gwas/funk_an/cojo_disc/20181023_gpc1_cojo_5e-8.txt', data.table = F)
st1 <- cojoFile$SNP
st1 <- as.data.frame(st1, stringsAsFactors=FALSE)
colnames(st1) <- c('SNP')
ind <- match(st1$SNP, backGWAS$rs_id)
backGWAS <- backGWAS[ind,]
st1 <- mutate(st1, CHR = cojoFile$Chr)
st1 <- mutate(st1, POS = cojoFile$bp)
st1 <- mutate(st1, REFERENCE_ALLELE = backGWAS$ra)
st1 <- mutate(st1, EFFECTIVE_ALLELE = backGWAS$ea)

# Creating a discovery-part of the table

gpcDisck <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC1/gpc1_output_done.csv', data.table = F)
ind <- match(st1$SNP, gpcDisck$rs_id)
gpcDisck <- gpcDisck[ind,]
st1 <- mutate(st1, MAF_DISC = pmin(gpcDisck$eaf, (1 - gpcDisck$eaf)))
st1 <- mutate(st1, EAF_DISC = gpcDisck$eaf)
st1 <- mutate(st1, IMPUTATION_QUALITY_DISC = backGWAS$info)
st1 <- mutate(st1, BETA_DISC = gpcDisck$beta)
st1 <- mutate(st1, SE_DISC = gpcDisck$se)
st1 <- mutate(st1, PVAL_DISC = gpcDisck$p)
st1 <- mutate(st1, PVAL_GC_DISC = pchisq(((st1$BETA / st1$SE)^2 / 1.016), df = 1, lower.tail = F))
st1 <- mutate(st1, HWE_DISC = NA)
st1 <- mutate(st1, N_DISC = gpcDisck$n)

# Creating a replication-EUR-part of the table

gpcEUR <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_EA/unification_gpc_results/GPC1/gpc1_output_done.csv', data.table = F)
ind <- match(st1$SNP, gpcEUR$rs_id)
gpcEUR <- gpcEUR[ind,]
backEUR <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/unification_results/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backEUR$rs_id)
backEUR <- backEUR[ind,]
st1 <- mutate(st1, BETA_EUR = gpcEUR$beta)
st1 <- mutate(st1, SE_EUR = gpcEUR$se)
st1 <- mutate(st1, PVAL_EUR = gpcEUR$p)
st1 <- mutate(st1, HWE_EUR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_EUR = backEUR$info)
st1 <- mutate(st1, EAF_EUR = gpcEUR$eaf)
st1 <- mutate(st1, N_EUR = gpcEUR$n)

# Creating a replication-AFR-part of the table

gpcAFR <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_AA/20181015_REPL_GPC1.txt')
tmp <- gpcAFR[1, ]
tmp <- tmp[-1, ]
na_row <- gpcAFR[1, ]
na_row[1, ] <- NA
for (row in c(1:nrow(st1))) {
  if (nrow(gpcAFR[chr == st1$CHR[row] & pos == st1$POS[row]]) == 0) {
  tmp <- rbind(tmp, na_row)
  } else {
  tmp <- rbind(tmp, gpcAFR[chr == st1$CHR[row] & pos == st1$POS[row]])
  }
 }
for (i in c(1:nrow(tmp))) {
  if (is.na(tmp[i, ]) == FALSE) {
    if (gpcDisck$ra[i] == tmp$A2[i]) {
	  if (gpcDisck$ea[i] != tmp$A1[i]) {
	  print(paste0('Third allele in row #', i))
	  }
	} else {
	  if (gpcDisck$ra[i] == tmp$A1[i]) {
	    if (gpcDisck$ea[i] == tmp$A2[i]) {
		  print(paste0('Flip in row #', i))
		  tmp$A2[i] <- gpcDisck$ra[i]
		  tmp$A1[i] <- gpcDisck$ea[i]
		  tmp$b[i] <- (-1)*tmp$b[i]
		  tmp$eaf[i]<- 1 - tmp$eaf[i]
		} else {
		  tmp$A1[i] <- tmp$A2[i]
          tmp$A2[i] <- gpcDisck$ra[i]		  
		  print(paste0('Third allele in row #', i))
	    } 
	  } else {
		print('Error!')
	  }
	}
  }
}
gpcAFR <- tmp 
backAFR <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/AA_10102018/raw/MV_Back_Repl.AA_gwas.BGEN.stats.txt')
tmp <- backAFR[1, ]
tmp <- tmp[-1, ]
na_row <- backAFR[1, ]
na_row[1, ] <- NA
for (row in c(1:nrow(st1))) {
  if (nrow(backAFR[CHR == st1$CHR[row] & BP == st1$POS[row]]) == 0) {
  tmp <- rbind(tmp, na_row)
  } else {
  tmp <- rbind(tmp, backAFR[CHR == st1$CHR[row] & BP == st1$POS[row]])
  }
 }
for (i in c(1:nrow(tmp))) {
  if (is.na(tmp[i, ]) == FALSE) {
    if (gpcDisck$ra[i] == tmp$ALLELE0[i]) {
	  if (gpcDisck$ea[i] != tmp$ALLELE1[i]) {
	  print(paste0('Third allele in row #', i))
	  }
	} else {
	  if (gpcDisck$ra[i] == tmp$ALLELE1[i]) {
	    if (gpcDisck$ea[i] == tmp$ALLELE0[i]) {
		  print(paste0('Flip in row #', i))
		  tmp$ALLELE0[i] <- gpcDisck$ra[i]
		  tmp$ALLELE1[i] <- gpcDisck$ea[i]
		  tmp$BETA[i] <- (-1)*tmp$BETA[i]
		  tmp$A1FREQ[i]<- 1 - tmp$A1FREQ[i]
		} else {
		  tmp$ALLELE1[i] <- tmp$ALLELE0[i]
          tmp$ALLELE0[i] <- gpcDisck$ra[i]		  
		  print(paste0('Third allele in row #', i))
	    } 
	  } else {
		print('Error!')
	  }
	}
  }
}
backAFR <- tmp
st1 <- mutate(st1, BETA_AFR = gpcAFR$b)
st1 <- mutate(st1, SE_AFR = gpcAFR$se)
st1 <- mutate(st1, PVAL_AFR = gpcAFR$p)
st1 <- mutate(st1, HWE_AFR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_AFR = backAFR$INFO)
st1 <- mutate(st1, EAF_AFR = gpcAFR$eaf)
st1 <- mutate(st1, N_AFR = gpcAFR$N)

# Creating a replication-ASI-part of the table

gpcASI <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/GPC_SA/20181015_REPL_GPC1.txt')
tmp <- gpcASI[1, ]
tmp <- tmp[-1, ]
na_row <- gpcASI[1, ]
na_row[1, ] <- NA
for (row in c(1:nrow(st1))) {
  if (nrow(gpcASI[chr == st1$CHR[row] & pos == st1$POS[row]]) == 0) {
  tmp <- rbind(tmp, na_row)
  } else {
  tmp <- rbind(tmp, gpcASI[chr == st1$CHR[row] & pos == st1$POS[row]])
  }
 }
for (i in c(1:nrow(tmp))) {
  if (is.na(tmp[i, ]) == FALSE) {
    if (gpcDisck$ra[i] == tmp$A2[i]) {
	  if (gpcDisck$ea[i] != tmp$A1[i]) {
	  print(paste0('Third allele in row #', i))
	  }
	} else {
	  if (gpcDisck$ra[i] == tmp$A1[i]) {
	    if (gpcDisck$ea[i] == tmp$A2[i]) {
		  print(paste0('Flip in row #', i))
		  tmp$A2[i] <- gpcDisck$ra[i]
		  tmp$A1[i] <- gpcDisck$ea[i]
		  tmp$b[i] <- (-1)*tmp$b[i]
		  tmp$eaf[i]<- 1 - tmp$eaf[i]
		} else {
		  tmp$A1[i] <- tmp$A2[i]
          tmp$A2[i] <- gpcDisck$ra[i]		  
		  print(paste0('Third allele in row #', i))
	    } 
	  } else {
		print('Error!')
	  }
	}
  }
}
gpcASI <- tmp
backASI <- fread(input = '/mnt/polyomica/projects/mv_gwas/data/chronic_replication/SA_10102018/raw/MV_Back_Repl.SA_gwas.BGEN.stats.txt')
tmp <- backASI[1, ]
tmp <- tmp[-1, ]
na_row <- backASI[1, ]
na_row[1, ] <- NA
for (row in c(1:nrow(st1))) {
  if (nrow(backASI[CHR == st1$CHR[row] & BP == st1$POS[row]]) == 0) {
  tmp <- rbind(tmp, na_row)
  } else {
  tmp <- rbind(tmp, backASI[CHR == st1$CHR[row] & BP == st1$POS[row]])
  }
 }
for (i in c(1:nrow(tmp))) {
  if (is.na(tmp[i, ]) == FALSE) {
    if (gpcDisck$ra[i] == tmp$ALLELE0[i]) {
	  if (gpcDisck$ea[i] != tmp$ALLELE1[i]) {
	  print(paste0('Third allele in row #', i))
	  }
	} else {
	  if (gpcDisck$ra[i] == tmp$ALLELE1[i]) {
	    if (gpcDisck$ea[i] == tmp$ALLELE0[i]) {
		  print(paste0('Flip in row #', i))
		  tmp$ALLELE0[i] <- gpcDisck$ra[i]
		  tmp$ALLELE1[i] <- gpcDisck$ea[i]
		  tmp$BETA[i] <- (-1)*tmp$BETA[i]
		  tmp$A1FREQ[i]<- 1 - tmp$A1FREQ[i]
		} else {
		  tmp$ALLELE1[i] <- tmp$ALLELE0[i]
          tmp$ALLELE0[i] <- gpcDisck$ra[i]		  
		  print(paste0('Third allele in row #', i))
	    } 
	  } else {
		print('Error!')
	  }
	}
  }
}
backASI <- tmp
st1 <- mutate(st1, BETA_ASI = gpcASI$b)
st1 <- mutate(st1, SE_ASI = gpcASI$se)
st1 <- mutate(st1, PVAL_ASI = gpcASI$p)
st1 <- mutate(st1, HWE_ASI = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_ASI = backASI$INFO)
st1 <- mutate(st1, EAF_ASI = gpcASI$eaf)
st1 <- mutate(st1, N_ASI = gpcASI$N)

# Creating an MA-replication-part of the table

gpcMA <- fread(input = '/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_repl/METAL_out/GPC1/gpc1_rep_raw1.tbl')
tmp <- NULL
for (row in c(1:nrow(st1))) {
  if (nrow(gpcMA[MarkerName == st1$SNP[row]]) == 0) {
  tmp <- rbind(tmp, gpcMA[MarkerName == paste0(st1$CHR[row], ':', st1$POS[row], '_', st1$REFERENCE_ALLELE[row], '_', st1$EFFECTIVE_ALLELE[row])])
  } else {
  tmp <- rbind(tmp, gpcMA[MarkerName == st1$SNP[row]])
  }
 }
for (i in c(1:nrow(tmp))) { 
  if (tolower(gpcDisck$ra[i]) == tmp$Allele2[i]) {
	if (tolower(gpcDisck$ea[i]) != tmp$Allele1[i]) {
	print(paste0('Third allele in row #', i))
	}
  } else {
    if (tolower(gpcDisck$ra[i]) == tmp$Allele1[i]) {
	  if (tolower(gpcDisck$ea[i]) == tmp$Allele2[i]) {
	  print(paste0('Flip in row #', i))
	  tmp$Allele2[i] <- tolower(gpcDisck$ra[i])
	  tmp$Allele1[i] <- tolower(gpcDisck$ea[i])
	  tmp$Effect[i] <- (-1)*tmp$Effect[i]
	  tmp$Freq1[i]<- 1 - tmp$Freq1[i]
	  } else {
		tmp$Allele1[i] <- tmp$Allele2[i]
        tmp$Allele2[i] <- tolower(gpcDisck$ra[i])		  
		print(paste0('Third allele in row #', i))
	  } 
	} else {
	  print('Error!')
	}
  }
}
gpcMA <- tmp
st1 <- mutate(st1, EAF_MA_REPL = gpcMA$Freq1)
st1 <- mutate(st1, EAF_SE = gpcMA$FreqSE)
st1 <- mutate(st1, BETA_MA_REPL = gpcMA$Effect)
st1 <- mutate(st1, SE_MA_REPL = gpcMA$StdErr)
st1 <- mutate(st1, PVAL_MA_REPL = gpcMA$'P-value')
st1 <- mutate(st1, Direction = gpcMA$Direction)
st1 <- mutate(st1, HetISq = gpcMA$HetISq)
st1 <- mutate(st1, HetChiSq = gpcMA$HetChiSq)
st1 <- mutate(st1, HetDf = gpcMA$HetDf)
st1 <- mutate(st1, HetPVal = gpcMA$HetPVal)
st1 <- mutate(st1, N_total = gpcMA$n_total)

fwrite(st1, file = '/mnt/polyomica/projects/mv_gwas/funk_an/ST1_GPC1.tsv', sep = '\t', dec = ',')









