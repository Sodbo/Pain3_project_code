for(gpc in 1:4){

	gpc_in <- paste0('20181011_DISC_GPC',gpc,'.txt')

	tmp <- data.table::fread(gpc_in)

	print(paste('Read complete GPC',gpc))

	tmp <- tmp[,c('SNP','A1','A2','eaf','b','se','p','N')]

	colnames(tmp) <- c('SNP','A1','A2','freq','b','se','p','N')

	gpc_out <- paste0('20181011_DISC_GPC',gpc,'.for_cojo')

	data.table::fwrite(tmp, file = gpc_out, sep = '\t')

	print(paste('Finished GPC',gpc))

}