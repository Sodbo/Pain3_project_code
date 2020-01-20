library(data.table)
library(dplyr)

gpc_disc_1 <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC1/gpc1_output_done.csv', data.table = F)
gpc_disc_2 <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC2/gpc2_output_done.csv', data.table = F)
gpc_disc_3 <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC3/gpc3_output_done.csv', data.table = F)
gpc_disc_4 <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_gpc_results/GPC4/gpc4_output_done.csv', data.table = F)

gpc_ma_eur_1 <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC1/gpc1_output_done.csv', data.table = F)
gpc_ma_eur_2 <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC2/gpc2_output_done.csv', data.table = F)
gpc_ma_eur_3 <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC3/gpc3_output_done.csv', data.table = F)
gpc_ma_eur_4 <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/unification_results/GPC4/gpc4_output_done.csv', data.table = F)

back_disc <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Back/Back_output_done.tsv', data.table = F)
neck_disc <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Neck/Neck_output_done.tsv', data.table = F)
knee_disc <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Knee/Knee_output_done.tsv', data.table = F)
hip_disc <- fread(file = '/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Hip/Hip_output_done.tsv', data.table = F)

back_ma_eur <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/unification_results/Back/Back_output_done.csv', data.table = F)
neck_ma_eur <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/unification_results/Neck/Neck_output_done.csv', data.table = F)
knee_ma_eur <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/unification_results/Knee/Knee_output_done.csv', data.table = F)
hip_ma_eur <- fread(file = '/mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/unification_results/Hip/Hip_output_done.csv', data.table = F)

# Working with ST1 and sorting other GWAS

st1_1 <- fread(file = '/mnt/polyomica/projects/mv_gwas/funk_an/19_loci_gpc.txt', data.table = F)

ind <- match(st1_1$SNP, gpc_disc_1$rs_id)
gpc_disc_1 <- gpc_disc_1[ind,]

ind <- match(st1_1$SNP, gpc_disc_2$rs_id)
gpc_disc_2 <- gpc_disc_2[ind,]

ind <- match(st1_1$SNP, gpc_disc_3$rs_id)
gpc_disc_3 <- gpc_disc_3[ind,]

ind <- match(st1_1$SNP, gpc_disc_4$rs_id)
gpc_disc_4 <- gpc_disc_4[ind,]

ind <- match(st1_1$SNP, gpc_ma_eur_1$rs_id)
gpc_ma_eur_1 <- gpc_ma_eur_1[ind,]

ind <- match(st1_1$SNP, gpc_ma_eur_2$rs_id)
gpc_ma_eur_2 <- gpc_ma_eur_2[ind,]

ind <- match(st1_1$SNP, gpc_ma_eur_3$rs_id)
gpc_ma_eur_3 <- gpc_ma_eur_3[ind,]

ind <- match(st1_1$SNP, gpc_ma_eur_4$rs_id)
gpc_ma_eur_4 <- gpc_ma_eur_4[ind,]

ind <- match(st1_1$SNP, back_disc$rs_id)
back_disc <- back_disc[ind,]

ind <- match(st1_1$SNP, neck_disc$rs_id)
neck_disc <- neck_disc[ind,]

ind <- match(st1_1$SNP, knee_disc$rs_id)
knee_disc <- knee_disc[ind,]

ind <- match(st1_1$SNP, hip_disc$rs_id)
hip_disc <- hip_disc[ind,]

ind <- match(st1_1$SNP, back_ma_eur$rs_id)
back_ma_eur <- back_ma_eur[ind,]

ind <- match(st1_1$SNP, neck_ma_eur$rs_id)
neck_ma_eur <- neck_ma_eur[ind,]

ind <- match(st1_1$SNP, knee_ma_eur$rs_id)
knee_ma_eur <- knee_ma_eur[ind,]

ind <- match(st1_1$SNP, hip_ma_eur$rs_id)
hip_ma_eur <- hip_ma_eur[ind,]

# Add information from GPC1 discovery
st1_1 <- mutate(st1_1, PVAL_DISC_GPC1 = gpc_disc_1$p)
st1_1 <- mutate(st1_1, BETA_DISC_GPC1 = gpc_disc_1$beta)
st1_1 <- mutate(st1_1, SE_DISC_GPC1 = gpc_disc_1$se)

# Add information from GPC2 discovery
st1_1 <- mutate(st1_1, PVAL_DISC_GPC2 = gpc_disc_2$p)
st1_1 <- mutate(st1_1, BETA_DISC_GPC2 = gpc_disc_2$beta)
st1_1 <- mutate(st1_1, SE_DISC_GPC2 = gpc_disc_2$se)

# Add information from GPC3 discovery
st1_1 <- mutate(st1_1, PVAL_DISC_GPC3 = gpc_disc_3$p)
st1_1 <- mutate(st1_1, BETA_DISC_GPC3 = gpc_disc_3$beta)
st1_1 <- mutate(st1_1, SE_DISC_GPC3 = gpc_disc_3$se)

# Add information from GPC4 discovery
st1_1 <- mutate(st1_1, PVAL_DISC_GPC4 = gpc_disc_4$p)
st1_1 <- mutate(st1_1, BETA_DISC_GPC4 = gpc_disc_4$beta)
st1_1 <- mutate(st1_1, SE_DISC_GPC4 = gpc_disc_4$se)

# Add information from GPC1 eur MA
st1_1 <- mutate(st1_1, PVAL_GPC1_MA_EUR = gpc_ma_eur_1$p)
st1_1 <- mutate(st1_1, BETA_GPC1_MA_EUR = gpc_ma_eur_1$beta)
st1_1 <- mutate(st1_1, SE_GPC1_MA_EUR = gpc_ma_eur_1$se)

# Add information from GPC2 eur MA
st1_1 <- mutate(st1_1, PVAL_GPC2_MA_EUR = gpc_ma_eur_2$p)
st1_1 <- mutate(st1_1, BETA_GPC2_MA_EUR = gpc_ma_eur_2$beta)
st1_1 <- mutate(st1_1, SE_GPC2_MA_EUR = gpc_ma_eur_2$se)

# Add information from GPC3 eur MA
st1_1 <- mutate(st1_1, PVAL_GPC3_MA_EUR = gpc_ma_eur_3$p)
st1_1 <- mutate(st1_1, BETA_GPC3_MA_EUR = gpc_ma_eur_3$beta)
st1_1 <- mutate(st1_1, SE_GPC3_MA_EUR = gpc_ma_eur_3$se)

# Add information from GPC4 eur MA
st1_1 <- mutate(st1_1, PVAL_GPC4_MA_EUR = gpc_ma_eur_4$p)
st1_1 <- mutate(st1_1, BETA_GPC4_MA_EUR = gpc_ma_eur_4$beta)
st1_1 <- mutate(st1_1, SE_GPC4_MA_EUR = gpc_ma_eur_4$se)

# Add information from Chronic Back pain discovery
st1_1 <- mutate(st1_1, PVAL_BACK_DISC = back_disc$p)
st1_1 <- mutate(st1_1, BETA_BACK_DISC = back_disc$beta)
st1_1 <- mutate(st1_1, SE_BACK_DISC = back_disc$se)

# Add information from Chronic Neck pain discovery
st1_1 <- mutate(st1_1, PVAL_NECK_DISC = neck_disc$p)
st1_1 <- mutate(st1_1, BETA_NECK_DISC = neck_disc$beta)
st1_1 <- mutate(st1_1, SE_NECK_DISC = neck_disc$se)

# Add information from Chronic Knee pain discovery
st1_1 <- mutate(st1_1, PVAL_KNEE_DISC = knee_disc$p)
st1_1 <- mutate(st1_1, BETA_KNEE_DISC = knee_disc$beta)
st1_1 <- mutate(st1_1, SE_KNEE_DISC = knee_disc$se)

# Add information from Chronic Hip pain discovery
st1_1 <- mutate(st1_1, PVAL_HIP_DISC = hip_disc$p)
st1_1 <- mutate(st1_1, BETA_HIP_DISC = hip_disc$beta)
st1_1 <- mutate(st1_1, SE_HIP_DISC = hip_disc$se)

# Add information from Chronic Back pain MA eur
st1_1 <- mutate(st1_1, PVAL_BACK_MA_EUR = back_ma_eur$p)
st1_1 <- mutate(st1_1, BETA_BACK_MA_EUR = back_ma_eur$beta)
st1_1 <- mutate(st1_1, SE_BACK_MA_EUR = back_ma_eur$se)

# Add information from Chronic Neck pain MA eur
st1_1 <- mutate(st1_1, PVAL_NECK_MA_EUR = neck_ma_eur$p)
st1_1 <- mutate(st1_1, BETA_NECK_MA_EUR = neck_ma_eur$beta)
st1_1 <- mutate(st1_1, SE_NECK_MA_EUR = neck_ma_eur$se)

# Add information from Chronic Knee pain MA eur
st1_1 <- mutate(st1_1, PVAL_KNEE_MA_EUR = knee_ma_eur$p)
st1_1 <- mutate(st1_1, BETA_KNEE_MA_EUR = knee_ma_eur$beta)
st1_1 <- mutate(st1_1, SE_KNEE_MA_EUR = knee_ma_eur$se)

# Add information from Chronic Hip pain MA eur
st1_1 <- mutate(st1_1, PVAL_HIP_MA_EUR = hip_ma_eur$p)
st1_1 <- mutate(st1_1, BETA_HIP_MA_EUR = hip_ma_eur$beta)
st1_1 <- mutate(st1_1, SE_HIP_MA_EUR = hip_ma_eur$se)

fwrite(st1_1, file = '/mnt/polyomica/projects/mv_gwas/funk_an/19_loci_extended.tsv', sep = '\t', dec = ',')
