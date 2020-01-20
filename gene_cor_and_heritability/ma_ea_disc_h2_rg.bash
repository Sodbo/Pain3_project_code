#!/bin/bash

for pain in 'Back' 'Neck' 'Knee' 'Hip' 'Head' 'Face' 'Stom'

do

/home/common/projects/ldsc_tool/munge_sumstats.py \
--sumstat /mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/01_METAL_out/"$pain"_chronic_EUR_MA_440K1.tbl \
--out /mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/02_gene_cor/"$pain"/20181012_"$pain"_EUR_MA \
--merge-alleles /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--N-col n_total \
--snp MarkerName \
--a1 Allele1 \
--a2 Allele2 \
--p P-value \
--signed-sumstats Effect,0 \
--frq Freq1 \
--ignore FreqSE,StdErr,Direction,HetISq,HetChiSq,HetDf,HetPVal

/home/common/projects/ldsc_tool/munge_sumstats.py \
--sumstat /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/"$pain"/"$pain"_output_done.tsv \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/"$pain"/20181012_"$pain"_DISC \
--merge-alleles /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/w_hm3.snplist \
--N-col n \
--snp rs_id \
--a1 ea \
--a2 ra \
--p p \
--signed-sumstats beta,0 \
--frq eaf \
--ignore gwas_id,snp_num,chr,bp,af_ref,se,z,info,af_outlier,pz_outlier

/home/common/projects/ldsc_tool/ldsc.py \
--h2 /mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/02_gene_cor/"$pain"/20181012_"$pain"_EUR_MA.sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/02_gene_cor/"$pain"/20181012_"$pain"_EUR_MA_h2

/home/common/projects/ldsc_tool/ldsc.py \
--h2 /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/"$pain"/20181012_"$pain"_DISC.sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/"$pain"/20181012_"$pain"_DISC_h2

/home/common/projects/ldsc_tool/ldsc.py \
--rg /mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/"$pain"/20181012_"$pain"_DISC.sumstats.gz,/mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/02_gene_cor/"$pain"/20181012_"$pain"_EUR_MA.sumstats.gz \
--ref-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--w-ld-chr /home/common/projects/ldsc_tool/1kg_eur/eur_w_ld_chr/ \
--out /mnt/polyomica/projects/mv_gwas/MA/chronic_ugwas_eur_ma/02_gene_cor/"$pain"/20181012_"$pain"_EUR_MA_DISC_rg

done
