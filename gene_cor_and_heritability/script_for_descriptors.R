x=fread("20190212_descriptors_full.csv")
table(x$collection)
CD_IBD              CVD            Cedar        Cytokines 
               2                8            22110               41 
         GTEx_v7          Lipidom     Metabolomics NMR_metabolomics 
         1225894              356              123              110 
           OLINK           PAIN-3   SomaLogic_2017     UKB_GenAtlas 
              82               76             1124               34 
    UKB_NealeLab      Westra_eQTL    glycomics-igg           osteop 
            2419             5948               77                3 
            pain   plasma glycome         varicose 
               2              113                6 



x[(x$collection==""),"collection"]="CD_IBD"

#list_of_traits=c("CD_IBD","Cytokines","CVD","Lipidom","Metabolomics",
#"NMR_metabolomics","OLINK","SomaLogic_2017","UKB_GenAtlas","UKB_NealeLab","glycomics-igg",
#"plasma glycome","pain","PAIN-3")

list_of_traits=c("CVD","Metabolomics",
"OLINK","SomaLogic_2017","UKB_GenAtlas","UKB_NealeLab")


ind=which(x$collection%in%list_of_traits)

x=x[ind,]



write.table(x=x,file="20180211_descriptors_v1.csv",col.names=T,row.names=F,sep="\t",quote=F)


colnames(x)
 [1] "gwas_id"                  "frequency_source"        
 [3] "license"                  "authors"                 
 [5] "organism"                 "trait_type"              
 [7] "phenotyping_technology"   "tissue"                  
 [9] "molecular_domain"         "data_doi"                
[11] "access_restrictions"      "trait_name"              
[13] "trait_abbreviation"       "study_year"              
[15] "population"               "download_link"           
[17] "study_website_link"       "reference_pmid"          
[19] "reference_doi"            "n_cases"                 
[21] "n_controls"               "n_people"                
[23] "genomic_build"            "genetic_model"           
[25] "association_type"         "association_metric"      
[27] "beta_units"               "model_association"       
[29] "imputation_reference"     "gwas_array"              
[31] "trait_transformations"    "collection"              
[33] "use_terms"                "trait_variance"          
[35] "n_of_significant_p_value" "date"                    
[37] "unified_snps_count"       "ldsc_lambda_gc"          
[39] "ldsc_intercept"           "ldsc_intercept_se"       
[41] "ldsc_mean_chi2"           "ldsc_h2"                 
[43] "ldsc_h2_se"               "ldsc_h2_pval" 

ind=which((x$n_people>=10000 & x$unified_snps_count>=1e6)&(is.na(x$n_cases)|(x$n_cases>=1000&x$n_controls>=1000)))
length(ind)
y=x[ind,]
table(is.na(y[,"ldsc_h2"]))
gwas_ids=y[is.na(y[,"ldsc_h2"]),"gwas_id"]

write.table(x=gwas_ids,file="20190212_gwasids.txt",col.names=F,row.names=F,quote=F)


#############

ind=grep(x$trait_name,pattern = "Non-cancer illness code, self-reported: varicose veins")
ind2=grep(x$trait_name,pattern = "I83")
ind=c(ind,ind2)
x=x[-ind,]

ind=which((x$collection!="UKB_NealeLab") | (x$collection=="UKB_NealeLab"&(is.na(x$n_cases)|(x$n_cases>=1000&x$n_controls>=1000))))
length(ind)

x=x[ind,]
table(x$collection)

gwas_ids=x$gwas_id

write.table(x=gwas_ids,file="../../../varicose_project/20190307/smr_heidi/20190309_filtered_2222_gwasids.txt",col.names=F,row.names=F,quote=F)


####
export PROD=TRUE

g2s=$(cat ~/20190212_gwasids.txt)

for g2 in $g2s
do
echo $g2

run_ldscore \
--h2 \
--gwas-id $g2

done
####



ind=which((x$n_people>=10000 & x$unified_snps_count>=1e6)&(is.na(x$n_cases)|(x$n_cases>=2000&x$n_controls>=2000)&
	(x$ldsc_h2)>0&x$ldsc_h2_pval<=pchisq(4,df=1,low=F)))


x=fread("20190212_descriptors_full.csv",header=T,stringsAsFactors=F,data.table=F)


x[(x$collection==""),"collection"]="CD_IBD"

list_of_traits=c("CD_IBD","Cytokines","CVD","Lipidom","Metabolomics",
"NMR_metabolomics","OLINK","SomaLogic_2017","UKB_GenAtlas","UKB_NealeLab","glycomics-igg",
"plasma glycome","pain","PAIN-3")

ind=which(x$collection%in%list_of_traits)

x=x[ind,]

ind=which((x$n_people>=10000 & x$unified_snps_count>=1e6)&(is.na(x$n_cases)|(x$n_cases>=1000&x$n_controls>=1000)))
length(ind)
y=x[ind,]
table(is.na(y[,"ldsc_h2"]))



ind=which(x$n_people>=10000 & x$unified_snps_count>=1e6 & (is.na(x$n_cases)|(x$n_cases>=2000 & x$n_controls>=2000))&
	(x$ldsc_h2)>0 & x$ldsc_h2_pval<=pchisq(4,df=1,low=F))
length(ind)
gwas_ids=x[ind,"gwas_id"]
write.table(x=gwas_ids,file="20190212_gwasids_for_ldsc.txt",col.names=F,row.names=F,quote=F)

write.table(x=x,file="20190212_descriptors_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")

list_of_traits=c("CD_IBD","Cytokines","CVD","Lipidom","Metabolomics",
"NMR_metabolomics","OLINK","SomaLogic_2017","UKB_GenAtlas","UKB_NealeLab","glycomics-igg",
"plasma glycome","pain","PAIN-3")

x=fread("20190212_descriptors_v2.txt",header=T,stringsAsFactors=F,data.table=F)
ind=with(x,which((collection=="UKB_NealeLab"&(is.na(n_cases)|(n_cases>=2000 & n_controls>=2000)))|collection!="UKB_NealeLab"))
length(ind)
gwas_ids=x[ind,"gwas_id"]
write.table(x=gwas_ids,file="20190213_gwasids_for_SMR_HEIDI.txt",col.names=F,row.names=F,quote=F)




