source("../../estimation_of_GIP_coefficients/00_core_functions.R")
load("../../estimation_of_GIP_coefficients/data/gcov_phe_matrices.RData")

h2=diag(rgs_cov)
rg=rgs_cov*(sqrt(1/h2)%o%sqrt(1/h2))


#back, 0.072775 (0.002636) 
#hip, 0.047745 (0.002522) 
#knee, 0.069480 (0.002602) 
#neck, 0.058585 (0.002576) 
gcta_h2_obs=c(0.072775,0.058585,0.069480,0.047745)

rgs_cov=rg*(sqrt(gcta_h2_obs%o%gcta_h2_obs))

four_pains=add_gpc(gcovm=as.matrix(rgs_cov),phem=as.matrix(phe),l=0)
four_pains$H2

#back        neck        knee         hip        GPC1        GPC2        GPC3        GPC4 
#0.072775000 0.058585000 0.069480000 0.047745000 0.117340745 0.041156916 0.016485035 0.008196399

#GIP1, 0.114066 (0.002779)
#GIP2, 0.032393 (0.002431)
#GIP3, 0.022172 (0.002370)
#GIP4, 0.012033 (0.002270)

gcta_h2_liab=c(0.15623962743293318,0.13155969100730291,0.15118066127892726,0.146830187183969)

rgs_cov=rg*(sqrt(gcta_h2_liab%o%gcta_h2_liab))

four_pains=add_gpc(gcovm=as.matrix(rgs_cov),phem=as.matrix(phe),l=0)
four_pains$H2

#back       neck       knee        hip       GPC1       GPC2       GPC3       GPC4 
#0.15623963 0.13155969 0.15118066 0.14683019 0.27727465 0.08942826 0.04392420 0.01936686 

#LDSC liability
ldsc_h2_liab=c(0.09051907578091768,0.07125456534275786,0.09413537666452507,0.07228033934749868)

rgs_cov=rg*(sqrt(ldsc_h2_liab%o%ldsc_h2_liab))

four_pains=add_gpc(gcovm=as.matrix(rgs_cov),phem=as.matrix(phe),l=0)
four_pains$H2

#back       neck       knee        hip       GPC1       GPC2       GPC3       GPC4 
#0.09051908 0.07125457 0.09413538 0.07228034 0.15423632 0.05412014 0.02322353 0.01046312 
