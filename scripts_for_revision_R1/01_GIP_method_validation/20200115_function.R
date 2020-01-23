######## FUNCTIONS


H2=function(a,covm=covm,phem=phem){
  h=sum(covm*(a%o%a))/sum(phem*(a%o%a))
  return(h)
}

gen=function(a,covm=covm){
  h=sum(covm*(a%o%a))
  return(h)
}

phen=function(a,phem=phem){
  h=sum(phem*(a%o%a))
  return(h)
}



# phenotipic correlation of linear combination with a coeeficient with i trait
cor_yi_alfa=function(a,i,phem=phem){
  cov_gi_sum_giai=sum(phem[i,]*a)
  var_gi=phem[i,i]
  var_sum_giai=sum(phem*(a%o%a)) #gen(a)
  
  cor_g_a=cov_gi_sum_giai/sqrt(var_gi*var_sum_giai)
  return(cor_g_a)
}
# phenotipic covariance of linear combination with a coeeficient with i trait
cov_yi_alfa=function(a,i,phem=phem){
  cov_gi_sum_giai=sum(phem[i,]*a)
  var_gi=phem[i,i]
  var_sum_giai=phen(a,phem=phem) #gen(a)
  
  cor_g_a=cov_gi_sum_giai
  return(cor_g_a)
}

# phenotipic covariance a1 a2
cov_yi_a1_a2=function(a1,a2,phem=phem){
  cov_gi_sum_g_a1_a2=sum(phem*(a1%o%a2))
  cor_g_a=cov_gi_sum_g_a1_a2
  return(cor_g_a)
}



# genetic correlation of linear combination with a coeeficient a1 and a coeeficient a2
cor_gi_a1_a2=function(a1,a2,covm=covm){
  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
  var_g_a1=sum(covm*(a1%o%a1))
  var_g_a2=sum(covm*(a2%o%a2))
  
  cor_g_a=cov_gi_sum_g_a1_a2/sqrt(var_g_a1*var_g_a2)
  return(cor_g_a)
}
# genetic correlation of linear combination with a coeeficient with i trait
cor_gi_alfa=function(a,i,covm=covm){
  cov_gi_sum_giai=sum(covm[i,]*a)
  var_gi=covm[i,i]
  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
  
  cor_g_a=cov_gi_sum_giai/sqrt(var_gi*var_sum_giai)
  return(cor_g_a)
}
# genetic covaiance of linear combination with a coeeficient a1 and a coeeficient a2
cov_gi_a1_a2=function(a1,a2,covm=covm){
  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
  var_g_a1=sum(covm*(a1%o%a1))
  var_g_a2=sum(covm*(a2%o%a2))
  
  cor_g_a=cov_gi_sum_g_a1_a2
  return(cor_g_a)
}

# genetic covariacne of linear combination with a coeeficient with i trait
cov_gi_alfa=function(a,i,covm=covm){
  cov_gi_sum_giai=sum(covm[i,]*a)
  var_gi=covm[i,i]
  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
  
  cor_g_a=cov_gi_sum_giai
  return(cor_g_a)
}

H2_cov=function(a,covm=covm){
  h=sum(covm*(a%o%a))
  return(h)
}


#function to add GPC
add_gpc=function(gcovm,phem,ind=1,l=0){
  
  gcovm_orig=gcovm
  
  #lambda=diag(nrow(gcovm))*(1-l)
  #gcovm=gcovm%*%lambda
  
  
  lambda=rep((1-l),nrow(gcovm))
  lambda=sqrt(lambda)%o%sqrt(lambda)
  diag(lambda)=1
  gcovm=gcovm*lambda
  
  
  
  eigens=eigen(gcovm)$vectors
  eigens=apply(eigens,MAR=2,FUN=function(x){if (x[ind]<0) {x=-x}; x})
  
  gcovm=gcovm_orig
  
  vars=apply(eigens,MAR=2,phen,phem=phem)
  vars[vars<=0]=1
  sds=sqrt(vars)
  i=1
  for (i in 1:ncol(eigens)){
    eigens[,i]=eigens[,i]/sds[i]
  }
  
  #gen covs
  h2=diag(gcovm)
  nl=nrow(gcovm)
  
  out=array(NA,c(nl*2,nl*2))
  colnames(out)=rep("tr",nl*2)
  colnames(out)[1:nl]=colnames(gcovm)
  colnames(out)[(1+nl):(nl*2)]=paste0("GPC",1:nl)
  rownames(out)=colnames(out)
  out[1:nl,1:nl]=gcovm
  
  h2_cor=NULL
  i=1
  for (i in 1:nl){
    av=rep(0,n=nl)
    av=eigens[,i]
    out[nl+i,nl+i]=H2_cov(av,covm=gcovm)
    h2_cor=c(h2_cor,H2(av,covm=gcovm,phem=phem))
    j=1
    for (j in 1:nl){
      out[nl+i,j]=out[j,nl+i]=cov_gi_alfa(av,j,covm=gcovm)
    }
    j=1
    for (j in 1:nl){
      out[nl+i,j+nl]=out[j+nl,nl+i]=cov_gi_a1_a2(a1=av,a2=eigens[,j],covm=gcovm)
    }
  }
  
  #h2_cov_tmp=diag(out)
  #out[(1+nl):(nl*2),(1+nl):(nl*2)]=0
  #diag(out)=h2_cov_tmp
  out_cg=out
  
  #phen covs
  out=array(NA,c(nl*2,nl*2))
  colnames(out)=rep("tr",nl*2)
  colnames(out)[1:nl]=colnames(phem)
  colnames(out)[(1+nl):(nl*2)]=paste0("GPC",1:nl)
  rownames(out)=colnames(out)
  out[1:nl,1:nl]=as.matrix(phem)
  
  i=1
  for (i in 1:nl){
    av=rep(0,n=nl)
    av=eigens[,i]
    out[nl+i,nl+i]=phen(av,phem=phem)
    j=1
    for (j in 1:nl){
      out[nl+i,j]=out[j,nl+i]=cov_yi_alfa(av,j,phem=phem)
    }
    j=1
    for (j in 1:nl){
      out[nl+i,j+nl]=out[j+nl,nl+i]=cov_yi_a1_a2(av,a2=eigens[,j],phem=phem)
    }
  }
  out_cy=out
  
  vare=diag(out_cy)
  z=vare
  vare[vare<=0]=1
  cor_y=out_cy*(sqrt(1/vare)%o%sqrt(1/vare))
  diag(cor_y)[z<=0]=0
  
  vare=diag(out_cg)
  diag(out_cg)[vare<0]=0
  vare[vare<=0]=1
  
  cor_g=out_cg*(sqrt(1/vare)%o%sqrt(1/vare))

  rownames(eigens)=colnames(gcovm)
  colnames(eigens)=paste0("GPC",1:ncol(eigens))
  
  H22=diag(out_cg)/diag(out_cy)
  H22[is.nan(H22)]=0
  
  big_out=list(cov_g=out_cg,cov_y=out_cy,cor_g=cor_g,cor_y=cor_y,H2=H22,eigens=eigens)
  return(big_out)
}

