##### simulation of snps

simulate_data_snp = function(n=300,p=300, prop=0.85){
  
  size_beta =p
  size_sample =n
  
  alpha = c(1.3,1,-1.5,-1.2,rnorm(size_beta-4,mean=0,sd=0.01))
  delta= c(-1, -1.4, -1.4, -2, rnorm(size_beta-4,mean=0,sd=0.01))
  beta_true= c(1.7,alpha,delta)
  beta_true
  x=cbind(1, matrix(rbinom(size_beta*size_sample,2,0.6), ncol=size_beta, nrow=size_sample))
  repl= function(x){ifelse(x==2,-1,x)}
  x= apply(x,2,repl)
  z= cbind(x, 1- abs(x[,2:(size_beta+1)]))
  y=ifelse(z %*% beta_true +rnorm(size_sample) >0,1,0)
  
  ## index
  index= sample(1:n,prop*n)
  ztrain =z[index,] ; ytrain = y[index]
  ztest = z[-index,]; ytest =y[-index]  
  
  
  return(list(z=z, y= y, 
              ztrain=ztrain, ytrain=ytrain,
              ztest=ztest,  ytest= ytest, beta_true=beta_true))
}


get_real_data_snp = function(var_sel=T){
  
  ##data
  zn= read_csv('../data/zn.csv') 
  zn=zn[,-1]
  zn= cbind(1,zn)
  n= nrow(zn)
  nc= ncol(zn)-1
  y_full= read_csv('../data/y.csv') %>% as.matrix()
  
  ### create z , 1-|z| 
  z_full=cbind(zn,1-abs(zn[,2:(nc+1)])) %>% as.matrix() ; rm(zn)
  
  
  ## real data that we split into train and test 
  ## split index 85/15
  index= sample(1:n,0.85*n)
  length(index)
  
  ## ztrain, ytrain (used for variable selection, 
  #but one could use the full dataset for more reliable inference.
  #Actually it's better)
  ztrain= z_full[index,]; ytrain = y_full[index]
  
  
  ## ztest and ytest used for prediction  
  ztest= z_full[-index,]
  ytest=y_full[-index] 
  
  return(list(z=z_full, y=y_full,
              ztrain=ztrain,
              ytrain=ytrain, 
              ztest=ztest,
              ytest= ytest ))
}

