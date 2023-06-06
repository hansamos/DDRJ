##### simulation of snps and rois

simulate_data_rois_snps = function(n=300,p=300, prop=0.85){
  
  size_alpha = p #snp
  size_mu = p #roi
  size_sample = n
  
  alpha = c(1.3,-1,1.5,1,rep(0,size_alpha-4))
  delta= c(-1.2, -1, -1.3, -2,  rep(0,size_alpha -4))
  mu_roi= c(1.3,0,1.5,rep(0,size_mu-5),1,0)
  beta_true=c(1,mu_roi,alpha,delta)

  
  roi= matrix(rnorm(size_mu*size_sample),nrow=size_sample,ncol=size_mu)
  snp= matrix(rbinom(size_alpha*size_sample,2,0.6), ncol=size_alpha, nrow=size_sample)
  repl= function(x){ifelse(x==2,-1,x)}
  snp= apply(snp,2,repl)
  z= cbind(snp, 1- abs(snp))
  xz=cbind(1,scale(roi),z)
  y=ifelse(xz %*% beta_true +rnorm(size_sample) >0,1,0)
  
  ## index
  index= sample(1:n,prop*n)
  xztrain =xz[index,] ; ytrain = y[index]
  xztest = xz[-index,]; ytest =y[-index]  
  
  
  return(list(xz=xz, y= y, 
              xztrain=xztrain, ytrain=ytrain,
              xztest=xztest,  xytest= ytest, beta_true=beta_true))
}


get_real_data_snp = function(var_sel=T){
  
  ##data
  #snp
  zn=  read_csv('../data/zn.csv')  %>%  as.matrix()
  zn=zn[,-1]
  nsnp=ncol(zn)
  z=cbind(zn,1-abs(zn))
  dim(z)
  #roi
  x=  read_csv('../data/x.csv.csv')  %>% as.matrix()
  x=scale(x)
  nroi=ncol(x)
  #y
  y_full= read_csv('../data/y.csv') %>%  as.matrix()
  # roi snp
  xz_full=cbind(1,x,z)
  ncol(xz_full)==(1+nroi +2*nsnp)
  #
  nc=ncol(xz_full)
  nr=nrow(xz_full)
  
  ## real data that we split into train and test 
  ## split index 85/15
  index= sample(1:nr,0.85*nr)
  length(index)
  
  ## ztrain, ytrain (used for variable selection, 
  #but one could use the full dataset for more reliable inference.
  #Actually it's better)
  xztrain= xz_full[index,]; ytrain = y_full[index]
  
  ## ztest and ytest used for prediction  
  ztest= z_full[-index,]
  ytest=y_full[-index]
  
  return(list(xz=xz_full, y=y_full,
              xztrain=xztrain,
              ytrain=ytrain, 
              xztest=xztest,
              ytest= ytest ))
}

