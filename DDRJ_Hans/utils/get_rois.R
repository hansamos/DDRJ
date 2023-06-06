##### simulation of rois

simulate_data_roi = function(n=300,p=300, prop=0.85){
  
  x_sim= rnorm(p*n) %>% matrix(nrow=n,ncol=p)
  mean_sim= apply(x_sim, 2, mean)
  sd_sim= apply(x_sim, 2, sd)
  nr=nrow(x_sim)
  nc=ncol(x_sim)
  beta_true = c(1,1.2,0.8,-1.5,-1,rnorm(nc-5,mean=0,sd=0.05),1.3)
  x_full=cbind(1,scale(x_sim))        
  y_full=ifelse(x_full %*% beta_true +rnorm(nr)> 0,1,0)
  index= sample(1:n,prop*n)
  xtrain= x_full[index,]; ytrain = y_full[index]
  xtest= x_full[-index,]; ytest=y_full[-index]
  
  return(sim_data = list(x=x_full, y= y_full, xtrain=xtrain,
                         ytrain=ytrain, xtest=xtest,
                         ytest= ytest, beta_true=beta_true))
}


get_real_data_roi = function(var_sel=T){
  
  x_full =  read_csv('../data/x.csv') %>% as.matrix() 
  y_full =  read_csv('../data/y.csv') %>% as.matrix()
  
  n = nrow(y_full)
  ## split index 85/15
  ## xtrain, ytrain (used for variable selection, 
  #but one could use the full dataset for more reliable inference.
  #Actually it's better)
  index= sample(1:n,0.85*n)
  xtrain= x_full[index,]; ytrain = y_full[index]
  mean_xtrain = apply(xtrain,2,mean); sd_xtrain=apply(xtrain, 2, sd)
  xtrain= cbind(1,xtrain)
  
  ## xtest and ytest used for prediction  
  xtest= x_full[-index,] %>% scale(center = mean_xtrain, scale = sd_xtrain)
  ytest=y_full[-index]; mean(ytest==1)
  xtest = cbind(1,xtest)
  
  ## 
  xfull = cbind(1,scale(x_full))
  
  return(list(x=x_full, y=y_full,
              xtrain=xtrain,
              ytrain=ytrain, 
              xtest=xtest,
              ytest= ytest ))
}

