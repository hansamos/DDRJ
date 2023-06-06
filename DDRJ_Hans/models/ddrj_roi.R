#############################################
#############################################
# DDRJ for ROIs selection.                 ##
# Provided a design matrix with continous  ##
# covariates and bernoulli outcome,        ##  
# this code performs covariate selection   ##
# using a data driven reversible jump using## 
# a Bayesian probit model.                 ##
# Prediction can also be done using Bayesian#
# model Averaging.                         ##
#############################################

#libraries
library(dplyr)
library(readr)
library(mvtnorm)
library(extraDistr) 
library(ROCR)
library(glmnet)

##utils function to simulate or get real data 
source("../utils/get_rois.R")

###################################

##data
data_sim = simulate_data_roi(n=300, p=300, prop=0.85)

## let's do variable selection, so we use the full data
x= data_sim$x
y= data_sim$y

## for prediction using .85 train test 
#x= data_sim$xtrain
#y= data_sim$ytrain
#xtest= data_sim$xtest
#ytest = data_sim$ytest

## true beta for comparison
beta_true= data_sim$beta_true 
sum(abs(beta_true)>0.5) ## significant effect + intercept 

###################################
# Functions
###################################



### draw a value from a  discrete distribution
#
rdiscrete<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val}


### Decide on birth, death 
#p=: length(current beta)-1 ( to exclude beta_0), nc: number of columns
# return probability and type of move
pbirth_death = function(p){
  if(p>=0  & p<= nc-1){ 
  move=ifelse(p==0,1,
              ifelse( p > 0 & p < nc-1, rbinom(1,1,0.5)+1,2))
  pbirthdeath=ifelse(p==0 | p==nc-1, 1, 0.5)
  list(pbirthdeath,move)}else{stop("p > col_number")}
  #move 1= birth, 2=death
}



### compute birth pbj probability and  candidate rois to be included
#x: matrix of rois
#beta: current coefficients beta
#latent : z 
birth_roi = function(x,beta,latent){
  pb_size_beta=rep(NA,nc)
  cand= setdiff(which(beta==0),1) #candidates
  res= latent - (x %*% beta) #compute current residual
  corr=cor(res,x[,cand]) %>%  abs()
  pbj_geral=corr/sum(corr)
  #pbj_geral= rep(1/length(cand),length(cand)) #rj with uniform move 
  pb_size_beta[cand]=pbj_geral
  roi_candidate_birth= rdiscrete(pbj_geral)
  insert_roi=cand[roi_candidate_birth]
  pbj= pb_size_beta[insert_roi]
  list(pbj,insert_roi,pb_size_beta)
}

### compute death pdj probability and candidate rois to be deleted

#x: matrix of rois
#beta: current coefficients beta
#latent: current latent variable z
death_roi = function(x, beta, latent){
  #browser()
  ## number of rois candidate to be deleted
  pd_size_beta=rep(NA,nc)
  cand=setdiff(which(beta!=0),1)
  xi=numeric(length(cand))
  for (i in 1:length(cand)) {
   xi[i]= 1/abs(beta[cand[i]])
  }
  pdj_geral= xi/sum(xi)
# pdj_geral= rep(1/length(cand),length(cand)) #rj with uniform move
  pd_size_beta[cand]=pdj_geral
  roi_candidate_death = rdiscrete(pdj_geral)
  remove_roi= cand[roi_candidate_death]
  pdj= pd_size_beta[remove_roi]
  list(pdj, remove_roi,pd_size_beta)
}


#latent: current latent variable
### Latent variable
ztr=function(x,y,beta){
  x=as.matrix(x)
 # beta=matrix(beta,ncol(x),1,byrow = T)
  return(ifelse(y==1,rtnorm(nr,x%*%beta,1,a=0),rtnorm(nr,x%*%beta,1,b=0)))
}



### acceptance probability of birth 
acceptance_birth_roi = function(x,y,beta_ant,beta_new,pbj,pdj){
  
## numerator 
#loglik
loglik_num= dmvnorm(x= as.numeric(latent[i-1,]- x%*%beta_new_full),mean=rep(0,nr),
                    sigma =diag(nr),log = TRUE)

#prior
prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                    sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)


#transition
pbirthdeath_num=pbirth_death(length(which(beta_new[-1]!=0)))[[1]]


roi_ant= which(beta_ant!=0)
vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(roi_ant)) + 
  t(x[,roi_ant])%*%x[,roi_ant] )
beta_post_mean_ant= vcov_beta_ant %*%t(x[,roi_ant])%*%latent[i-1,]
transition_num=log(pbirthdeath_num) + log(pdj) +
  dmvnorm(x=beta_ant[roi_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 

#  
num= loglik_num +prior_num +transition_num
  
## denominator
#loglik
loglik_denom= dmvnorm(x=as.numeric( latent[i-1,]- x%*%beta_ant),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)

#prior
prior_denom = dmvnorm(x=beta_ant[roi_ant],mean=rep(0, length(roi_ant)),
                    sigma=diag(sigma_beta^2,length(roi_ant)), log=TRUE)

#transition
pbirthdeath_denom= pbirth_death(length(which(beta_ant[-1]!=0)))[[1]]
transition_denom=log(pbirthdeath_denom) + log(pbj) + dmvnorm(x=beta_new, mean = beta_post_mean, 
                                                     sigma = vcov_beta,log=TRUE) 

##  
denom= loglik_denom + prior_denom +transition_denom
  
### decide 
a= min(exp(num-denom),1)
rdiscrete(c(a,1-a))
}


####################################################
### acceptance probability of death 
acceptance_death_roi = function(x,y,beta_ant,beta_new,pbj,pdj){
  
  ## numerator 
  #loglik
  loglik_num= dmvnorm(x= as.numeric(latent[i-1,]- x%*%beta_new_full),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)
  
  #prior
  prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                      sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)
  
  
  #transition
  pbirthdeath_num=pbirth_death(length(which(beta_new[-1]!=0)))[[1]]
  
  
  roi_ant= which(beta_ant!=0)
  vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(roi_ant)) + 
                         t(x[,roi_ant])%*%x[,roi_ant] )
  beta_post_mean_ant= vcov_beta_ant %*%t(x[,roi_ant])%*%latent[i-1,]
  transition_num=log(pbirthdeath_num) + log(pbj) +
    dmvnorm(x=beta_ant[roi_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 
  
  #  
  num= loglik_num +prior_num +transition_num
  
  ## denominator
  #loglik
  loglik_denom= dmvnorm(x= as.numeric(latent[i-1,]- x%*%beta_ant),mean=rep(0,nr),
                        sigma =diag(nr),log = TRUE)
  
  #prior
  prior_denom = dmvnorm(x=beta_ant[roi_ant],mean=rep(0, length(roi_ant)),
                        sigma=diag(sigma_beta^2,length(roi_ant)), log=TRUE)
  
  #transition
  pbirthdeath_denom= pbirth_death(length(which(beta_ant[-1]!=0)))[[1]]
  transition_denom=log(pbirthdeath_denom) + log(pdj) + dmvnorm(x=beta_new, mean = beta_post_mean, 
                                                               sigma = vcov_beta,log=TRUE) 
  
  ##  
  denom= loglik_denom + prior_denom +transition_denom
  
  ### decide 
  a= min(exp(num-denom),1)
  rdiscrete(c(a,1-a))
}



#### INITIALIZE TO SELECT ROIs####
nc=ncol(x); p= nc-1
nr=n=nrow(x)
sigma_beta= 100 
iter=5000

beta= matrix(0,nrow=iter, ncol=nc) 
convergence_likelihood= numeric(iter)
convergence_likelihood[1]=0
latent= matrix(NA,nrow = iter,ncol = nr)
latent[1,]=replicate(3 ,ztr(x=x,y=y,beta=beta[1,]),F)[[3]] 
moves=numeric(iter)
moves[1]=0  

##convergence matrix to run two chains
nchains=1
convergence_likelihood_matrix = matrix(NA,iter,nchains)                                                                                                                                                                                                                                                                                                                                                                                       

##start
pb = txtProgressBar(min = 0, max = iter, initial = 0) 
for( j in 1:nchains){
beta[1,1]= rnorm(1,j,sd=2)
for(i in 2:iter){
  #browser()
  
  beta_length= length(which(beta[i-1,-1]!=0))
  # current roi 
  current_roi = setdiff(which(beta[i-1,]!=0),1)
  ## decide on birth  or death 
  decide=pbirth_death(beta_length)
  pbirthdeath= decide[[1]]
  move= decide[[2]]
  moves[i]=move
  if(move==1){
    ## birth
    ## which roi to include
    candidate_pbj= birth_roi(x,beta[i-1,],latent[i-1,])
    candidate=candidate_pbj[[2]]
    new_set_roi= sort(c(1,current_roi,candidate)) # 1 for beta_0
    pbj= candidate_pbj[[1]]
    
    ### generate beta
    vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set_roi)) + 
      t(x[,new_set_roi])%*%x[,new_set_roi] )
    beta_post_mean= vcov_beta %*%t(x[,new_set_roi])%*%latent[i-1,]
    beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
    
    ### compute pdj for the acceptance probability
    beta_new_full=numeric(nc);beta_new_full[new_set_roi]=beta_new
    pdj=death_roi(x=x,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
    
    
    ###test acceptance  
    alpha_acc_birth=acceptance_birth_roi(x=x, y=y,beta_ant=beta[i-1,], 
                                          beta_new=beta_new, pbj=pbj,pdj=pdj)
    if(alpha_acc_birth==1){
      beta[i,new_set_roi]= beta_new
      latent[i,]= replicate(3,ztr(x=x,y=y,beta=beta[i,]),F)[[3]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- x%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
    }else{ 
      ## Don't stay static move the chain 
    vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current_roi)+1) +
                           t(x[,c(1,current_roi)]) %*% x[,c(1,current_roi)])
    beta_static_mean= vcov_beta_static %*%t(x[,c(1,current_roi)])%*%latent[i-1,]
    beta_static= rmvnorm(10, mean = beta_static_mean, sigma = vcov_beta_static)[10,]
    beta[i,c(1,current_roi)]= beta_static
    latent[i,]= replicate(10,ztr(x=x,y=y,beta=beta[i,]),F)[[10]] 
    convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- x%*%beta[i,]),mean=rep(0,nr),
                                        sigma =diag(nr),log = TRUE)
    }
    
  }else{ 
    # death
    ## which roi to include
    candidate_pdj= death_roi(x,beta[i-1,],latent[i-1,])
    candidate=candidate_pdj[[2]]
    new_set_roi= c(1,setdiff(current_roi,candidate)) # 1 for beta_0
    pdj= candidate_pdj[[1]]
    
    ### generate beta
    vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set_roi)) + 
                        t(x[,new_set_roi])%*%x[,new_set_roi] )
    beta_post_mean= vcov_beta %*%t(x[,new_set_roi])%*%latent[i-1,]
    beta_new= rmvnorm(10,mean = beta_post_mean, sigma = vcov_beta)[10,]
    
    ### compute pbj for the acceptance probability
    beta_new_full=numeric(nc);beta_new_full[new_set_roi]=beta_new
    pbj=birth_roi(x=x,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
  
    
    
    ###test acceptance
    alpha_acc_death=acceptance_death_roi(x=x, y=y,beta_ant=beta[i-1,], 
                                         beta_new=beta_new, pbj=pbj,pdj=pdj)
    if(alpha_acc_death==1){
      beta[i,new_set_roi]= beta_new
      latent[i,]= replicate(10,ztr(x=x,y=y,beta=beta[i,]),F)[[10]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- x%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
    }else{ 
      ## Don't stay static move the chain 
      vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current_roi)+1) +
                                t(x[,c(1,current_roi)]) %*% x[,c(1,current_roi)])
      beta_static_mean= vcov_beta_static %*%t(x[,c(1,current_roi)])%*%latent[i-1,]
      beta_static= rmvnorm(10, mean = beta_static_mean, sigma = vcov_beta_static)[10,]
      beta[i,c(1,current_roi)]= beta_static
      latent[i,]=replicate(10, ztr(x=x,y=y,beta=beta[i,]),F )[[10]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- x%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
      }
    
    
  }
  setTxtProgressBar(pb,i)

}

convergence_likelihood_matrix[,j]= convergence_likelihood  
}

## log_posterior chains
plot(convergence_likelihood_matrix[,1],type='l',ylab='log posterior', xlab='iterations', col='blue')
lines(convergence_likelihood_matrix[,2], col='red')

### run at least 2 chains so  we can compute R_hat
rstan::Rhat(convergence_likelihood_matrix)

## model size
apply(beta[,-1],1, function(x){length(which(x[-1]!=0))}) %>% plot(ylab="P: number of ROIs", xlab="iterations")

##check acf
burn_in=floor(iter/3)
acf(convergence_likelihood_matrix[,1],50)

## select 1 sample from 10
comp=seq(burn_in,iter,10)

##posterior probability of model size
apply(beta,1, function(x){length(which(x[-1]!=0))}) %>% table() %>%  
  prop.table() %>% round(digits=4) 

#marginal ppi
plot(1:p,apply(beta[burn_in:iter,2:(p+1)], 2, function(x)mean(x!=0)),pch=3,xlab = 'ROIs',
 ylab="Marginal Posterior Probability of Inclusion",ylim=c(0,1));abline(h=0.5,col=2)

nao_nulo_estimado=apply(beta[burn_in:iter,1:(p+1)], 2, function(x)mean(x!=0))


### threshold for the posterior marginal posterior probability of inclusion(mppi)
tr= 0.5

selecionado= which(apply(beta[burn_in:iter,1:(p+1)], 2, function(x)mean(x!=0)>=tr))
selected_coefs = selecionado[-1] -1

## true versus selected 
true_beta= which(abs(beta_true)>0.5)

cbind(true_beta, selecionado)

##  posterior mean and 95% CI
tibble(selected_covariates = c('intercept', paste0('ROI ' , selecionado[-1]-1)),
       mppi= nao_nulo_estimado[selecionado],
       real = beta_true[true_beta],
post_mean = apply(beta[comp,selecionado],2,mean) ,
CI_95= do.call(rbind.data.frame,apply(beta[comp,selecionado],2,function(x)bayestestR::ci(x,method="HDI")) ))


## most visited models with posterior probability
models=apply(beta[burn_in:iter,], 1,
                     function(x){ which(x!=0) %>%as.character() %>% paste(collapse = ',')})

most_visited_models=tibble(models=table(models) %>% names(),
          mppi=table(models) %>% prop.table()%>% as.numeric() ) %>% 
  arrange(desc(mppi))

most_visited_models %>% head(5)

## bayesian model averaging for prediction 
latent_pred_matrix=matrix(NA,nrow = nrow(most_visited_models),ncol = length(ytest))

pb = txtProgressBar(min = 0, max = nrow(most_visited_models), initial = 0) 
for(i in 1:nrow(most_visited_models)){ 
  vars_model=strsplit(most_visited_models[i,]$models[1],',') %>% unlist() %>% as.numeric()
  beta_estimado=apply(beta[burn_in:iter,vars_model,drop=F], 2, mean)
  latent_pred_matrix[i,]=as.matrix(xtest[,vars_model],ncol=length(ytest))%*% beta_estimado

  setTxtProgressBar(pb,i)
}


## prediction
latent_pred_matrix = latent_pred_matrix*most_visited_models$mppi
latent_test = apply(latent_pred_matrix, 2, sum)
prob_test = pnorm(latent_test)
y_pred= ifelse(latent_test>0,1,0)
table(y_pred, ytest)

## confusion matrix 
mce= mean(y_pred!=ytest)
mce

## auc curve 
pred <- prediction( prob_test, ytest)
perf <- performance(pred,"tpr","fpr")
plot(perf)

## precision/recall curve (x-axis: recall, y-axis: precision)
perf1 <- performance(pred, "auc")

perf1@y.values[[1]] 


### lasso 
mod.lasso= cv.glmnet(x[,-1],y, family="binomial",type.measure = "auc",standardize=F)
plot(mod.lasso)
y_pred_lasso= predict(mod.lasso, newx = xtest[,-1],s = mod.lasso$lambda.1se,type='class') %>% as.numeric()
y_pred_lasso
table(y_pred_lasso, ytest)
mce_lasso= mean(y_pred_lasso!=ytest)
mce_lasso

pred <- prediction(y_pred_lasso %>% as.numeric(), ytest)
perf <- performance(pred,"tpr","fpr")
plot(perf)

## precision/recall curve (x-axis: recall, y-axis: precision)
perf1 <- performance(pred, "auc")
perf1@y.values[[1]]


