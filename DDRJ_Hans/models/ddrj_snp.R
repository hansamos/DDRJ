#############################################
#############################################
# DDRJ for SNPs selection.                 ##
# Provided a  matrix with categorical      ##
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
library(Rcpp)

#source c++ code for Kruskal wallis
sourceCpp("../utils/kruskal.cpp")

##source to simulate ou get real data
source("../utils/get_snps.R")

###################################
set.seed(123)
###################################


##simulated data 
data_sim = simulate_data_snp(n=300, p=300, prop=0.85)

## let's do variable selection, so we use the full data
z= data_sim$z
y= data_sim$y

## for prediction using .85 train test 
z= data_sim$ztrain
y= data_sim$ytrain
ztest= data_sim$ztest
ytest = data_sim$ytest

## real beta for comparison
beta_true= data_sim$beta_true 
1+ (sum(abs(beta_true)>0.5) -1)/2  ## significant effect + intercept 



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


### Decide on birth, death or 
#p=: length(current beta )-1    ( -1 to exclude beta_0 )
# return probability and type of move
pbirth_death = function(p){
  if(p>=0  & p<= 2*nsnp){ 
    move=ifelse(p==0,1,
                ifelse( p > 0 & p < 2*nsnp, rbinom(1,1,0.5)+1,2))
    pbirthdeath=ifelse(p==0 | p== 2*nsnp, 1, 0.5)
    list(pbirthdeath,move)}else{stop("p > col_number")}
  #move 1= birth, 2=death
}



### compute birth pbj probability for candidate snps
#z: matrix of snps
#beta: current coefficients beta
#latent: current latent variable
birth_snp = function(z,beta,latent){
  pb_size_beta=rep(NA,nsnp)
  cand= setdiff(which(beta[1:(nsnp+1)]==0),1)
  res= latent - (z %*% beta) #compute current residual
  krusk=kruskalcpp(cbind(res,z[,cand]))
  pbj_geral= krusk/sum(krusk)
  #pbj_geral= rep(1/length(cand),length(cand)) #rj with uniform move
  pb_size_beta[cand]=pbj_geral
  snp_candidate_birth= rdiscrete(pbj_geral)
  insert_snp=cand[snp_candidate_birth]
  pbj= pb_size_beta[insert_snp]
  insert_snp=c(insert_snp, insert_snp+nsnp)
  
  list(pbj,insert_snp,pb_size_beta)
}

### compute death pdj probability for candidate snps

#z: matrix of snps
#beta: current coefficients beta
#latent: current latent variable
death_snp = function(z, beta, latent){
  
  ## number of snps candidate to be deleted
  pd_size_beta=rep(NA,nsnp)
  cand=setdiff(which(beta[1:(nsnp+1)]!=0),1)
  alpha=beta[cand]
  delta=beta[cand+nsnp]
  alpha_delta= 1/ ( abs(alpha) + abs(delta))
 
  pdj_geral= alpha_delta/sum(alpha_delta)
  #pdj_geral= rep(1/length(cand),length(cand)) #rj with uniform move
  pd_size_beta[cand]=pdj_geral
  snp_candidate_death = rdiscrete(pdj_geral)
  remove_snp= cand[snp_candidate_death]
  pdj= pd_size_beta[remove_snp]
  remove_snp = c(remove_snp,remove_snp+nsnp)
  list(pdj, remove_snp,pd_size_beta)
}

### Latent variable
ztr=function(z,y,beta){
  beta=matrix(beta,ncol(z),1,byrow = T)
  return(ifelse(y==1,rtnorm(nr,z%*%beta,1,a=0),rtnorm(nr,z%*%beta,1,b=0)))
}



### acceptance probability of birth 
acceptance_birth_snp = function(z,y,beta_ant,beta_new,pbj,pdj){
  
  ## numerator 
  #loglik
  loglik_num= dmvnorm(x= as.numeric(latent[i-1,]- z%*%beta_new_full),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)
  
  #prior
  prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                      sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)
  
  
  #transition
  pbirthdeath_num=pbirth_death(length(which(beta_new[-1]!=0)))[[1]]
  
  
  snp_ant= which(beta_ant!=0) ## take everybody beta_0 included
  vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(snp_ant)) + 
                         t(z[,snp_ant])%*%z[,snp_ant] )
  beta_post_mean_ant= vcov_beta_ant %*%t(z[,snp_ant])%*%latent[i-1,]
  transition_num=log(pbirthdeath_num) + log(pdj) +
    dmvnorm(x=beta_ant[snp_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 
  
  #  
  num= loglik_num +prior_num +transition_num
  
  ## denominator
  #loglik
  loglik_denom= dmvnorm(x= as.numeric(latent[i-1,]- z%*%beta_ant),mean=rep(0,nr),
                        sigma =diag(nr),log = TRUE)
  
  #prior
  prior_denom = dmvnorm(x=beta_ant[snp_ant],mean=rep(0, length(snp_ant)),
                        sigma=diag(sigma_beta^2,length(snp_ant)), log=TRUE)
  
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
acceptance_death_snp = function(z,y,beta_ant,beta_new,pbj,pdj){
  
  ## numerator 
  #loglik
  loglik_num= dmvnorm(x=as.numeric( latent[i-1,]- z%*%beta_new_full),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)
  
  #prior
  prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                      sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)
  
  
  #transition
  pbirthdeath_num=pbirth_death(length(which(beta_new[-1]!=0)))[[1]]
  
  
  snp_ant= which(beta_ant!=0)
  vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(snp_ant)) + 
                         t(z[,snp_ant])%*%z[,snp_ant] )
  beta_post_mean_ant= vcov_beta_ant %*%t(z[,snp_ant])%*%latent[i-1,]
  transition_num=log(pbirthdeath_num) + log(pbj) +
    dmvnorm(x=beta_ant[snp_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 
  
  #  
  num= loglik_num +prior_num +transition_num
  
  ## denominator
  #loglik
  loglik_denom= dmvnorm(x= as.numeric(latent[i-1,]- z[,snp_ant]%*%beta_ant[snp_ant]),mean=rep(0,nr),
                        sigma =diag(nr),log = TRUE)
  
  #prior
  prior_denom = dmvnorm(x=beta_ant[snp_ant],mean=rep(0, length(snp_ant)),
                        sigma=diag(sigma_beta^2,length(snp_ant)), log=TRUE)
  
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

#### INITIALIZE TO SELECT SNP ####
nc=ncol(z)
nr=nrow(z)
nsnp=(nc-1)/2
sigma_beta=100


#### DDRJ ####
iter=5000
burn_in=floor(iter/5)
comp= seq(burn_in,iter,10)
beta= matrix(0,nrow=iter, ncol=nc) 
beta[1,1]= rnorm(1,2)
latent= matrix(NA,nrow = iter,ncol = nr)
latent[1,]=replicate(3,ztr(z=z,y=y,beta=beta[1,]),F)[[3]] 
convergence_likelihood= numeric(iter)
convergence_likelihood[1]=dmvnorm(x= as.numeric(latent[1,]- z%*%beta[1,]),mean=rep(0,nr),
                                  sigma =diag(nr),log = TRUE)
moves=numeric(iter)
moves[1]=0


##convergence matrix to run two chains
nchains=1
convergence_likelihood_matrix = matrix(NA,iter,nchains)                                                                                                                                                                                                                                                                                                                                                                                       

pb = txtProgressBar(min = 0, max = iter, initial = 0) 
for( j in 1:nchains){
  beta[1,1]= rnorm(1,j,sd=2)
for(i in 2:iter){

#  browser()
  beta_length= length(which(beta[i-1,2:nc]!=0))
  ## current snp 
  current_snp = setdiff(which(beta[i-1,]!=0),1)
  ## decide on birth  or death 
  decide =pbirth_death(beta_length)
  pbirthdeath= decide[[1]]
  move= decide[[2]]
  moves[i]=move
  if(move==1){
    ## birth
    ## which snp to include
    candidate_pbj= birth_snp(z,beta[i-1,],latent[i-1,])
    candidate=candidate_pbj[[2]]
    new_set_snp= sort(c(1,current_snp,candidate)) # 1 for beta_0
    pbj= candidate_pbj[[1]]
    
    ### generate beta
    vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set_snp)) + 
                        t(z[,new_set_snp])%*%z[,new_set_snp] )
    beta_post_mean= vcov_beta %*%t(z[,new_set_snp])%*%latent[i-1,]
    beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
    
    ### compute pdj for the acceptance probability
    beta_new_full=numeric(nc);beta_new_full[new_set_snp]=beta_new
    pdj=death_snp(z=z,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
    
    ###test acceptance  
    alpha_acc_birth=acceptance_birth_snp(z=z, y=y,beta_ant=beta[i-1,], 
                                         beta_new=beta_new, pbj=pbj,pdj=pdj)
    if(alpha_acc_birth==1){
      beta[i,new_set_snp]= beta_new
      latent[i,]= replicate(3,ztr(z=z,y=y,beta=beta[i,]),F)[[3]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- z%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
    }else{ 
      ## Don't stay static move the chain 
      vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current_snp)+1) +
                                t(z[,c(1,current_snp)]) %*% z[,c(1,current_snp)])
      beta_static_mean= vcov_beta_static %*%t(z[,c(1,current_snp)])%*%latent[i-1,]
      beta_static= rmvnorm(3, mean = beta_static_mean, sigma = vcov_beta_static)[3,]
      beta[i,c(1,current_snp)]= beta_static
      latent[i,]= replicate(3,ztr(z=z,y=y,beta=beta[i,]),F)[[3]] 
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- z%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
    }
    
  }else{ 
    # death
    ## birth
    ## which snp to include
    candidate_pdj= death_snp(z,beta[i-1,],latent[i-1,])
    candidate=candidate_pdj[[2]]
    new_set_snp= sort(c(1,setdiff(current_snp,candidate))) # 1 for beta_0
    pdj= candidate_pdj[[1]]
    
    ### generate beta
    vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set_snp)) + 
                        t(z[,new_set_snp])%*%z[,new_set_snp] )
    beta_post_mean= vcov_beta %*%t(z[,new_set_snp])%*%latent[i-1,]
    beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
    
    ### compute pbj for the acceptance probability
    beta_new_full=numeric(nc);beta_new_full[new_set_snp]=beta_new
    pbj=birth_snp(z=z,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
    
    ###test acceptance
    alpha_acc_death=acceptance_death_snp(z=z, y=y,beta_ant=beta[i-1,], 
                                         beta_new=beta_new, pbj=pbj,pdj=pdj)
    if(alpha_acc_death==1){
      beta[i,new_set_snp]= beta_new
      latent[i,]= replicate(3,ztr(z=z,y=y,beta=beta[i,]),F)[[3]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- z%*%beta[i,]),mean=rep(0,nr),
                                          sigma =diag(nr),log = TRUE)
    }else{ 
      ## Don't stay static move the chain 
      vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current_snp)+1) +
                                t(z[,c(1,current_snp)]) %*% z[,c(1,current_snp)])
      beta_static_mean= vcov_beta_static %*%t(z[,c(1,current_snp)])%*%latent[i-1,]
      beta_static= rmvnorm(10, mean = beta_static_mean, sigma = vcov_beta_static)[10,]
      beta[i,c(1,current_snp)]= beta_static
      latent[i,]=replicate(10, ztr(z=z,y=y,beta=beta[i,]),F )[[10]]
      convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- z%*%beta[i,]),mean=rep(0,nr),
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
apply(beta,1, function(x){length(which(x[-1]!=0))/2}) %>% 
  plot(ylab="K: number of SNPs", xlab="iterations")

##check acf
burn_in=floor(iter/5)
acf(convergence_likelihood,50)

## select 1 sample from 10
comp=seq(burn_in,iter,10)

##posterior probability of model size
apply(beta,1, function(x){length(which(x[-1]!=0))/2}) %>% table() %>%  
  prop.table() %>% round(digits=4) 

#marginal ppi
p=nsnp
plot(1:p,apply(beta[burn_in:iter,2:(p+1)], 2, function(x)mean(x!=0)),pch=3,xlab = 'SNPs',
     ylab="Marginal Posterior Probability of Inclusion",ylim=c(0,1));abline(h=0.5,col=2)

nao_nulo_estimado=apply(beta[burn_in:iter,1:(p+1)], 2, function(x)mean(x!=0))


### threshold for the posterior marginal posterior probability of inclusion(mppi)
tr= 0.5

selecionado= which(apply(beta[burn_in:iter,1:(2*nsnp+1)], 2, function(x)mean(x!=0)>=tr))


## true versus selected 
true_beta= which(abs(beta_true)>0.5)

cbind(true_beta, selecionado)


##  posterior mean and 95% CI
tibble(selected_covariates = c('intercept', paste0('SNP ' , selecionado[-1]-1)),
       mppi= nao_nulo_estimado[selecionado],
       true_coef = beta_true[true_beta],
       post_mean_alpha = apply(beta[comp,selecionado],2,mean) ,
       CI_95= do.call(rbind.data.frame,
                      apply(beta[comp,selecionado],2,function(x)bayestestR::ci(x,method="HDI")) ),
       
     
       )


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
  latent_pred_matrix[i,]=as.matrix(ztest[,vars_model],ncol=length(ytest))%*% beta_estimado
  
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
mod.lasso= cv.glmnet(z[,2:(p+1)],y, family="binomial",type.measure = "auc",standardize=F)
plot(mod.lasso)
y_pred_lasso= predict(mod.lasso, newx = ztest[,2:(p+1)],s = mod.lasso$lambda.1se,type='class') %>% as.numeric()
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


