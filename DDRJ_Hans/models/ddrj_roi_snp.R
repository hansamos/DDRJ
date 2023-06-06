#############################################
#############################################
# DDRJ for joint ROIs and SNPs selection.  ##           
# Provided a design matrix with both       ##
#continous and categorical                 ##
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

### source to simulate or get real data 
source('../utils/get_rois_snps.R')

###################################
set.seed(123)
###################################

##data
p=200;n=300
data_sim = simulate_data_rois_snps(n=n, p=p, prop=0.85)

## let's do variable selection, so we use the full data
xz= data_sim$xz
y= data_sim$y


## for prediction using .85 train test 
#xz= data_sim$xztrain
#y= data_sim$ytrain
#xztest= data_sim$xztest
#ytest = data_sim$ytest

## true beta for comparison
beta_true= data_sim$beta_true 
sum(abs(beta_true)>0.5) -1   ## significant effect + intercept 


###################################
# Functions
###################################

### draw a value of a discrete distribution
#
rdiscrete<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val}

## choose to jump on model snp or roi
jump_roi_snp=function(nroi,nsnp){
  rbinom(1,1,nroi/(nsnp+nroi))} # 1 roi 0 snp  


### Decide on birth, death for snp
#p=: length current snp alpha and delta
# return probability and type of move
pbirth_death_snp = function(p){
  if(p>=0  & p<= 2*nsnp){ 
    move=ifelse(p==0,1,
                ifelse( p > 0 & p < 2*nsnp, rbinom(1,1,0.5)+1,2))
    pbirthdeath=ifelse(p==0 | p== 2*nsnp, 1, 0.5)
    list(pbirthdeath,move)}else{stop("p > col_number")}
  #move 1= birth, 2=death
}

### Decide on birth, death for rois
#k=: length current rois
pbirth_death_roi = function(p){
  if(p>=0  & p<= nroi){ 
    move=ifelse(p==0,1,
                ifelse( p > 0 & p < nroi, rbinom(1,1,0.5)+1,2))
    pbirthdeath=ifelse(p==0 | p==nroi, 1, 0.5)
    list(pbirthdeath,move)}else{stop("p > col_number")}
  #move 1= birth, 2=death
}



### compute birth pbj probability for candidate snps
#xz: matrix of roi-snps
#beta: current coefficients beta
#latent: current latent variable
birth_snp = function(xz,beta,latent){
  pb_size_beta=rep(NA,nc)
  cand= setdiff(which(beta[1:(1+nroi+nsnp)]==0),c(1:(nroi+1)))
  res= latent - (xz %*% beta) #compute current residual
  krusk=kruskalcpp(cbind(res,xz[,cand]))
  pbj_geral= krusk/sum(krusk)
  #pbj_geral = rep(1/length(cand), length(cand))# rj 
  pb_size_beta[cand]=pbj_geral
  snp_candidate_birth= rdiscrete(pbj_geral)
  insert_snp=cand[snp_candidate_birth]
  pbj= pb_size_beta[insert_snp]
  insert_snp=c(insert_snp, insert_snp+nsnp)
  
  list(pbj,insert_snp,pb_size_beta)
}


### compute birth pbj probability for candidate rois
#x: matrix of rois
#beta: current coefficients beta
birth_roi = function(xz,beta,latent){
  pb_size_beta=rep(NA,nc)
  cand= setdiff(which(beta[1:(1+nroi)]==0),1)
  res= latent - (xz %*% beta) #compute current residual
  corr=cor(res,xz[,cand]) %>%  abs()
  pbj_geral=corr/sum(corr)
  #pbj_geral = rep(1/length(cand), length(cand)) # rj
  pb_size_beta[cand]=pbj_geral
  roi_candidate_birth= rdiscrete(pbj_geral)
  insert_roi=cand[roi_candidate_birth]
  pbj= pb_size_beta[insert_roi]
  list(pbj,insert_roi,pb_size_beta)
}


### compute death pdj probability for candidate snps

#xz: matrix of rois-snps
#beta: current coefficients beta
#latent: current latent variable
death_snp = function(xz, beta, latent){
  
  ## number of snps candidate to be deleted
  pd_size_beta=rep(NA,nsnp)
  
  cand= setdiff(which(beta[1:(1+nroi+nsnp)]!=0),c(1:(nroi+1)))
  alpha=beta[cand]
  delta=beta[cand+nsnp]
  alpha_delta= 1/ ( abs(alpha) + abs(delta))
  
 pdj_geral= alpha_delta/sum(alpha_delta)
# pdj_geral = rep(1/length(cand), length(cand)) # rj with uniform jump
  pd_size_beta[cand]=pdj_geral
  snp_candidate_death = rdiscrete(pdj_geral)
  remove_snp= cand[snp_candidate_death]
  pdj= pd_size_beta[remove_snp]
  remove_snp = c(remove_snp,remove_snp+nsnp)
  list(pdj, remove_snp,pd_size_beta)
}



### compute death pdj probability for candidate rois
#xz: matrix of rois-snps
#beta: current coefficients beta
#latent: current latent variable
death_roi = function(xz, beta, latent){
  #browser()
  ## number of rois candidate to be deleted
  pd_size_beta=rep(NA,nc)
  cand= setdiff(which(beta[1:(1+nroi)]!=0),1)
  res= latent - xz%*% beta
  xi=numeric(length(cand))
  for (i in 1:length(cand)) {
    xi[i]= 1/abs(beta[cand[i]])
  }
  pdj_geral= xi/sum(xi)
  #pdj_geral=rep(1/length(cand), length(cand)) #rj with uniform jump
  pd_size_beta[cand]=pdj_geral
  roi_candidate_death = rdiscrete(pdj_geral)
  remove_roi= cand[roi_candidate_death]
  pdj= pd_size_beta[remove_roi]
  list(pdj, remove_roi,pd_size_beta)
}

### Latent variable
ztr=function(xz,y,beta){
  #beta=matrix(beta,ncol(xz),1,byrow = T)
  return(ifelse(y==1,rtnorm(nr,xz%*%beta,1,a=0),rtnorm(nr,xz%*%beta,1,b=0)))
}



### acceptance probability of birth 
acceptance_birth = function(xz,y,beta_ant,beta_new,pbj,pdj){
  
  ## numerator 
  #loglik
  loglik_num= dmvnorm(x= as.numeric(latent[i-1,]- xz%*%beta_new_full),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)
  
  #prior
  prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                      sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)
  
  
  #transition
  pbirthdeath_num=ifelse(jump==0, pbirth_death_snp(length(c(current_snp,candidate)))[[1]], 
                          pbirth_death_roi(length(c(current_roi,candidate)))[[1]])
  
  
  var_ant=current ## take everybody beta_0 included
  vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(var_ant)) + 
                         t(xz[,var_ant])%*%xz[,var_ant] )
  beta_post_mean_ant= vcov_beta_ant %*%t(xz[,var_ant])%*%latent[i-1,]
  transition_num=log(pbirthdeath_num) + log(pdj) +
    dmvnorm(x=beta_ant[var_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 
  
  #  
  num= loglik_num +prior_num +transition_num
  
  ## denominator
  #loglik
  loglik_denom= dmvnorm(x= as.numeric(latent[i-1,]- xz%*%beta_ant),mean=rep(0,nr),
                        sigma =diag(nr),log = TRUE)
  
  #prior
  prior_denom = dmvnorm(x=beta_ant[var_ant],mean=rep(0, length(var_ant)),
                        sigma=diag(sigma_beta^2,length(var_ant)), log=TRUE)
  
  #transition
  pbirthdeath_denom= ifelse(jump==0, pbirth_death_snp(length(current_snp))[[1]], 
                            pbirth_death_roi(length(current_roi))[[1]])
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
acceptance_death = function(xz,y,beta_ant,beta_new,pbj,pdj){
  
  ## numerator 
  #loglik
  loglik_num= dmvnorm(x=as.numeric( latent[i-1,]- xz%*%beta_new_full),mean=rep(0,nr),
                      sigma =diag(nr),log = TRUE)
  
  #prior
  prior_num = dmvnorm(x=beta_new,mean=rep(0, length(beta_new)),
                      sigma=diag(sigma_beta^2,length(beta_new)), log=TRUE)
  
  
  #transition
  pbirthdeath_num=ifelse(jump==0, pbirth_death_snp(length(setdiff(current_snp,candidate)))[[1]], 
                         pbirth_death_roi(length(setdiff(current_roi,candidate)))[[1]])
  
  
  var_ant= current
  vcov_beta_ant= solve(diag(1/(sigma_beta^2),length(var_ant)) + 
                         t(xz[,var_ant])%*%xz[,var_ant] )
  beta_post_mean_ant= vcov_beta_ant %*%t(xz[,var_ant])%*%latent[i-1,]
  transition_num=log(pbirthdeath_num) + log(pbj) +
    dmvnorm(x=beta_ant[var_ant], mean = beta_post_mean_ant, sigma = vcov_beta_ant,log=TRUE) 
  
  #  
  num= loglik_num +prior_num +transition_num
  
  ## denominator
  #loglik
  loglik_denom= dmvnorm(x= as.numeric(latent[i-1,]- xz[,var_ant]%*%beta_ant[var_ant]),mean=rep(0,nr),
                        sigma =diag(nr),log = TRUE)
  
  #prior
  prior_denom = dmvnorm(x=beta_ant[var_ant],mean=rep(0, length(var_ant)),
                        sigma=diag(sigma_beta^2,length(var_ant)), log=TRUE)
  
  #transition
  pbirthdeath_denom= ifelse(jump==0, pbirth_death_snp(length(current_snp))[[1]], 
                            pbirth_death_roi(length(current_roi))[[1]])
  transition_denom=log(pbirthdeath_denom) + log(pdj) + dmvnorm(x=beta_new, mean = beta_post_mean, 
                                                               sigma = vcov_beta,log=TRUE) 
  
  ##  
  denom= loglik_denom + prior_denom +transition_denom
  
  ### decide 
  a= min(exp(num-denom),1)
  rdiscrete(c(a,1-a))
}



#### INITIALIZE TO SELECT SNP and ROI ####
nc=ncol(xz)
nr=nrow(xz)
nsnp=p
nroi=p
iter=3000
burn_in=floor(iter/7)
sigma_beta=100
beta= matrix(0,nrow=iter, ncol=nc) 
beta[1,1]= rnorm(1,sd=2)
latent= matrix(NA,nrow = iter,ncol = nr)
latent[1,]=replicate(3,ztr(xz=xz,y=y,beta=beta[1,]),F)[[3]] 
convergence_likelihood= numeric(iter)
convergence_likelihood[1]=dmvnorm(x= as.numeric(latent[1,]- xz%*%beta[1,]),mean=rep(0,nr),
                                  sigma =diag(nr),log = TRUE)
rej_birth=rej_death=rep(NA,iter)
rej_birth[1]=rej_death[1]=2
moves=jumps=numeric(iter)
moves[1]=0
jumps[1]=0


##convergence matrix to run two chains
nchains=2
convergence_likelihood_matrix = matrix(NA,iter,nchains)                                                                                                                                                                                                                                                                                                                                                                                       


pb = txtProgressBar(min = 0, max = iter, initial = 0) 
for( j in 1:nchains){
  beta[1,1]= rnorm(1,j,sd=2)
for(i in 2:iter){
 jump=jump_roi_snp(nroi,nsnp)
 jumps[i]=jump
 #jump =1 roi , jump=0 snp
 if (jump==0) {
   #  browser()
   
   ## how many non zero snp
   beta_length_snp= length(which(beta[i-1,(nroi+2):nc]!=0))
   ##current roi 
   current_roi = setdiff(which(beta[i-1,1:(nroi+1)]!=0),1)
   ## current snp 
   current_snp = setdiff(which(beta[i-1,]!=0),1:(nroi+1))
   #current covs
   current= sort(c(1, current_roi,current_snp))
   ## decide on birth  or death 
   decide =pbirth_death_snp(beta_length_snp)
   pbirthdeath= decide[[1]]
   move= decide[[2]]
   moves[i]=move
   if(move==1){
     ## birth
     ## which snp to include
     candidate_pbj= birth_snp(xz,beta[i-1,],latent[i-1,])
     candidate=candidate_pbj[[2]]
     new_set= sort(c(1,current_roi,current_snp,candidate )) # 1 for beta_0
     pbj= candidate_pbj[[1]]
     
     ### generate beta
     vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set)) + 
                         t(xz[,new_set])%*%xz[,new_set] )
     beta_post_mean= vcov_beta %*%t(xz[,new_set])%*%latent[i-1,]
     beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
     
     ### compute pdj for the acceptance probability
     beta_new_full=numeric(nc);beta_new_full[new_set]=beta_new
     pdj=death_snp(xz=xz,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
     
     ###test acceptance  
     alpha_acc_birth=acceptance_birth(xz=xz, y=y,beta_ant=beta[i-1,], 
                                          beta_new=beta_new, pbj=pbj,pdj=pdj)
     rej_birth[i]=alpha_acc_birth
     if(alpha_acc_birth==1){
       beta[i,new_set]= beta_new
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }else{ 
       ## Don't stay static move the chain 
       vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current)) +
                                 t(xz[,current]) %*% xz[,current])
       beta_static_mean= vcov_beta_static %*%t(xz[,current])%*%latent[i-1,]
       beta_static= rmvnorm(3, mean = beta_static_mean, sigma = vcov_beta_static)[3,]
       beta[i,current]= beta_static
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]] 
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }
     
   }else{ 
     # death
     ## birth
     ## which snp to exclude
     candidate_pdj= death_snp(xz,beta[i-1,],latent[i-1,])
     candidate=candidate_pdj[[2]]
     new_set= sort(c(1,current_roi,setdiff(current_snp,candidate))) # 1 for beta_0
     pdj= candidate_pdj[[1]]
     
     ### generate beta
     vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set)) + 
                         t(xz[,new_set])%*%xz[,new_set] )
     beta_post_mean= vcov_beta %*%t(xz[,new_set])%*%latent[i-1,]
     beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
     
     ### compute pbj for the acceptance probability
     beta_new_full=numeric(nc);beta_new_full[new_set]=beta_new
     pbj=birth_snp(xz=xz,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
     
     ###test acceptance
     alpha_acc_death=acceptance_death(xz=xz, y=y,beta_ant=beta[i-1,], 
                                          beta_new=beta_new, pbj=pbj,pdj=pdj)
     rej_death[i]=alpha_acc_death
     if(alpha_acc_death==1){
       beta[i,new_set]= beta_new
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }else{ 
       ## Don't stay static move the chain 
       vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current)) +
                                 t(xz[,current]) %*% xz[,current])
       beta_static_mean= vcov_beta_static %*%t(xz[,current])%*%latent[i-1,]
       beta_static= rmvnorm(3, mean = beta_static_mean, sigma = vcov_beta_static)[3,]
       beta[i,current]= beta_static
       latent[i,]=replicate(3, ztr(xz=xz,y=y,beta=beta[i,]),F )[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }
     
     
   }
   
   
 }else{ ## jump into roi space
   
   ## how many non zero roi
   beta_length_roi= length(which(beta[i-1,2:(nroi+1)]!=0))
   ##current roi 
   current_roi = setdiff(which(beta[i-1,1:(nroi+1)]!=0),1)
   ## current snp 
   current_snp = setdiff(which(beta[i-1,]!=0),1:(nroi+1))
   #current covs
   current= sort(c(1, current_roi,current_snp))
   ## decide on birth  or death roi
   decide =pbirth_death_roi(beta_length_roi)
   pbirthdeath= decide[[1]]
   move= decide[[2]]
   moves[i]=move
   if(move==1){
     ## birth
     ## which roi to include
     candidate_pbj= birth_roi(xz,beta[i-1,],latent[i-1,])
     candidate=candidate_pbj[[2]]
     new_set= sort(c(1,current_roi,current_snp,candidate )) # 1 for beta_0
     pbj= candidate_pbj[[1]]
     
     ### generate beta
     vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set)) + 
                         t(xz[,new_set])%*%xz[,new_set] )
     beta_post_mean= vcov_beta %*%t(xz[,new_set])%*%latent[i-1,]
     beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
     
     ### compute pdj for the acceptance probability
     beta_new_full=numeric(nc);beta_new_full[new_set]=beta_new
     pdj=death_roi(xz=xz,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
     
     ###test acceptance  
     alpha_acc_birth=acceptance_birth(xz=xz, y=y,beta_ant=beta[i-1,], 
                                      beta_new=beta_new, pbj=pbj,pdj=pdj)
     rej_birth[i]=alpha_acc_birth
     if(alpha_acc_birth==1){
       beta[i,new_set]= beta_new
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }else{ 
       ## Don't stay static move the chain 
       vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current)) +
                                 t(xz[,current]) %*% xz[,current])
       beta_static_mean= vcov_beta_static %*%t(xz[,current])%*%latent[i-1,]
       beta_static= rmvnorm(3, mean = beta_static_mean, sigma = vcov_beta_static)[3,]
       beta[i,current]= beta_static
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]] 
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }
     
   }else{ 
     # death
     
     ## which roi to exclude
     candidate_pdj= death_roi(xz,beta[i-1,],latent[i-1,])
     candidate=candidate_pdj[[2]]
     new_set= sort(c(1,setdiff(current_roi,candidate),current_snp)) # 1 for beta_0
     pdj= candidate_pdj[[1]]
     
     ### generate beta
     vcov_beta= solve( diag(1/(sigma_beta^2),length(new_set)) + 
                         t(xz[,new_set])%*%xz[,new_set] )
     beta_post_mean= vcov_beta %*%t(xz[,new_set])%*%latent[i-1,]
     beta_new= rmvnorm(3,mean = beta_post_mean, sigma = vcov_beta)[3,]
     
     ### compute pbj for the acceptance probability
     beta_new_full=numeric(nc);beta_new_full[new_set]=beta_new
     pbj=birth_roi(xz=xz,beta=beta_new_full,latent = latent[i-1,])[[3]][candidate[1]]
     
     ###test acceptance
     alpha_acc_death=acceptance_death(xz=xz, y=y,beta_ant=beta[i-1,], 
                                      beta_new=beta_new, pbj=pbj,pdj=pdj)
     rej_death[i]=alpha_acc_death
     if(alpha_acc_death==1){
       beta[i,new_set]= beta_new
       latent[i,]= replicate(3,ztr(xz=xz,y=y,beta=beta[i,]),F)[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }else{ 
       ## Don't stay static move the chain 
       vcov_beta_static= solve(diag(1/(sigma_beta^2),length(current)) +
                                 t(xz[,current]) %*% xz[,current])
       beta_static_mean= vcov_beta_static %*%t(xz[,current])%*%latent[i-1,]
       beta_static= rmvnorm(3, mean = beta_static_mean, sigma = vcov_beta_static)[3,]
       beta[i,current]= beta_static
       latent[i,]=replicate(3, ztr(xz=xz,y=y,beta=beta[i,]),F )[[3]]
       convergence_likelihood[i] = dmvnorm(x= as.numeric(latent[i,]- xz%*%beta[i,]),mean=rep(0,nr),
                                           sigma =diag(nr),log = TRUE)
     }
     
     
   }
   
 }#end else jump 
 setTxtProgressBar(pb,i)
}#end for
  convergence_likelihood_matrix[,j]= convergence_likelihood  
}  

## log_posterior chains
plot(convergence_likelihood_matrix[,1],type='l',ylab='log posterior', xlab='iterations', col='blue')
lines(convergence_likelihood_matrix[,2], col='red')

### run at least 2 chains so  we can compute R_hat
rstan::Rhat(convergence_likelihood_matrix)

## model size
apply(beta,1, function(x){length(which(x[-1]!=0))}) %>% plot()

##posterior probability of model size
apply(beta,1, function(x){length(which(x[-1]!=0))}) %>% table() %>%  
  prop.table() %>% round(digits=4) 

### selected features
comp= seq(burn_in,iter,10)
par(mfrow=c(1,2))
rois_space = setdiff(which(jumps!=0) , 1:burn_in)
snps_space = setdiff(which(jumps==0),1:burn_in)
plot(1:nroi,apply(beta[rois_space,2:(nroi+1)], 2, function(x)mean(x!=0)),pch=3,xlab = 'ROIs',
     ylab="Posterior Probability of Inclusion",ylim=c(0,1),main = "full dataset")

abline(h=0.5,col=2)
plot(1:nsnp,apply(beta[snps_space,(nroi+2):(ncol(beta)-nsnp)], 2, function(x)mean(x!=0)),pch=3,xlab = 'SNPs',
     ylab="Posterior Probabilty of Inclusion",ylim=c(0,1),main="full dataset")
abline(h=0.5,col=2)

#rois
tr = 0.5
which(apply(beta[burn_in:iter,2:(nroi+1)], 2, function(x)mean(x!=0))>tr)
#snps
which(apply(beta[burn_in:iter,(nroi+1):(ncol(beta)-nsnp)], 2, function(x)mean(x!=0))>=tr)

nao_nulo_estimado = apply(beta[burn_in:iter,1:nc], 2, function(x)mean(x!=0))
selecionado=which(apply(beta[burn_in:iter,1:nc], 2, function(x)mean(x!=0))>=tr)
selecionado
### true beta 
true_beta= which(abs(beta_true)>0.5)
true_beta

cbind(selecionado, true_beta)


### posterior mean  and 95 CI
tibble(selected_covariates = c('intercept', paste0('feature ' , selecionado[-1]-1)),
       mppi= nao_nulo_estimado[selecionado],
       true_coef = beta_true[selecionado],
       post_mean = apply(beta[comp,selecionado],2,mean) ,
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
  latent_pred_matrix[i,]=as.matrix(xztest[,vars_model],ncol=length(ytest))%*% beta_estimado
  
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


#### test data
est=apply(beta[comp,selecionado], 2, mean)
latent_est = xz[,selecionado]%*%est
latent_real=xz %*% beta_true[pre_vars] +rnorm(nrow(xz))
plot(latent_real,latent_est)
cor(latent_est,latent_real)
ypred=ifelse(latent_est>0,1,0)
table(y, ypred)
mean(ypred!=y)

#save.image("new_sim_roi_snp/sim/ddrj_300.RData")

for (i in 1:ncol(beta)) {
  browser()
  print(i)
  #par(mfrow=c(1,2))
  plot(beta[,i],type='l',main=paste0("beta",sep="_",i-1))#;abline(h=beta_true[i],col=2)
  #plot(beta[,i] %>% density); abline(v=beta_true[i],col=2)
}

##density
for (i in 1:ncol(beta)) {
  browser()
  print(i)
  #par(mfrow=c(1,2))
  #plot(beta[,i],type='l');abline(h=beta_true[i],col=2)
  plot(beta[,i] %>% density,main=paste0("beta",sep="_",i-1))#; abline(v=beta_true[i],col=2)
}


for (i in 1:nrow(beta)) {
  browser()
  plot(beta[i,],xlab='roi and snp', ylab= 'coef');abline(h=0,col=2);abline(v=117,col=3)
  
  #par(mfrow=c(1,2))
  #plot(beta[,i],type='l');abline(h=beta_true[i],col=2)
  #plot(beta[,i] %>% density,main=paste0("beta",sep="_",i-1)); abline(v=beta_true[i],col=2)
}
#rm(latent)
#load("quali/10/train1.RData")

rm(list=ls())

