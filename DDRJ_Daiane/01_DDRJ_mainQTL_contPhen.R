###################################
###################################
# This is an R code of the
# DDRJ procedure for identifying main QTLs when the available genetic map is not dense enough
# and the method try to identify QTLs and predict their genotype between the available markers.
# It assumes a continuous phenotype as response variable.
# Assumptions: 1 - only one QTL in a marker interval
#              2 - independent individuals from a F2 population
#              3 - only additive and dominance effect of QTLs influence the phenotype
#              4 - Uniform prior distribution for the number of QTLs, Uniform prior distribution for QTLs'
# location, Normal prior distribution for the general mean and QTLs' effects and Gamma prior distribution for
# precision (1/error variance) parameter, as defined in the main paper published at Genetics(2016).
#
###################################
###################################
#
###################################
# Functions
###################################
#
### Cramér association metric
#
cv.test = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(as.numeric(CV))}
#
################
### sample a discrete value from the set {1,...,K}, where p is the probability vector of each value
#
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 return(val)}
#
################
### compute the recombination rate through Haldane function for a specific genetic distance (dist) in Morgan (M) metric
#
Haldane<-function(dist){
 recomb<-0.5*(1-exp(-2*dist))
 return(recomb)}
#
################
### compute the probability of selecting a specific marker using Kruskal-wallis and sample a marker
#
prob.selec.marc<-function(dados, res1, qtls, loc.marc){
  if (length(qtls)==0) krusk<-kruskalcpp(dados) else krusk<-kruskalcpp(cbind(res1,dados[,2:ncol(dados)]))
  if (sum(qtls<=loc.marc[2])>0) krusk[1]<-0
  for (i in 2:(length(loc.marc)-1)) if ((sum(qtls>=loc.marc[i-1] & qtls<=loc.marc[i])>0)&(sum(qtls>=loc.marc[i] & qtls<=loc.marc[i+1])>0)) krusk[i]<-0
  if (sum(qtls>=loc.marc[length(loc.marc)-1])>0) krusk[length(loc.marc)]<-0
  prob<-krusk/sum(krusk)
  marc<-rDiscreta(prob)
  list(marc,log(prob[marc]),krusk)}
#
################
### compute the probability of excluding a specific marker from the current model and sample a marker to be excluded
#
prob.excl.marc<-function(num.QTLs,vet.coef){
 efeitos<-NULL
 for (i in 1:num.QTLs) efeitos[i]<-1/sum(c(abs(vet.coef[(2*i),1]),abs(vet.coef[((2*i)+1),1])))
 proba<-efeitos/sum(efeitos)
 qtl<-rDiscreta(proba)
 list(qtl,log(proba[qtl]))}
#
################
### sample the position of the new QTL around the chosen marker
#
pos.qtl<-function(qtls,loc.marc,marc,krusk){
 if (marc==1 | sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])>0){
  marc.ini<-marc
  marc.fim<-marc+1} else {
  if (marc==length(loc.marc) | sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])>0){
   marc.ini<-marc-1
   marc.fim<-marc} else {
   if ((sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])==0)&(sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])==0)){
    marc.ini<-marc-1
    marc.fim<-marc+1}}}
 marc.pos<-(loc.marc[marc.ini:marc.fim]-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 estas<-krusk[marc.ini:marc.fim]
 mi.pos<-(t(marc.pos)%*%estas)/sum(estas)
 parA<-mi.pos/(1-mi.pos)
 parB<-1
 uger<-rbeta(1,parA,parB)
 loc.qtl<-loc.marc[marc.ini]+(loc.marc[marc.fim]-loc.marc[marc.ini])*uger
 dens<-dbeta(uger,parA,parB,log=TRUE)
 list(loc.qtl,dens)}
#
dens.pos.qtl<-function(qtls,loc.marc,marc,krusk,QTLmg){
 if (marc==1 | sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])>0){
  marc.ini<-marc
  marc.fim<-marc+1} else {
  if (marc==length(loc.marc) | sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])>0){
   marc.ini<-marc-1
   marc.fim<-marc} else {
   if ((sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])==0)&(sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])==0)){
    marc.ini<-marc-1
    marc.fim<-marc+1}}}
 marc.pos<-(loc.marc[marc.ini:marc.fim]-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 estas<-krusk[marc.ini:marc.fim]
 mi.pos<-(t(marc.pos)%*%estas)/sum(estas)
 parA<-mi.pos/(1-mi.pos)
 parB<-1
 uger<-(QTLmg-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 dens<-dbeta(uger,parA,parB)
 return(dens)}
#
################
### compute the prior density of QTLs' location
#
priori.loc.qtls<-function(num.qtls,loc.marc){
 loc.qtls<-NULL
 num.marc<-length(loc.marc)
 prob.loc<-0
 if (num.qtls>0){
 for (i in 1:num.qtls){
   loc.qtls[i]<-runif(1,min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)])
   prob.loc<-prob.loc+dunif(loc.qtls[i],min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)],log=TRUE)}}
 list(loc.qtls,prob.loc)}
#
dens.priori.loc.qtls<-function(loc.qtls,loc.marc){
 num.qtls<-length(loc.qtls)
 num.marc<-length(loc.marc)
 prob.loc<-0
 if (num.qtls>0){
 for (i in 1:num.qtls){
   prob.loc<-prob.loc+dunif(loc.qtls[i],min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)],log=TRUE)}}
 return(prob.loc)}
#
################
### compute the probability of QTL's genotype based on the flanking markers' genotype
#
# gen = 1 if dominant homozygous
# gen = 0 if heterozygous
# gen = -1 if recessive homozygous
gen.igual<-function(recomb,gen1) if (gen1==0) {(recomb*recomb)+(1-recomb)**2} else {(1-recomb)**2}
gen.dif1<-function(recomb) 2*(1-recomb)*recomb
gen.dif2<-function(recomb) recomb*recomb
#
calc.prob.gen<-function(gen1,gen2,recomb){
 if (abs(gen1-gen2)==0) {prob<-gen.igual(recomb,gen1)} else {
  if (abs(gen1-gen2)==1) {if (gen1==0) {prob<-gen.dif1(recomb)/2} else {prob<-gen.dif1(recomb)}}
  else {prob<-gen.dif2(recomb)}}
 return(prob)}
#
################
### sample from the posterior distribution of the error variance
#
poster.sigma2<-function(neta.a,neta.b,residuos){
 alpha<-(length(residuos)/2)+neta.a
 beta<-(sum(residuos^2)/2)+neta.b
 sigma2<-1/(rgamma(1,alpha,beta))
 dens<-dgamma((1/sigma2),alpha,beta,log=TRUE)
 list(sigma2,dens)}
#
################
### sample from the posterior distribution of the general mean
#
poster.mi<-function(media.mi,sigma2.mi,sigma2,residuos,mi.anterior){
 res.mi<-residuos+mi.anterior
 quo<-(length(residuos)/sigma2)+(1/sigma2.mi)
 media<-((sum(res.mi)/sigma2)+(media.mi/sigma2.mi))/quo
 variancia<-1/quo
 mi<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(mi,media,sqrt(variancia),log=TRUE)
 list(mi,dens)}
#
################
### sample from the posterior distribution of the additive effect of each QTL
#
poster.alpha<-function(media.alpha,sigma2.alpha,sigma2,residuos,alpha.anterior,gen.QTL){
 res.alpha<-residuos+(alpha.anterior*gen.QTL)
 quo<-(sum(gen.QTL^2)/sigma2)+(1/sigma2.alpha)
 media<-((sum(gen.QTL*res.alpha)/sigma2)+(media.alpha/sigma2.alpha))/quo
 variancia<-1/quo
 alpha<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(alpha,media,sqrt(variancia),log=TRUE)
 list(alpha,dens)}
#
################
### sample from the posterior distribution of the dominance effect of each QTL
#
poster.delta<-function(media.delta,sigma2.delta,sigma2,residuos,delta.anterior,gen.QTL.dom){
 res.delta<-residuos+(delta.anterior*gen.QTL.dom)
 quo<-(sum(gen.QTL.dom^2)/sigma2)+(1/sigma2.delta)
 media<-((sum(gen.QTL.dom*res.delta)/sigma2)+(media.delta/sigma2.delta))/quo
 variancia<-1/quo
 delta<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(delta,media,sqrt(variancia),log=TRUE)
 list(delta,dens)}
#
################
### choose between the movement of death or birth of a QTL
#
dec.sp.mg<-function(num.QTLs,num.marc){
 if (num.QTLs==0) {psplit<-1; pmerge<-0} else {if (num.QTLs==(num.marc-1)) {psplit<-0;pmerge<-1} else {psplit<-pmerge<-1/2}}
 prob<-c(psplit,pmerge)
 ind.sp.mg<-rDiscreta(prob)
 list(ind.sp.mg,log(prob))}
#
#################
### sample a candidate of birth and compute its transition function
#
gera.inclusao.QTL<-function(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,alpha.vig,delta.vig,media.mi,sigma2.mi,media.alpha,sigma2.alpha,media.delta,sigma2.delta,neta.a,neta.b){
 marcador<-prob.selec.marc(dados,residuos,pos.qtls,loc.marc)
 qtl<-pos.qtl(pos.qtls,loc.marc,marcador[[1]],marcador[[3]])
 #
 dQTL<-qtl[[1]]
 Marc1<-sum(loc.marc<=dQTL)
 Marc2<-length(loc.marc)-(sum(loc.marc>=dQTL)-1)
 dM1<-loc.marc[Marc1]
 dM2<-loc.marc[Marc2]
 r12<-Haldane(abs(dM2-dM1))
 r1<-Haldane(abs(dQTL-dM1))
 r2<-Haldane(abs(dQTL-dM2))
 #
 matriz.prob.gen.QTL<-NULL
 probQTL<-numeric(3)
 for (j in 1:nrow(dados)){
  gen<-c(-1,0,1)
  for (i in 1:3) probQTL[i]<-(calc.prob.gen(dados[j,(Marc1+1)],gen[i],r1)*calc.prob.gen(gen[i],dados[j,(Marc2+1)],r2))/calc.prob.gen(dados[j,(Marc1+1)],dados[j,(Marc2+1)],r12)
 matriz.prob.gen.QTL<-rbind(matriz.prob.gen.QTL,probQTL)}
 #
 dados.QTL<-matrix(0,nrow(dados),1)
 for (i in 1:nrow(dados)) dados.QTL[i,1]<-gen[rDiscreta(matriz.prob.gen.QTL[i,])]
 mat.delinea<-cbind(mat.delinea,dados.QTL)
 #
 alpha<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,alpha.vig,dados.QTL[,1])
 vet.coef<-rbind(vet.coef,alpha[[1]])
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 dados.QTL<-cbind(dados.QTL,(1-abs(dados.QTL[,1])))
 delta<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,delta.vig,dados.QTL[,2])
 mat.delinea<-cbind(mat.delinea,dados.QTL[,2])
 vet.coef<-rbind(vet.coef,delta[[1]])
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])
 vet.coef[1,1]<-mi[[1]]
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
 list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],marcador[[2]],qtl[[2]],matriz.prob.gen.QTL,alpha[[2]],delta[[2]],mi[[2]],sigma2[[2]],c(pos.qtls,qtl[[1]]))}
#
#################
### sample a candidate of death and compute its transition function
#
gera.exclusao.QTL<-function(pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.mi,sigma2.mi,neta.a,neta.b){
 num.QTLs<-length(pos.qtls)
 qtl<-prob.excl.marc(num.QTLs,vet.coef)
 #
 pos.qtls<-pos.qtls[-qtl[[1]]]
 mat.delinea<-matrix(mat.delinea[,-c((2*qtl[[1]]),(2*qtl[[1]]+1))],nrow=nrow(dados))
 vet.coef<-matrix(vet.coef[-c((2*qtl[[1]]),(2*qtl[[1]]+1))],ncol=1)
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])
 vet.coef[1,1]<-mi[[1]]
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
 list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],qtl[[2]],mi[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],pos.qtls)}
#
#################
### sample a candidate of merge and compute its transition function
#
gera.juncao.QTL<-function(dados,num.QTLs,mat.delinea,vet.coef,pos.qtls,media.alpha,sigma2.alpha,sigma2.vig,media.delta,sigma2.delta,media.mi,sigma2.mi,neta.a,neta.b){
 cramer<-numeric(num.QTLs-1)
 for (i in 1:(num.QTLs-1)) cramer[i]<-abs(cv.test(mat.delinea[,2*i],mat.delinea[,2*(i+1)])) # cramer entre os QTLs vizinhos
 prob_par<-cramer/sum(cramer)
 junta<-rDiscreta(prob_par)
 par_QTL<-c(junta,junta+1)
 efeitos<-c(1/sum(abs(vet.coef[2*junta,1]),abs(vet.coef[2*junta+1,1])),1/sum(abs(vet.coef[2*(junta+1),1]),abs(vet.coef[2*(junta+1)+1,1])))
 prob<-efeitos/sum(efeitos)
 gera<-rDiscreta(prob)
 qtl<-par_QTL[gera]
 #
 pos.qtls.c<-pos.qtls[-qtl]
 mat.delinea.c<-matrix(mat.delinea[,-c((2*qtl),(2*qtl+1))],nrow=nrow(dados))
 vet.coef.c<-matrix(vet.coef[-c((2*qtl),(2*qtl+1))],ncol=1)
 predito<-mat.delinea.c%*%vet.coef.c
 residuos<-dados[,1]-predito
 #
 alpha<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,vet.coef.c[(2*par_QTL[[1]]),1],mat.delinea.c[,(2*par_QTL[[1]])])
 vet.coef.c[(2*par_QTL[[1]]),1]<-alpha[[1]]
 residuos<-dados[,1]-(mat.delinea.c%*%vet.coef.c)
 delta<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,vet.coef.c[(2*par_QTL[[1]])+1,1],mat.delinea.c[,(2*par_QTL[[1]])+1])
 vet.coef.c[(2*par_QTL[[1]])+1,1]<-delta[[1]]
 residuos<-dados[,1]-(mat.delinea.c%*%vet.coef.c)
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef.c[1,1])
 vet.coef.c[1,1]<-mi[[1]]
 predito<-mat.delinea.c%*%vet.coef.c
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 prob_merge<-log(prob[gera])+log(prob_par[junta])
 #
 list(qtl,mat.delinea.c,vet.coef.c,sigma2[[1]],prob_merge,mi[[2]],sigma2[[2]],alpha[[2]],delta[[2]],par_QTL[which(par_QTL!=qtl)],par_QTL,pos.qtls.c)}
#
#################
### compute the acceptance probability of a birth move
#
prob.aceitacao<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphasp,media.alpha,sigma2.alpha,deltasp,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos){
 vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
 priori<-dnorm(misp,media.mi,sqrt(sigma2.mi),log=TRUE)-dnorm(mi,media.mi,sqrt(sigma2.mi),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         dnorm(alphasp,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp,media.delta,sqrt(sigma2.delta),log=TRUE)+
         pri.pos.sp-pri.pos
 trans<-pmerge+plambdamg+postmi+postsigma2-psplit-pmarc-plambdasp-postalphasp-postdeltasp-postmisp-postsigma2sp
 prob.ace<-exp(vero+priori+trans)
 return(prob.ace)}
#
#################
### compute the acceptance probability of a split move
#
prob.aceitacao.split<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphamg,alphasp1,alphasp2,media.alpha,sigma2.alpha,deltamg,deltasp1,deltasp2,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postalpha,postdelta,postmi,postsigma2,postalphasp1,postalphasp2,postdeltasp1,postdeltasp2,postmisp,postsigma2sp,pri.pos.sp,pri.pos){
 vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
 priori<-dnorm(misp,media.mi,sqrt(sigma2.mi),log=TRUE)-dnorm(mi,media.mi,sqrt(sigma2.mi),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         dnorm(alphasp1,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp1,media.delta,sqrt(sigma2.delta),log=TRUE)+
         dnorm(alphasp2,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp2,media.delta,sqrt(sigma2.delta),log=TRUE)-
         (dnorm(alphamg,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltamg,media.delta,sqrt(sigma2.delta),log=TRUE))+
         pri.pos.sp-pri.pos
 trans<-pmerge+plambdamg+postalpha+postdelta+postmi+postsigma2-psplit-pmarc-plambdasp-postalphasp1-postalphasp2-postdeltasp1-postdeltasp2-postmisp-postsigma2sp
 prob.ace<-exp(vero+priori+trans)
 return(prob.ace)}
#
#################
### run the MCMC procedure starting without QTLs - main function
#
DDRJ_mainQTL_nondensemap<-function(dados,loc.marc=seq(from=0,to=(ncol(dados)-2),by=1),directory=getwd(),finsample=1000,burnin=1000,jumps=1, neta.a=0.1,neta.b=0.1,media.mi=0,sigma2.mi=100,sigma2.alpha=100,media.alpha=0,sigma2.delta=100,media.delta=0
){
	# dados = the data matrix with the phenotype variable in the first column and markers' genotype in the following columns, all numeric variables. It should not have missing value and the genotypes must be coded as -1 for recessive homozygotes, 0 for heterozygotes and 1 for dominant homozygotes;
	# loc.marc = a vector with the position (in Morgans) of each marker and the first element should be equal to zero. If it is not available, the algorithm assumes 1 M between the markers which represents almost no association between the markers' genotype. However, we recommend to inform markers' position (in Morgans) for better results and faster convergence;
	# directory = a directory in your computer where the MCMC results will be recorded. If it is not available, the algorithm assumes the default directory in R;
	# finsample = MCMC sample size after burnin and jumps. The default is assumed as 1000. However, we recommend to assume a larger sample for more confident results;
	# burnin = number of iterations that will be discarded in the beginning of the MCMC. The default is assumed as 1000. However, we recommend to assume a larger burnin period (5000 for instance) for better convergence since it is a complex model, specially when the number of available markers is great.
	# jumps = number of iterations that will be discarded between two recorded iterations.The default is assumed as 1 which represents no jump;
	# all remaining parameters of this function are hyperparameters' value of prior Gamma and Normal prior distributions assumed for the precision parameter, general mean, additive and dominance effects. The default option assumes vague prior distribution with large variance for these parameters.
	#
	# initializing the model without QTLs
	#
	num.marc<-ncol(dados)-1
	residuos<-dados[,1]
	sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
	#
	mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,0)[[1]]
	#
	pos.qtls<-NULL # QTLs' positon in the initial model
	num.QTLs<-length(pos.qtls)
	mat.delinea<-matrix(1,nrow(dados),1)
	vet.coef<-matrix(mi,1,1)
	predito<-mat.delinea%*%vet.coef
	residuos<-dados[,1]-predito
	sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
	#
	AmostrasTotal<-burnin+finsample*jumps
	#
	# starting the MCMC iterations
	#
	enableJIT(3)
	for (int in (1:AmostrasTotal)){
		cat('\n', int, 'of', AmostrasTotal, 'iterations. Number of QTLs=', num.QTLs)
		cand.sp.mg<-dec.sp.mg(num.QTLs,num.marc)
		indSpMgtotal<-cand.sp.mg[[1]]
		#	indSpMgtotal= 1 - birth candidate; 2 - death candidate
		#
		##############
		###### birth of a QTL
		##############
		#
		if (indSpMgtotal==1) {candidato<-gera.inclusao.QTL(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,0,0,media.mi,sigma2.mi,media.alpha,sigma2.alpha,media.delta,sigma2.delta,neta.a,neta.b)
			residuossp<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
			sigma2<-sigma2.vig
			sigma2sp<-candidato[[4]]
			mi<-vet.coef[1,1]
			misp<-candidato[[3]][1,1]
			alphasp<-candidato[[3]][(nrow(candidato[[3]])-1),1]
			deltasp<-candidato[[3]][(nrow(candidato[[3]])),1]
			num.QTLssp<-num.QTLs+1
			psplit<-cand.sp.mg[[2]][1]
			pmerge<-dec.sp.mg(num.QTLssp,num.marc)[[2]][2]
			pmarc<-candidato[[5]]
			plambdasp<-candidato[[6]]
			#
			efeito<-NULL
			if (num.QTLs==0) efeito<-0
			if (num.QTLs>0) for (i in 1:num.QTLs) efeito[i]<-1/(abs(vet.coef[2*i,1])+abs(vet.coef[(2*i)+1,1]))
			efeitosp<-1/(abs(alphasp)+abs(deltasp))
			plambdamg<-log(efeitosp/(efeitosp+sum(efeito)))
			#
			res.mi<-residuos+mi
			quo<-(length(residuos)/sigma2sp)+(1/sigma2.mi)
			media<-((sum(res.mi)/sigma2sp)+(media.mi/sigma2.mi))/quo
			variancia<-1/quo
			postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
			#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
			#
			postalphasp<-candidato[[8]]
			postdeltasp<-candidato[[9]]
			postmisp<-candidato[[10]]
			postsigma2sp<-candidato[[11]]
			pos.qtls.sp<-sort(candidato[[12]])
			pri.pos.sp<-dens.priori.loc.qtls(pos.qtls.sp,loc.marc)
			pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
			#
			probace<-prob.aceitacao(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphasp,media.alpha,sigma2.alpha,deltasp,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos)}
		#
		##############
		###### death of a QTL
		##############
		#
		if (indSpMgtotal==2) {candidato<-gera.exclusao.QTL(pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.mi,sigma2.mi,neta.a,neta.b)
			residuosmg<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
			sigma2<-sigma2.vig
			sigma2mg<-candidato[[4]]
			mi<-vet.coef[1,1]
			mimg<-candidato[[3]][1,1]
			num.QTLsmg<-num.QTLs-1
			pmerge<-cand.sp.mg[[2]][2]
			psplit<-dec.sp.mg(num.QTLsmg,num.marc)[[2]][1]
			postmimg<-candidato[[6]]
			postsigma2mg<-candidato[[7]]
			pos.qtls.mg<-sort(candidato[[12]])
			pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
			pri.pos.mg<-dens.priori.loc.qtls(pos.qtls.mg,loc.marc)
			plambdamg<-candidato[[5]]
			alpha<-vet.coef[(2*candidato[[1]]),1]
			delta<-vet.coef[(2*candidato[[1]]+1),1]
			#
			krusk<-prob.selec.marc(dados,residuosmg,pos.qtls.mg,loc.marc)[[3]]
			prob<-krusk/sum(krusk)
			QTLmg<-pos.qtls[candidato[[1]]]
			Marc1<-sum(loc.marc<=QTLmg)
			Marc2<-num.marc-sum(loc.marc>=QTLmg)+1
			loc.Marc1<-dens.pos.qtl(pos.qtls.mg,loc.marc,Marc1,krusk,QTLmg)
			loc.Marc2<-dens.pos.qtl(pos.qtls.mg,loc.marc,Marc2,krusk,QTLmg)
			plambdasp<-log(prob[Marc1]*loc.Marc1+prob[Marc2]*loc.Marc2)
			pmarc<-0
			#
			quo<-(sum(mat.delinea[,2*candidato[[1]]]^2)/sigma2mg)+(1/sigma2.alpha)
			media<-((sum(mat.delinea[,2*candidato[[1]]]*residuosmg)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
			variancia<-1/quo
			postalpha<-dnorm(alpha,media,sqrt(variancia),log=TRUE)
			#
			quo<-(sum(mat.delinea[,2*candidato[[1]]+1]^2)/sigma2mg)+(1/sigma2.delta)
			media<-((sum(mat.delinea[,2*candidato[[1]]+1]*residuosmg)/sigma2mg)+(media.delta/sigma2.delta))/quo
			variancia<-1/quo
			postdelta<-dnorm(delta,media,sqrt(variancia),log=TRUE)
			#
			res.mi<-residuos+vet.coef[1,1]
			quo<-(length(residuos)/sigma2mg)+(1/sigma2.mi)
			media<-((sum(res.mi)/sigma2mg)+(media.mi/sigma2.mi))/quo
			variancia<-1/quo
			postmi<-dnorm(mi,media,sqrt(variancia),log=TRUE)
			#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2<-dgamma((1/sigma2),aa1,bb1,log=TRUE)
			#
			probace<-1/prob.aceitacao(residuos,residuosmg,sigma2mg,sigma2,mi,mimg,media.mi,sigma2.mi,neta.a,neta.b,alpha,media.alpha,sigma2.alpha,delta,media.delta,sigma2.delta,
loc.marc,num.QTLsmg,num.QTLs,psplit,pmerge,pmarc,plambdasp,plambdamg,postmimg,postsigma2mg,postalpha,postdelta,postmi,postsigma2,pri.pos,pri.pos.mg)}
		#
		#####################
		#  update model's parameters
		#####################
		#
		# update the number of QTLs (accepting or rejecting the candidate of birth or death of a QTL)
		#
		probacetotal<-probace
		aux2<-runif(1)
		if (aux2<probace){
			indrejtotal<-0
			pos.qtls<-candidato[[12]]
			num.QTLs<-length(pos.qtls)
			mat.delinea<-candidato[[2]]
			vet.coef<-candidato[[3]]
			sigma2.vig<-candidato[[4]]
			if (num.QTLs>0){
				posicao<-order(pos.qtls)
				pos.qtls<-pos.qtls[posicao]
				posicao2<-1
				for (i in 1:length(posicao)) posicao2<-c(posicao2,posicao[i]*2,posicao[i]*2+1)
				mat.delinea<-mat.delinea[,posicao2]
				vet.coef<-matrix(vet.coef[posicao2,],ncol(mat.delinea),1)}
			residuos<-dados[,1]-(mat.delinea%*%vet.coef)}
		#
		if (aux2>=probace) indrejtotal<-1
		#
		#### evaluate a merge move of two consecutive QTLs
		#
		if (num.QTLs>1){
			candidato<-gera.juncao.QTL(dados,num.QTLs,mat.delinea,vet.coef,pos.qtls,media.alpha,sigma2.alpha,sigma2.vig,media.delta,sigma2.delta,media.mi,sigma2.mi,neta.a,neta.b)
			residuosmg<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
			sigma2<-sigma2.vig
			sigma2mg<-candidato[[4]]
			mi<-vet.coef[1,1]
			mimg<-candidato[[3]][1,1]
			num.QTLsmg<-num.QTLs-1
			pmerge<-0
			psplit<-0
			postmimg<-candidato[[6]]
			postsigma2mg<-candidato[[7]]
			postalphamg<-candidato[[8]]
			postdeltamg<-candidato[[9]]
			pos.qtls.mg<-sort(candidato[[12]])
			pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
			pri.pos.mg<-dens.priori.loc.qtls(pos.qtls.mg,loc.marc)
			plambdamg<-candidato[[5]]
			alphamg<-candidato[[3]][(2*candidato[[11]][1]),1]
			deltamg<-candidato[[3]][(2*candidato[[11]][1]+1),1]
			alphasp1<-vet.coef[(2*candidato[[1]]),1]
			deltasp1<-vet.coef[(2*candidato[[1]]+1),1]
			alphasp2<-vet.coef[(2*candidato[[10]]),1]
			deltasp2<-vet.coef[(2*candidato[[10]]+1),1]
			#
			krusk<-prob.selec.marc(dados,residuosmg,pos.qtls.mg,loc.marc)[[3]]
			prob<-krusk/sum(krusk)
			QTLmg<-pos.qtls[candidato[[1]]]
			Marc1<-sum(loc.marc<=QTLmg)
			Marc2<-num.marc-sum(loc.marc>=QTLmg)+1
			loc.Marc1<-dens.pos.qtl(pos.qtls.mg,loc.marc,Marc1,krusk,QTLmg)
			loc.Marc2<-dens.pos.qtl(pos.qtls.mg,loc.marc,Marc2,krusk,QTLmg)
			plambdasp<-log(prob[Marc1]*loc.Marc1+prob[Marc2]*loc.Marc2)
			pmarc<-0
			#
			quo<-(sum(mat.delinea[,2*candidato[[1]]]^2)/sigma2mg)+(1/sigma2.alpha)
			media<-((sum(mat.delinea[,2*candidato[[1]]]*residuosmg)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
			variancia<-1/quo
			postalpha1<-dnorm(alphasp1,media,sqrt(variancia),log=TRUE)
			#
			quo<-(sum(mat.delinea[,2*candidato[[1]]+1]^2)/sigma2mg)+(1/sigma2.delta)
			media<-((sum(mat.delinea[,2*candidato[[1]]+1]*residuosmg)/sigma2mg)+(media.delta/sigma2.delta))/quo
			variancia<-1/quo
			postdelta1<-dnorm(deltasp1,media,sqrt(variancia),log=TRUE)
			#
			vet_coef_parc<-rbind(candidato[[3]],vet.coef[candidato[[1]]*2,1],vet.coef[candidato[[1]]*2+1,1])
			mat_delinea_parc<-cbind(candidato[[2]],mat.delinea[,candidato[[1]]*2],mat.delinea[,candidato[[1]]*2+1])
			vet_coef_parc<-matrix(vet_coef_parc[-(2*candidato[[11]][1]),],ncol=1)
			mat_delinea_parc<-mat_delinea_parc[,-(2*candidato[[11]][1])]
			residuosparc<-dados[,1]-mat_delinea_parc%*%vet_coef_parc
			quo<-(sum(mat.delinea[,2*candidato[[10]]]^2)/sigma2mg)+(1/sigma2.alpha)
			media<-((sum(mat.delinea[,2*candidato[[10]]]*residuosparc)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
			variancia<-1/quo
			postalpha2<-dnorm(alphasp2,media,sqrt(variancia),log=TRUE)
			#
			vet_coef_parc<-matrix(vet.coef[-(2*candidato[[10]]+1)],ncol=1)
			mat_delinea_parc<-mat.delinea[,-(2*candidato[[10]]+1)]
			residuosparc<-dados[,1]-mat_delinea_parc%*%vet_coef_parc
			quo<-(sum(mat.delinea[,2*candidato[[10]]+1]^2)/sigma2mg)+(1/sigma2.delta)
			media<-((sum(mat.delinea[,2*candidato[[10]]+1]*residuosparc)/sigma2mg)+(media.delta/sigma2.delta))/quo
			variancia<-1/quo
			postdelta2<-dnorm(deltasp2,media,sqrt(variancia),log=TRUE)
			#
			res.mi<-residuos+vet.coef[1,1]
			quo<-(length(residuos)/sigma2mg)+(1/sigma2.mi)
			media<-((sum(res.mi)/sigma2mg)+(media.mi/sigma2.mi))/quo
			variancia<-1/quo
			postmi<-dnorm(mi,media,sqrt(variancia),log=TRUE)
			#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2<-dgamma((1/sigma2),aa1,bb1,log=TRUE)
			#
			probace<-1/prob.aceitacao.split(residuos,residuosmg,sigma2mg,sigma2,mi,mimg,media.mi,sigma2.mi,neta.a,neta.b,alphamg,alphasp1,alphasp2,media.alpha,sigma2.alpha,deltamg,deltasp1,deltasp2,media.delta,
sigma2.delta,loc.marc,num.QTLsmg,num.QTLs,psplit,pmerge,pmarc,plambdasp,plambdamg,postalphamg,postdeltamg,postmimg,postsigma2mg,postalpha1,postalpha2,postdelta1,postdelta2,postmi,postsigma2,pri.pos,pri.pos.mg)
			#
			probacemerge<-probace
			aux2<-runif(1)
			if (aux2<probace){
				pos.qtls<-candidato[[12]]
				num.QTLs<-length(pos.qtls)
				mat.delinea<-candidato[[2]]
				vet.coef<-candidato[[3]]
				sigma2.vig<-candidato[[4]]
				residuos<-dados[,1]-(mat.delinea%*%vet.coef)}}
		#
		#### update QTLs' position using a Metropolis-Hastings algorithm
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
				M.esq<-sum(loc.marc<pos.qtls[i])
				M.dir<-num.marc-(sum(loc.marc>pos.qtls[i]))+1
				loc.cand<-runif(1,min=loc.marc[M.esq],max=loc.marc[M.dir])
				#
				vet.delinea<-mat.delinea
				r1c<-Haldane(abs(loc.cand-loc.marc[M.esq]))
				r2c<-Haldane(abs(loc.cand-loc.marc[M.dir]))
				r12<-Haldane(abs(loc.marc[M.dir]-loc.marc[M.esq]))
				r1<-Haldane(abs(pos.qtls[i]-loc.marc[M.esq]))
				r2<-Haldane(abs(pos.qtls[i]-loc.marc[M.dir]))
				probQTL<-numeric(3)
 				probQTLat<-numeric(3)
				logdens<-numeric(3)
				logdensat<-numeric(3)
				densacumat<-probacumat<-probacum<-densacum<-0
				gen<-c(-1,0,1)
				for (j in 1:nrow(dados)){
					for (l in 1:3){
						probQTL[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1c)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2c))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
						probQTLat[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
						dom<-1-abs(gen[l])
						vet.delinea[j,2*i]<-gen[l]
						vet.delinea[j,(2*i)+1]<-dom
						logdens[l]<-dnorm(dados[j,1],vet.delinea[j,]%*%vet.coef,sqrt(sigma2.vig),log=TRUE)
						logdensat[l]<-dnorm(dados[j,1],vet.delinea[j,]%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
					prob<-exp(probQTL+logdens-max(probQTL+logdens))/sum(exp(probQTL+logdens-max(probQTL+logdens)))
					probat<-exp(probQTLat+logdensat-max(probQTLat+logdensat))/sum(exp(probQTLat+logdensat-max(probQTLat+logdensat)))
					ger.gen<-rDiscreta(prob)
					vet.delinea[j,2*i]<-gen[ger.gen]
					vet.delinea[j,(2*i)+1]<-1-abs(vet.delinea[j,2*i])
					probacum<-probacum+log(prob[ger.gen])
					densacum<-densacum+logdens[ger.gen]
					gen.at<-sum(gen<=mat.delinea[j,2*i])
					probacumat<-probacumat+log(probat[gen.at])
					densacumat<-densacumat+logdensat[gen.at]}
				#
				numer<-denom<-0
				for (j in 1:nrow(dados)){
					numer<-numer+log((calc.prob.gen(dados[j,(M.esq+1)],vet.delinea[j,2*i],r1c)*calc.prob.gen(vet.delinea[j,2*i],dados[j,(M.dir+1)],r2c))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
					denom<-denom+log((calc.prob.gen(dados[j,(M.esq+1)],mat.delinea[j,2*i],r1)*calc.prob.gen(mat.delinea[j,2*i],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))}
				#
				paceit<-exp(densacum+numer+probacumat-densacumat-denom-probacum)
				aux3<-runif(1)
				if (aux3<paceit){
					pos.qtls[i]<-loc.cand
					mat.delinea[,2*i]<-vet.delinea[,2*i]
					mat.delinea[,(2*i)+1]<-vet.delinea[,(2*i)+1]}}}
		#
		#### update QTLs' genotype using a Gibbs sampling step
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
				M.esq<-sum(loc.marc<pos.qtls[i])
				M.dir<-num.marc-(sum(loc.marc>pos.qtls[i]))+1
				r12<-Haldane(abs(loc.marc[M.dir]-loc.marc[M.esq]))
				r1<-Haldane(abs(pos.qtls[i]-loc.marc[M.esq]))
				r2<-Haldane(abs(pos.qtls[i]-loc.marc[M.dir]))
				probQTL<-numeric(3)
				logdens<-numeric(3)
				gen<-c(-1,0,1)
				for (j in 1:nrow(dados)){
					for (l in 1:3){
						probQTL[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
						dom<-1-abs(gen[l])
						vet.delinea<-mat.delinea[j,]
						vet.delinea[2*i]<-gen[l]
						vet.delinea[(2*i)+1]<-dom
						logdens[l]<-dnorm(dados[j,1],vet.delinea%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
					prob<-exp(probQTL+logdens-max(probQTL+logdens))/sum(exp(probQTL+logdens-max(probQTL+logdens)))
					mat.delinea[j,2*i]<-gen[rDiscreta(prob)]
					mat.delinea[j,(2*i)+1]<-1-abs(mat.delinea[j,2*i])}}}
		residuos<-dados[,1]-(mat.delinea%*%vet.coef)
		#
		#### update general mean - Gibbs sampling
		#
		vet.coef[1,1]<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])[[1]]
		residuos<-dados[,1]-(mat.delinea%*%vet.coef)
		#
		#### update the additive and dominance effects - Gibbs sampling
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
				vet.coef[(2*i),1]<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,vet.coef[(2*i),1],mat.delinea[,(2*i)])[[1]]
				residuos<-dados[,1]-(mat.delinea%*%vet.coef)
				vet.coef[(2*i)+1,1]<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,vet.coef[(2*i)+1,1],mat.delinea[,(2*i)+1])[[1]]
				residuos<-dados[,1]-(mat.delinea%*%vet.coef)}}
		#
		#### update the error variance - Gibbs sampling
		#
		sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
		#
		############################## record acceptance's information
		#
		cat('',indSpMgtotal,file=paste(directory,"/ind_birdeath.txt",sep=""),append=T)
		cat('',indrejtotal,file=paste(directory,"/ind_rej_birdeath.txt",sep=""),append=T)
		cat('',round(probacetotal,2),file=paste(directory,"/prob_accep_birdeath.txt",sep="")
,append=T)
		#
		############################## record the sample after burn-in and jumps
		#
		if (int>burnin & int%%jumps==0){
			cat('',num.QTLs,file=paste(directory,"/number_QTLs.txt",sep=""),append=T)
			cat('',pos.qtls,file=paste(directory,"/QTLs_position.txt",sep=""),append=T)
			cat('',vet.coef,file=paste(directory,"/QTLs_effects.txt",sep=""),append=T)
			cat('',sigma2.vig,file=paste(directory,"/error_variance.txt",sep=""),append=T)}}
}
#
###############
###### Running the synthetic example of the original paper with 450 markers
###############
#
library(Rcpp)
library(compiler)
#
directory<-"/Users/Daiane/Downloads" # include a directory in your computer where the data files are recorded
mapa<-scan(paste(directory,"/mapa_locus.txt",sep=""))
gen.marc<-matrix(scan(paste(directory,"/gen_marc.txt",sep="")),nrow=300,byrow=TRUE)
loc.marc<-mapa/100
#
###### QTLs at position 15, 82, 299, 362, 390 ######
#
set.seed(51)
mi.verd<-20
sigma2.verd<-0.5
#
gen.QTL<-gen.marc[,c(15,82,299,362,390)]
res<-rnorm(nrow(gen.marc),0,sqrt(sigma2.verd))
#
ef.adi.verd<-c(-0.60,0.90,0.25,-0.40,0.40)
ef.dom.verd<-c(0.30,0.05,-0.25,0.15,-0.15)
#
fenot<-mi.verd+
ef.adi.verd[1]*gen.QTL[,1]+ef.dom.verd[1]*(1-abs(gen.QTL[,1]))+
ef.adi.verd[2]*gen.QTL[,2]+ef.dom.verd[2]*(1-abs(gen.QTL[,2]))+
ef.adi.verd[3]*gen.QTL[,3]+ef.dom.verd[3]*(1-abs(gen.QTL[,3]))+
ef.adi.verd[4]*gen.QTL[,4]+ef.dom.verd[4]*(1-abs(gen.QTL[,4]))+
ef.adi.verd[5]*gen.QTL[,5]+ef.dom.verd[5]*(1-abs(gen.QTL[,5]))+res
#
dados<-cbind(fenot,gen.marc[,-c(15,82,299,362,390)])
loc.marc<-loc.marc[-c(15,82,299,362,390)]
#
sourceCpp(paste(directory,"/cppKruskall.cpp",sep=""))
#set.seed(100)
DDRJ_mainQTL_nondensemap(dados=dados,loc.marc=loc.marc,directory=directory)
#
#### Analysing the results
#
num.qtls<-scan(file=paste(directory,"/number_QTLs.txt",sep=""))
pos.qtls<-scan(file=paste(directory,"/QTLs_position.txt",sep=""))
coef<-scan(file=paste(directory,"/QTLs_effects.txt",sep=""))
sigma2<-scan(file=paste(directory,"/error_variance.txt",sep=""))
#
finsample<-length(num.qtls)
mat.loc<-matrix(0,nrow=finsample,ncol=max(num.qtls))
mat.coef<-matrix(0,nrow=finsample,ncol=(max(num.qtls)*2+1)) # the first column of this matrix has the MCMC sample for the general mean, the second column has simulated values for the additive effect of the first QTL, the third column has simulated values for the dominance effect of the first QTL, the fourth column has simulated values for the additive effect of the second QTL, the fifth column has simulated values for the dominance effect of the second QTL and so on.
obs<-1
obs2<-1
for (i in 1:finsample){
  if (num.qtls[i]>0){
    mat.loc[i,1:num.qtls[i]]<-pos.qtls[obs:(obs+num.qtls[i]-1)]
    mat.coef[i,1:(num.qtls[i]*2+1)]<-round(coef[obs2:(obs2+(num.qtls[i]*2+1)-1)],3)}
  obs<-obs+num.qtls[i]
  obs2<-obs2+(num.qtls[i]*2+1)}
#
plot(num.qtls)
plot(pos.qtls)
plot(density(pos.qtls))
plot(mat.loc[mat.loc[,1]!=0,1],type='l')
plot(mat.loc[mat.loc[,2]!=0,2],type='l')
plot(mat.loc[mat.loc[,3]!=0,3],type='l')
