###################################
###################################
# This is an R code of the
# DDRJ procedure for identifying main and epistatic QTLs when the available genetic map is not dense enough
# and the method try to identify QTLs and predict their genotype between the available markers.
# It assumes a continuous phenotype as response variable.
# Assumptions: 1 - only one QTL in a marker interval
#              2 - independent individuals from a F2 population
#              3 - additive, dominance and epistasis effect of QTLs influence the phenotype
#              4 - Uniform prior distribution for the number of QTLs, Uniform prior distribution for QTLs'
# location, Normal prior distribution for the general mean and QTLs' effects and Gamma prior distribution for
# precision (1/error variance) parameter, as defined in the main paper published at IEEE/ACM Transactions on Computational Biology and Bioinformatics (2021).
#
###################################
###################################
#
###################################
# Functions
###################################
###################################
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
### compute the probability of selecting an epistasis effect to be added in the model through Kruskal-wallis
#
prob.selec.marc.epist<-function(ep_AA,ep_AD,ep_DA,ep_DD,res1,qtls,loc.marc,krusk_matrix){
  cont<-nrow(krusk_matrix)+1
  Amat<-cbind(res1,ep_AA,ep_AD,ep_DA,ep_DD)
  krusk<-matrix(kruskalcpp(Amat),ncol=4)  
  krusk<-cbind(krusk_matrix,krusk)
  prob<-c(krusk[,3:6])^3/sum(c(krusk[,3:6])^3)
  epis<-rDiscreta(prob)
  tipo<-sum(epis>c(cont-1,(cont-1)*2,(cont-1)*3,(cont-1)*4))
  marc1<-krusk[epis%%(cont-1),1] # pega o resto da divisão
  marc2<-krusk[epis%%(cont-1),2]
  if (epis%%(cont-1)==0){
	marc1<-krusk[(cont-1),1]
	marc2<-krusk[(cont-1),2]}
  list(marc1,marc2,prob[epis],epis,tipo,krusk)}
#
################
### compute the probability of excluding QTLs without an epistasis effect and select a QTL to be dropped
#
prob.excl.marc.epis<-function(num.QTLs,vet.coef,epis.ef,pos.qtls){
	efeitos<-NULL
	for (i in 1:num.QTLs) efeitos[i]<-1/sum(c(abs(vet.coef[(2*i),1]),abs(vet.coef[((2*i)+1),1])))
	QTL_epist<-which(pos.qtls%in%epis.ef[,1:2])
	if (length(QTL_epist)>0) efeitos[QTL_epist]<-0
	proba<-efeitos/sum(efeitos)
	qtl<-rDiscreta(proba)
	list(qtl,log(proba[qtl]))}
#
################
### compute the probability of excluding epistasis effect and select an epistasis to be dropped
#
prob.excl.epist<-function(vet.coef.epist,epis.ef){
	efeito<-1/c(abs(vet.coef.epist))
	proba<-efeito/sum(efeito)
	epist<-rDiscreta(proba)
	list(epist,log(proba[epist]))}
#
################
### calculate the position of the new QTL around the chosen marker
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
 dens<-dbeta(uger,parA,parB,log = TRUE)
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
### compute the position of QTLs involved in the chosen epistasis
#
pos_qtl_epist<-function(qtls,loc.marc,marc1,marc2){
 num.marc<-length(loc.marc)
 num.qtls<-length(qtls)
 dens<-0
 if (marc1==1){
 	if (sum(qtls<loc.marc[2] & qtls>loc.marc[1])==0){
 		lambda_1<-runif(1,min=loc.marc[1],max=loc.marc[2])
 		dens<-dens+dunif(lambda_1,min=loc.marc[1],max=loc.marc[2],log=TRUE)} else lambda_1<-qtls[1]} else{
		if (marc1==num.marc-1){
			if (sum(qtls<loc.marc[num.marc-1] & qtls>loc.marc[num.marc-2])==0){
				lambda_1<-runif(1,min=loc.marc[num.marc-2],max=loc.marc[num.marc-1])
				dens<-dens+dunif(lambda_1,min=loc.marc[num.marc-2],max=loc.marc[num.marc-1],log=TRUE)} else lambda_1<-qtls[which(qtls<loc.marc[num.marc-1] & qtls>loc.marc[num.marc-2])]} else{
 		aux<-rDiscreta(c(0.5,0.5))
 		dens<-dens+log(0.5)
 		if (aux==1){
			if (sum(qtls<loc.marc[marc1] & qtls>loc.marc[marc1-1])==0){
				lambda_1<-runif(1,min=loc.marc[marc1-1],max=loc.marc[marc1])
				dens<-dens+dunif(lambda_1,min=loc.marc[marc1-1],max=loc.marc[marc1],log=TRUE)} else lambda_1<-qtls[which(qtls<loc.marc[marc1] & qtls>loc.marc[marc1-1])]} else{
 			if (sum(qtls<loc.marc[marc1+1] & qtls>loc.marc[marc1])==0){
 				lambda_1<-runif(1,min=loc.marc[marc1],max=loc.marc[marc1+1])
 				dens<-dens+dunif(lambda_1,min=loc.marc[marc1],max=loc.marc[marc1+1],log=TRUE)} else lambda_1<-qtls[which(qtls<loc.marc[marc1+1] & qtls>loc.marc[marc1])]}}}
#
 if (marc2==num.marc){
	if (sum(qtls<loc.marc[num.marc] & qtls>loc.marc[num.marc-1])==0){
		lambda_2<-runif(1,min=loc.marc[num.marc-1],max=loc.marc[num.marc])
		dens<-dens+dunif(lambda_2,min=loc.marc[num.marc-1],max=loc.marc[num.marc],log=TRUE)} else lambda_2<-qtls[num.qtls]} else{
		if (lambda_1<loc.marc[marc2] & lambda_1>loc.marc[marc2-1]){
			if (sum(qtls<loc.marc[marc2+1] & qtls>loc.marc[marc2])==0){
				lambda_2<-runif(1,min=loc.marc[marc2],max=loc.marc[marc2+1])
				dens<-dens+dunif(lambda_2,min=loc.marc[marc2],max=loc.marc[marc2+1],log=TRUE)} else lambda_2<-qtls[which(qtls<loc.marc[marc2+1] & qtls>loc.marc[marc2])]} else{
 			aux<-rDiscreta(c(0.5,0.5))
 			dens<-dens+log(0.5)
 			if (aux==1){
				if (sum(qtls<loc.marc[marc2] & qtls>loc.marc[marc2-1])==0){
					lambda_2<-runif(1,min=loc.marc[marc2-1],max=loc.marc[marc2])
					dens<-dens+dunif(lambda_2,min=loc.marc[marc2-1],max=loc.marc[marc2],log=TRUE)} else lambda_2<-qtls[which(qtls<loc.marc[marc2] & qtls>loc.marc[marc2-1])]} else{
 				if (sum(qtls<loc.marc[marc2+1] & qtls>loc.marc[marc2])==0){
 					lambda_2<-runif(1,min=loc.marc[marc2],max=loc.marc[marc2+1])
 					dens<-dens+dunif(lambda_2,min=loc.marc[marc2],max=loc.marc[marc2+1],log=TRUE)} else lambda_2<-qtls[which(qtls<loc.marc[marc2+1] & qtls>loc.marc[marc2])]}}}
 list(lambda_1,lambda_2,dens)}
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
 dens<-dgamma((1/sigma2),alpha,beta,log = TRUE)
 list(sigma2,dens)}
#
################
### sample from the posterior distribution of effects (mean, additive, dominance and epistasis effects)
#
poster_effect<-function(mean_prior,sigma2_prior,sigma2_cur,residuos,current_effect,covariate){
 res_effect<-residuos+(current_effect*covariate)
 quo<-(sum(covariate^2)/sigma2_cur)+(1/sigma2_prior)
 media<-((sum(covariate*res_effect)/sigma2_cur)+(mean_prior/sigma2_prior))/quo
 variancia<-1/quo
 effect<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(effect,media,sqrt(variancia),log = TRUE)
 list(effect,dens)}
#
################
### choose between an inclusion or exclusion of a QTL or an epistasis
#
dec.inc.exc.ef<-function(num.QTLs,num.marc,epis.ef,pos.qtls){
	dens<-0
	if (num.QTLs==0) {pincl<-1; pexcl<-0} else {if (num.QTLs==(num.marc-1) & (length(epis.ef)/3)==(sum(1:(num.marc-1))*4)) {pincl<-0;pexcl<-1} else {pincl<-pexcl<-1/2}}
	prob<-c(pincl,pexcl)
	ind.inc.exc<-rDiscreta(prob)
	dens<-dens+log(prob[ind.inc.exc])
	if (ind.inc.exc==1){
		if (num.QTLs==(num.marc-1)){ pQTL<-0; pepis<-1} else { pQTL<-1/2; pepis<-1/2}} else {
	if (length(epis.ef)==0) {pQTL<-1;pepis<-0} else {if (sum(pos.qtls%in%epis.ef[,1:2])==num.QTLs) {pQTL<-0;pepis<-1} else {pQTL<-1/2;pepis<-1/2}}} 
	prob<-c(pQTL,pepis) 
	ind.qtl.epis<-rDiscreta(prob) 
	dens<-dens+log(prob[ind.qtl.epis]) 
 list(ind.inc.exc,ind.qtl.epis,dens)}
#
#################
### sample a candidate of birth of a main QTL and compute its transition function
#
gera.inclusao.QTL.2<-function(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef){ 
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
	dados.QTL<-matrix(0,nrow(dados),1)
	probQTL<-numeric(3)
	gen<-c(-1,0,1)
	for (j in 1:nrow(dados)){
		for (i in 1:3) probQTL[i]<-(calc.prob.gen(dados[j,(Marc1+1)],gen[i],r1)*calc.prob.gen(gen[i],dados[j,(Marc2+1)],r2))/calc.prob.gen(dados[j,(Marc1+1)],dados[j,(Marc2+1)],r12)
		dados.QTL[j,1]<-gen[rDiscreta(probQTL)]}
	mat.delinea<-cbind(mat.delinea,dados.QTL)
	#
	alpha<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,1])
	vet.coef<-rbind(vet.coef,alpha[[1]])
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	#
	dados.QTL<-cbind(dados.QTL,(0.5-abs(dados.QTL[,1])))
	delta<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,2])
	mat.delinea<-cbind(mat.delinea,dados.QTL[,2])
	vet.coef<-rbind(vet.coef,delta[[1]])
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	#
	mi<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[1,1],mat.delinea[,1])
	vet.coef[1,1]<-mi[[1]]
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	sigma2<-poster.sigma2(neta.a,neta.b,residuos)
	#
	list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],marcador[[2]],qtl[[2]],0,alpha[[2]],delta[[2]],mi[[2]],sigma2[[2]],c(pos.qtls,qtl[[1]]),mat.delinea.epist,vet.coef.epist,epis.ef)}
#
#################
### draw an epistasis inclusion and calculate its acceptance probability
#
gera.inclusao.epist<-function(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,epis.ef,mat.delinea.epist,vet.coef.epist,ep_AA,ep_AD,ep_DA,ep_DD,krusk_matrix){
	ptransition<-0
	epist<-prob.selec.marc.epist(ep_AA,ep_AD,ep_DA,ep_DD,residuos,pos.qtls,loc.marc,krusk_matrix)
	esc_QTLs_posi<-pos_qtl_epist(pos.qtls,loc.marc,epist[[1]],epist[[2]])
	cand<-c(esc_QTLs_posi[[1]],esc_QTLs_posi[[2]],epist[[5]])
	ptransition<-ptransition+log(epist[[3]])+esc_QTLs_posi[[3]]
	equal<-0
	if (length(epis.ef)>0){
		for (i in 1:nrow(epis.ef)){
			if (sum(epis.ef[i,]==cand)==3) equal<-1}}
	if (equal==1) list(1) else{
		epis.ef<-rbind(epis.ef,cand)
		if (!(esc_QTLs_posi[[1]] %in% pos.qtls)){
			pos.qtls<-c(pos.qtls,esc_QTLs_posi[[1]])
			dQTL<-esc_QTLs_posi[[1]]
			Marc1<-sum(loc.marc<=dQTL)
			Marc2<-length(loc.marc)-(sum(loc.marc>=dQTL)-1)
			dM1<-loc.marc[Marc1]
			dM2<-loc.marc[Marc2]
			r12<-Haldane(abs(dM2-dM1))
			r1<-Haldane(abs(dQTL-dM1))
			r2<-Haldane(abs(dQTL-dM2))
			#
			dados.QTL<-matrix(0,nrow(dados),1)
			probQTL<-numeric(3)
			gen<-c(-1,0,1)
			for (j in 1:nrow(dados)){
				for (i in 1:3) probQTL[i]<-(calc.prob.gen(dados[j,(Marc1+1)],gen[i],r1)*calc.prob.gen(gen[i],dados[j,(Marc2+1)],r2))/calc.prob.gen(dados[j,(Marc1+1)],dados[j,(Marc2+1)],r12)
				dados.QTL[j,1]<-gen[rDiscreta(probQTL)]}
			mat.delinea<-cbind(mat.delinea,dados.QTL)
			#
			alpha<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,1])
			ptransition<-ptransition+alpha[[2]]
			vet.coef<-rbind(vet.coef,alpha[[1]])
			predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
			residuos<-dados[,1]-predito
			#
			dados.QTL<-cbind(dados.QTL,(0.5-abs(dados.QTL[,1])))
			delta<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,2])
			ptransition<-ptransition+delta[[2]]
			mat.delinea<-cbind(mat.delinea,dados.QTL[,2])
			vet.coef<-rbind(vet.coef,delta[[1]])
			predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
			residuos<-dados[,1]-predito}
		#
		if (!(esc_QTLs_posi[[2]] %in% pos.qtls)){
			pos.qtls<-c(pos.qtls,esc_QTLs_posi[[2]])
			dQTL<-esc_QTLs_posi[[2]]
			Marc1<-sum(loc.marc<=dQTL)
			Marc2<-length(loc.marc)-(sum(loc.marc>=dQTL)-1)
			dM1<-loc.marc[Marc1]
			dM2<-loc.marc[Marc2]
			r12<-Haldane(abs(dM2-dM1))
			r1<-Haldane(abs(dQTL-dM1))
			r2<-Haldane(abs(dQTL-dM2))
			#
			dados.QTL<-matrix(0,nrow(dados),1)
			probQTL<-numeric(3)
			gen<-c(-1,0,1)
			for (j in 1:nrow(dados)){
				for (i in 1:3) probQTL[i]<-(calc.prob.gen(dados[j,(Marc1+1)],gen[i],r1)*calc.prob.gen(gen[i],dados[j,(Marc2+1)],r2))/calc.prob.gen(dados[j,(Marc1+1)],dados[j,(Marc2+1)],r12)
				dados.QTL[j,1]<-gen[rDiscreta(probQTL)]}
			mat.delinea<-cbind(mat.delinea,dados.QTL)
			#
			alpha<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,1])
			ptransition<-ptransition+alpha[[2]]
			vet.coef<-rbind(vet.coef,alpha[[1]])
			predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
			residuos<-dados[,1]-predito
			#
			dados.QTL<-cbind(dados.QTL,(0.5-abs(dados.QTL[,1])))
			delta<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,dados.QTL[,2])
			ptransition<-ptransition+delta[[2]]
			mat.delinea<-cbind(mat.delinea,dados.QTL[,2])
			vet.coef<-rbind(vet.coef,delta[[1]])
			predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
			residuos<-dados[,1]-predito}
		#
		if (epist[[5]]==0) covar_epist<-mat.delinea[,2*which(pos.qtls==cand[1])]*mat.delinea[,2*which(pos.qtls==cand[2])]
		if (epist[[5]]==1) covar_epist<-mat.delinea[,2*which(pos.qtls==cand[1])]*mat.delinea[,(2*which(pos.qtls==cand[2]))+1]
		if (epist[[5]]==2) covar_epist<-mat.delinea[,(2*which(pos.qtls==cand[1]))+1]*mat.delinea[,2*which(pos.qtls==cand[2])]
		if (epist[[5]]==3) covar_epist<-mat.delinea[,(2*which(pos.qtls==cand[1]))+1]*mat.delinea[,(2*which(pos.qtls==cand[2]))+1]
		mat.delinea.epist<-cbind(mat.delinea.epist,covar_epist)
		beta_epist<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,covar_epist)
		ptransition<-ptransition+beta_epist[[2]]
		vet.coef.epist<-rbind(vet.coef.epist,beta_epist[[1]])
		predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
		residuos<-dados[,1]-predito
		#
		mi<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[1,1],mat.delinea[,1])
		vet.coef[1,1]<-mi[[1]]
		ptransition<-ptransition+mi[[2]]
		predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
		residuos<-dados[,1]-predito
		sigma2<-poster.sigma2(neta.a,neta.b,residuos)
		ptransition<-ptransition+sigma2[[2]]
		list(0,mat.delinea,vet.coef,sigma2[[1]],ptransition,0,0,0,0,0,0,pos.qtls,mat.delinea.epist,vet.coef.epist,epis.ef)
		}}
#
#################
### drop a QTL without an epistasis from the model
#
gera.exclusao.QTL.2<-function(dados,pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef){
	num.QTLs<-length(pos.qtls)
	qtl<-prob.excl.marc.epis(num.QTLs,vet.coef,epis.ef,pos.qtls)
 #
	pos.qtls<-pos.qtls[-qtl[[1]]]
	mat.delinea<-matrix(mat.delinea[,-c((2*qtl[[1]]),(2*qtl[[1]]+1))],nrow=nrow(mat.delinea))
	vet.coef<-matrix(vet.coef[-c((2*qtl[[1]]),(2*qtl[[1]]+1))],ncol=1)
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
 #
 	mi<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[1,1],mat.delinea[,1])
	vet.coef[1,1]<-mi[[1]]
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
	list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],qtl[[2]],mi[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],pos.qtls,mat.delinea.epist,vet.coef.epist,epis.ef)}
#
#################
### drop an epistasis effect from the model
#
gera.exclusao.epist<-function(dados,pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef){
	epist<-prob.excl.epist(vet.coef.epist,epis.ef)
 #
	epis.ef<-matrix(epis.ef[-epist[[1]],],ncol=3)
	mat.delinea.epist<-matrix(mat.delinea.epist[,-epist[[1]]],nrow=nrow(mat.delinea.epist))
	vet.coef.epist<-matrix(vet.coef.epist[-epist[[1]]],ncol=1)
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
 #
 	mi<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[1,1],mat.delinea[,1])
	vet.coef[1,1]<-mi[[1]]
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
	list(epist[[1]],mat.delinea,vet.coef,sigma2[[1]],epist[[2]],mi[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],pos.qtls,mat.delinea.epist,vet.coef.epist,epis.ef)}
#
#################
### acceptance probability for a QTL birth
#
prob.aceitacao.2<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.prior,sigma2.prior,neta.a,neta.b,alphasp,deltasp,loc.marc,num.QTLs,num.QTLssp,pbirthQTL,pdeathQTL,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos){
 vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
 priori<-dnorm(misp,media.prior,sqrt(sigma2.prior),log=TRUE)-dnorm(mi,media.prior,sqrt(sigma2.prior),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         dnorm(alphasp,media.prior,sqrt(sigma2.prior),log=TRUE)+dnorm(deltasp,media.prior,sqrt(sigma2.prior),log=TRUE)+
         pri.pos.sp-pri.pos
 trans<-pdeathQTL+plambdamg+postmi+postsigma2-pbirthQTL-pmarc-plambdasp-postalphasp-postdeltasp-postmisp-postsigma2sp
 prob.ace<-exp(vero+priori+trans)
 return(prob.ace)}
#
#################
### acceptance probability for an epistasis birth
#
prob.aceitacao.epis<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.prior,sigma2.prior,neta.a,neta.b,pbirthepis,pdeathepis,plambdamg,postmi,postsigma2,ptransi,pri.pos.sp,pri.pos,vet.coef.sp,vet.coef.epist.sp,pos.qtls,pos.qtls.sp){	
	vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
	trans<-pdeathepis+plambdamg+postmi+postsigma2-pbirthepis-ptransi
	QTLsamais<-length(pos.qtls.sp)-length(pos.qtls)
	priori<-dnorm(misp,media.prior,sqrt(sigma2.prior),log=TRUE)-dnorm(mi,media.prior,sqrt(sigma2.prior),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         pri.pos.sp-pri.pos+
         dnorm(vet.coef.epist.sp[nrow(vet.coef.epist.sp),1],media.prior,sqrt(sigma2.prior),log=TRUE)
	if (QTLsamais>0) priori<-priori+sum(dnorm(vet.coef.sp[(nrow(vet.coef.sp)-(QTLsamais*2)+1):nrow(vet.coef.sp),1],media.prior,sqrt(sigma2.prior),log=TRUE))
	prob.ace<-exp(vero+priori+trans)
	return(prob.ace)}
#
#################
### acceptance probability for an epistasis exclusion
#
prob.exclu.epis<-function(residuosmg,residuos,sigma2,sigma2mg,mimg,mi,media.prior,sigma2.prior,neta.a,neta.b,pbirthepis,pdeathepis,plambdasp,postmi,postmimg,postbetamg,postsigma2,postsigma2mg,betaepist,pepist){
	vero<-sum(dnorm(residuosmg,0,sqrt(sigma2mg),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
	trans<-pbirthepis+plambdasp+postbetamg+postmimg+postsigma2mg-pdeathepis-pepist-postmi-postsigma2
	priori<-dnorm(mimg,media.prior,sqrt(sigma2.prior),log=TRUE)-dnorm(mi,media.prior,sqrt(sigma2.prior),log=TRUE)+
         dgamma((1/sigma2mg),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)-
         dnorm(betaepist,media.prior,sqrt(sigma2.prior),log=TRUE)
	prob.ace<-exp(vero+priori+trans)
	return(prob.ace)}
#
#################
### draw the candidate model and calculate its acceptance probability
#
gera.candidato<-function(chos.move,dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef,ep_AA,ep_AD,ep_DA,ep_DD,krusk_matrix){
	num.QTLs<-length(pos.qtls)
	num.marc<-ncol(dados)-1
	if (chos.move[[1]]==1){
		if (chos.move[[2]]==1){
			candidato<-gera.inclusao.QTL.2(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef)
			residuossp<-dados[,1]-(cbind(candidato[[2]],candidato[[13]])%*%rbind(candidato[[3]],candidato[[14]]))
			sigma2<-sigma2.vig
			sigma2sp<-candidato[[4]]
			mi<-vet.coef[1,1]
			misp<-candidato[[3]][1,1]
			alphasp<-candidato[[3]][(nrow(candidato[[3]])-1),1]
			deltasp<-candidato[[3]][(nrow(candidato[[3]])),1]
			num.QTLssp<-num.QTLs+1
			pos.qtls.sp<-sort(candidato[[12]])
			pbirthQTL<-chos.move[[3]]
			if (num.QTLssp==(num.marc-1) & (length(epis.ef)/3)==(sum(1:(num.marc-1))*4)) pexcl<-1 else pexcl<-1/2   
			if (length(epis.ef)==0) pQTL<-1 else { if (sum(pos.qtls.sp%in%epis.ef[,1:2])==num.QTLssp) pQTL<-0 else pQTL<-1/2}
			pdeathQTL<-log(pexcl)+log(pQTL)
			pmarc<-candidato[[5]]
			plambdasp<-candidato[[6]]
		#
			efeito<-NULL
			if (num.QTLs==0) efeito<-0
			if (num.QTLs>0) for (i in 1:num.QTLs) efeito[i]<-1/(abs(vet.coef[2*i,1])+abs(vet.coef[(2*i)+1,1]))
			efeitosp<-1/(abs(alphasp)+abs(deltasp))
			plambdamg<-log(efeitosp/(efeitosp+sum(efeito[-which(pos.qtls%in%epis.ef[,1:2])])))
		#
			res.mi<-residuos+mi
			quo<-(length(residuos)/sigma2sp)+(1/sigma2.prior)
			media<-((sum(res.mi)/sigma2sp)+(media.prior/sigma2.prior))/quo
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
			pri.pos.sp<-dens.priori.loc.qtls(pos.qtls.sp,loc.marc)
			pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
		#
			probace<-prob.aceitacao.2(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.prior,sigma2.prior,neta.a,neta.b,alphasp,deltasp,loc.marc,num.QTLs,num.QTLssp,pbirthQTL,pdeathQTL,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos)} else {
		#
			candidato<-gera.inclusao.epist(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,epis.ef,mat.delinea.epist,vet.coef.epist,ep_AA,ep_AD,ep_DA,ep_DD,krusk_matrix)
			if (length(candidato)>1){
			#
			residuossp<-dados[,1]-(cbind(candidato[[2]],candidato[[13]])%*%rbind(candidato[[3]],candidato[[14]]))
			sigma2<-sigma2.vig
			sigma2sp<-candidato[[4]]
			mi<-vet.coef[1,1]
			misp<-candidato[[3]][1,1]
			pos.qtls.sp<-sort(candidato[[12]])   
			num.QTLssp<-length(pos.qtls.sp)
			pbirthepis<-chos.move[[3]]
			epis.efsp<-candidato[[15]]
			if (num.QTLssp==(num.marc-1) & (length(epis.efsp)/3)==(sum(1:(num.marc-1))*4)) pexcl<-1 else pexcl<-1/2
			if (length(epis.efsp)==0) pepis<-0 else {if (sum(pos.qtls.sp%in%epis.efsp[,1:2])==num.QTLssp) pepis<-1 else pepis<-1/2}
			pdeathepis<-log(pexcl)+log(pepis)
		#
			if (length(epis.efsp)==0) efeito<-0
			if (length(epis.efsp)>0) efeito<-1/c(abs(candidato[[14]]))
			plambdamg<-log(efeito/sum(efeito))[length(efeito)]
		#
			res.mi<-residuos+mi
			quo<-(length(residuos)/sigma2sp)+(1/sigma2.prior)
			media<-((sum(res.mi)/sigma2sp)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
		#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
		#
			pri.pos.sp<-dens.priori.loc.qtls(pos.qtls.sp,loc.marc)
			pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
		#
			probace<-prob.aceitacao.epis(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.prior,sigma2.prior,neta.a,neta.b,pbirthepis,pdeathepis,plambdamg,postmi,postsigma2,candidato[[5]],pri.pos.sp,pri.pos,candidato[[3]],candidato[[14]],pos.qtls,pos.qtls.sp)} else probace<-0}} else {
		#
		if (chos.move[[2]]==1){
			candidato<-gera.exclusao.QTL.2(dados,pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef)
			residuosmg<-dados[,1]-(cbind(candidato[[2]],candidato[[13]])%*%rbind(candidato[[3]],candidato[[14]]))		
			sigma2<-sigma2.vig
			sigma2mg<-candidato[[4]]
			mi<-vet.coef[1,1]
			mimg<-candidato[[3]][1,1]
			pos.qtls.mg<-sort(candidato[[12]])
			num.QTLsmg<-length(pos.qtls.mg)
			pdeathQTL<-chos.move[[3]]
			if (num.QTLsmg==0) pincl<-1 else pincl<-1/2   
			pbirthQTL<-log(pincl)+log(0.5)
			postmimg<-candidato[[6]]
			postsigma2mg<-candidato[[7]]
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
			quo<-(sum(mat.delinea[,2*candidato[[1]]]^2)/sigma2mg)+(1/sigma2.prior)
			media<-((sum(mat.delinea[,2*candidato[[1]]]*residuosmg)/sigma2mg)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postalpha<-dnorm(alpha,media,sqrt(variancia),log = TRUE)
		#
			quo<-(sum(mat.delinea[,2*candidato[[1]]+1]^2)/sigma2mg)+(1/sigma2.prior)
			media<-((sum(mat.delinea[,2*candidato[[1]]+1]*residuosmg)/sigma2mg)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postdelta<-dnorm(delta,media,sqrt(variancia),log = TRUE)
		#
			res.mi<-residuos+vet.coef[1,1]
			quo<-(length(residuos)/sigma2mg)+(1/sigma2.prior)
			media<-((sum(res.mi)/sigma2mg)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
		#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
		#
			probace<-1/prob.aceitacao.2(residuos,residuosmg,sigma2mg,sigma2,mi,mimg,media.prior,sigma2.prior,neta.a,neta.b,alpha,delta,
loc.marc,num.QTLsmg,num.QTLs,pbirthQTL,pdeathQTL,pmarc,plambdasp,plambdamg,postmimg,postsigma2mg,postalpha,postdelta,postmi,postsigma2,pri.pos,pri.pos.mg)} else{
		#
			candidato<-gera.exclusao.epist(dados,pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef)
			residuosmg<-dados[,1]-(cbind(candidato[[2]],candidato[[13]])%*%rbind(candidato[[3]],candidato[[14]]))
			sigma2<-sigma2.vig
			sigma2mg<-candidato[[4]]
			mi<-vet.coef[1,1]
			mimg<-candidato[[3]][1,1]
			pos.qtls.mg<-sort(candidato[[12]])
			num.QTLsmg<-length(pos.qtls.mg)
			pdeathepis<-chos.move[[3]]
			epis.efmg<-candidato[[15]]		
			if (num.QTLsmg==0) pincl<-1 else pincl<-1/2
			if (num.QTLsmg==(num.marc-1)) pepis<-1 else pepis<-1/2
			pbirthepis<-log(pincl)+log(pepis)
		#
			krusk<-matrix(kruskalcpp(cbind(residuosmg,ep_AA,ep_AD,ep_DA,ep_DD)),ncol=4)  
			krusk<-cbind(krusk_matrix,krusk)
			cons_norma<-sum(c(krusk[,3:6])^3)
			tipo<-epis.ef[candidato[[1]],3]
			qtl1<-epis.ef[candidato[[1]],1]
			qtl2<-epis.ef[candidato[[1]],2]
			marcadores<-c(sum(loc.marc<=qtl1),num.marc-(sum(loc.marc>=qtl1))+1,sum(loc.marc<=qtl2),num.marc-(sum(loc.marc>=qtl2))+1)
			prob_soma<-0
			for (i in 1:2){
				for (j in 3:4){
					p1<-0
					p2<-0
					if (marcadores[i]<marcadores[j]){
						if (marcadores[i]==1 |marcadores[i]==(num.marc-1)) p1<-1 else p1<-0.5
						if (marcadores[j]==num.marc | (marcadores[j]==(marcadores[i]+1) & marcadores[j]==marcadores[j-1])) p2<-1 else p2<-0.5
						prob_soma<-prob_soma+(p1*p2*(krusk[which(marcadores[i]==krusk[,1] & marcadores[j]==krusk[,2]),tipo+3]^3/cons_norma))}}}
			plambdasp<-log(prob_soma)
		#
			quo<-(sum(mat.delinea.epist[,candidato[[1]]]^2)/sigma2mg)+(1/sigma2.prior)
			media<-((sum(mat.delinea.epist[,candidato[[1]]]*residuosmg)/sigma2mg)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postbetamg<-dnorm(vet.coef.epist[candidato[[1]],1],media,sqrt(variancia),log = TRUE)
		#
			res.mi<-residuos+vet.coef[1,1]
			quo<-(length(residuos)/sigma2mg)+(1/sigma2.prior)
			media<-((sum(res.mi)/sigma2mg)+(media.prior/sigma2.prior))/quo
			variancia<-1/quo
			postmimg<-dnorm(mi,media,sqrt(variancia),log = TRUE)
		#
			aa1<-(length(residuos)/2)+neta.a
			bb1<-(sum(residuos^2)/2)+neta.b
			postsigma2mg<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
		#
			probace<-prob.exclu.epis(residuosmg,residuos,sigma2,sigma2mg,mimg,mi,media.prior,sigma2.prior,neta.a,neta.b,pbirthepis,pdeathepis,plambdasp,candidato[[6]],postmimg,postbetamg,candidato[[7]],postsigma2mg,vet.coef.epist[candidato[[1]],1],candidato[[5]])}}
		#
	list(candidato,probace)}
#
#################
### run the MCMC procedure starting without QTLs and epistasis - main function
#
DDRJ_mainepistQTL_nondensemap<-function(dados,loc.marc=seq(from=0,to=(ncol(dados)-2),by=1),directory=getwd(),finsample=1000,burnin=1000,jumps=1,neta.a=0.1,neta.b=0.1,sigma2.prior=100,media.prior=0){
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
	enableJIT(3)
	krusk_matrix<-matrix(0,ncol=2,nrow=sum(1:(num.marc-1)))
	ep_AA<-NULL
	ep_AD<-NULL
	ep_DA<-NULL
	ep_DD<-NULL
	cont<-1
	for (i in 2:(ncol(dados)-1)){
		for (j in (i+1):ncol(dados)){
			krusk_matrix[cont,1]<-(i-1)
			krusk_matrix[cont,2]<-(j-1)
			ep_AA<-cbind(ep_AA,dados[,i]*dados[,j])
  			ep_AD<-cbind(ep_AD,dados[,i]*(0.5-abs(dados[,j])))
			ep_DA<-cbind(ep_DA,(0.5-abs(dados[,i]))*dados[,j])
			ep_DD<-cbind(ep_DD,(0.5-abs(dados[,i]))*(0.5-abs(dados[,j])))
			cont<-cont+1}}
	ep_AD<-ep_AD/0.5
	ep_DA<-ep_DA/0.5
	ep_DD<-ep_DD/0.25
	#
	residuos<-dados[,1]
	sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
	pos.qtls<-NULL
	num.QTLs<-length(pos.qtls)
	mat.delinea<-matrix(1,nrow(dados),1)
	mat.delinea.epist<-NULL
	vet.coef.epist<-NULL
	epis.ef<-NULL
	mi<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,0,mat.delinea[,1])[[1]]
	vet.coef<-matrix(mi,1,1)
	predito<-cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist)
	residuos<-dados[,1]-predito
	sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
	#
	AmostrasTotal<-burnin+finsample*jumps
	#
	enableJIT(3)
	for (int in (1:AmostrasTotal)){
		chos.move<-dec.inc.exc.ef(num.QTLs,num.marc,epis.ef,pos.qtls)
		cat('',chos.move[[1]],file=paste(directory,"/inclu_exclu.txt",sep=""),append=T)
		cat('',chos.move[[2]],file=paste(directory,"/qtl_epist.txt",sep=""),append=T)
		candidatof<-gera.candidato(chos.move,dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,media.prior,sigma2.prior,neta.a,neta.b,mat.delinea.epist,vet.coef.epist,epis.ef,ep_AA,ep_AD,ep_DA,ep_DD,krusk_matrix)
		if (length(candidatof[[1]])==1){
			cat('',99,file=paste(directory,"/accep_prob.txt",sep=""),append=T)} else {
				cat('',round(candidatof[[2]],3),file=paste(directory,"/accep_prob.txt",sep=""),append=T)
				#
				aux2<-runif(1)
				if (aux2<candidatof[[2]]){
					cat('',1,file=paste(directory,"/ind_accep.txt",sep=""),append=T)
					pos.qtls<-candidatof[[1]][[12]]
					num.QTLs<-length(pos.qtls)
					mat.delinea<-candidatof[[1]][[2]]
					vet.coef<-candidatof[[1]][[3]]
					vet.coef.epist<-candidatof[[1]][[14]]
					sigma2.vig<-candidatof[[1]][[4]]
					epis.ef<-candidatof[[1]][[15]]
					if (num.QTLs>0){
						posicao<-order(pos.qtls)
						pos.qtls<-pos.qtls[posicao]
						posicao2<-1
						for (i in 1:length(posicao)) posicao2<-c(posicao2,posicao[i]*2,posicao[i]*2+1)
    					mat.delinea<-mat.delinea[,posicao2]
    					vet.coef<-matrix(vet.coef[posicao2,],ncol(mat.delinea),1)
    					if (length(epis.ef)>0){
							posicao<-order(epis.ef[,1],epis.ef[,2],epis.ef[,3])
							epis.ef<-matrix(epis.ef[posicao,],ncol=3)
							vet.coef.epist<-matrix(vet.coef.epist[posicao,],ncol=1)
							mat.delinea.epist<-NULL
							for (i in 1:nrow(epis.ef)){
								if (epis.ef[i,3]==0) mat.delinea.epist<-cbind(mat.delinea.epist,(mat.delinea[,2*which(pos.qtls==epis.ef[i,1])]*mat.delinea[,2*which(pos.qtls==epis.ef[i,2])]))
								if (epis.ef[i,3]==1) mat.delinea.epist<-cbind(mat.delinea.epist,(mat.delinea[,2*which(pos.qtls==epis.ef[i,1])]*mat.delinea[,(2*which(pos.qtls==epis.ef[i,2]))+1]))
								if (epis.ef[i,3]==2) mat.delinea.epist<-cbind(mat.delinea.epist,(mat.delinea[,(2*which(pos.qtls==epis.ef[i,1]))+1]*mat.delinea[,2*which(pos.qtls==epis.ef[i,2])]))
								if (epis.ef[i,3]==3) mat.delinea.epist<-cbind(mat.delinea.epist,(mat.delinea[,(2*which(pos.qtls==epis.ef[i,1]))+1]*mat.delinea[,(2*which(pos.qtls==epis.ef[i,2]))+1]))}} else mat.delinea.epist<-NULL}
					residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))} else cat('',0,file=paste(directory,"/ind_accep.txt",sep=""),append=T)}
		#
		cat('\n', int, 'of', AmostrasTotal, 'iterations. Number of QTLs=', num.QTLs, '. Number of epistasis=', nrow(epis.ef))
		#
		#### uptade QTLs positions - Metropolis Hasting
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
				esta_epist<-sort(c(which(pos.qtls[i]%in%epis.ef[,1]),which(pos.qtls[i]%in%epis.ef[,2])))
				M.esq<-sum(loc.marc<pos.qtls[i])
				M.dir<-num.marc-(sum(loc.marc>pos.qtls[i]))+1
				if (M.dir<num.marc & sum(pos.qtls >loc.marc[M.dir] & pos.qtls <loc.marc[M.dir+1])==0) M.dir<-M.dir+1
				if (M.esq>1 & sum(pos.qtls<loc.marc[M.esq] & pos.qtls>loc.marc[M.esq-1])==0) M.esq<-M.esq-1
				loc.cand<-runif(1,min=loc.marc[M.esq],max=loc.marc[M.dir])
				ptrans<-dunif(loc.cand,min=loc.marc[M.esq],max=loc.marc[M.dir],log=TRUE)
				#
				vet.delinea<-mat.delinea
				vet.delinea.epist<-mat.delinea.epist
				M.esq2<-sum(loc.marc<loc.cand)
				M.dir2<-num.marc-(sum(loc.marc>loc.cand))+1
				M.esq<-sum(loc.marc<pos.qtls[i])
				M.dir<-num.marc-(sum(loc.marc>pos.qtls[i]))+1
				r1c<-Haldane(abs(loc.cand-loc.marc[M.esq2]))
				r2c<-Haldane(abs(loc.cand-loc.marc[M.dir2]))
				r12c<-Haldane(abs(loc.marc[M.dir2]-loc.marc[M.esq2]))
				r1<-Haldane(abs(pos.qtls[i]-loc.marc[M.esq]))
				r2<-Haldane(abs(pos.qtls[i]-loc.marc[M.dir]))
				r12<-Haldane(abs(loc.marc[M.dir]-loc.marc[M.esq]))
				probQTL<-numeric(3)
				probQTLat<-numeric(3)
				logdens<-numeric(3)
				logdensat<-numeric(3)
				densacumat<-0
				probacumat<-0
				probacum<-0
				densacum<-0
				gen<-c(-1,0,1)
				for (j in 1:nrow(dados)){
					for (l in 1:3){
						probQTL[l]<-log((calc.prob.gen(dados[j,(M.esq2+1)],gen[l],r1c)*calc.prob.gen(gen[l],dados[j,(M.dir2+1)],r2c))/calc.prob.gen(dados[j,(M.esq2+1)],dados[j,(M.dir2+1)],r12c))
						probQTLat[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
						dom<-0.5-abs(gen[l])
						vet.delinea[j,2*i]<-gen[l]
						vet.delinea[j,(2*i)+1]<-dom
						if (length(esta_epist)>0){
							for(ep in esta_epist){
								if (epis.ef[ep,3]==0) vet.delinea.epist[j,ep]<-vet.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]
								if (epis.ef[ep,3]==3) vet.delinea.epist[j,ep]<-vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
								if (epis.ef[ep,3]==1) vet.delinea.epist[j,ep]<-vet.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
								if (epis.ef[ep,3]==2) vet.delinea.epist[j,ep]<-vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]}}
						logdens[l]<-logdensat[l]<-dnorm(dados[j,1],c(vet.delinea[j,],vet.delinea.epist[j,])%*%rbind(vet.coef,vet.coef.epist),sqrt(sigma2.vig),log=TRUE)}
					aux4<-exp(probQTL+logdens-max(probQTL+logdens))
					prob<-aux4/sum(aux4)
					aux5<-exp(probQTLat+logdensat-max(probQTLat+logdensat))
					probat<-aux5/sum(aux5)
					ger.gen<-rDiscreta(prob)
					vet.delinea[j,2*i]<-gen[ger.gen]
					vet.delinea[j,(2*i)+1]<-0.5-abs(vet.delinea[j,2*i])
					if (length(esta_epist)>0){
						for(ep in esta_epist){
							if (epis.ef[ep,3]==0) vet.delinea.epist[j,ep]<-vet.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]
							if (epis.ef[ep,3]==3) vet.delinea.epist[j,ep]<-vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
							if (epis.ef[ep,3]==1) vet.delinea.epist[j,ep]<-vet.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
							if (epis.ef[ep,3]==2) vet.delinea.epist[j,ep]<-vet.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]}}
					probacum<-probacum+log(prob[ger.gen])
					densacum<-densacum+logdens[ger.gen]
					gen.at<-sum(gen<=mat.delinea[j,2*i])
					probacumat<-probacumat+log(probat[gen.at])
					densacumat<-densacumat+logdensat[gen.at]}
					#
				numer<-denom<-0
				for (j in 1:nrow(dados)){
					numer<-numer+log((calc.prob.gen(dados[j,(M.esq2+1)],vet.delinea[j,2*i],r1c)*calc.prob.gen(vet.delinea[j,2*i],dados[j,(M.dir2+1)],r2c))/calc.prob.gen(dados[j,(M.esq2+1)],dados[j,(M.dir2+1)],r12c))
					denom<-denom+log((calc.prob.gen(dados[j,(M.esq+1)],mat.delinea[j,2*i],r1)*calc.prob.gen(mat.delinea[j,2*i],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))}
				#
				pos.qtls.aux<-pos.qtls
				pos.qtls.aux[i]<-loc.cand
				if (M.dir2<num.marc & sum(pos.qtls.aux >loc.marc[M.dir2] & pos.qtls.aux <loc.marc[M.dir2+1])==0) M.dir2<-M.dir2+1
				if (M.esq2>1 & sum(pos.qtls.aux<loc.marc[M.esq2] & pos.qtls.aux>loc.marc[M.esq2-1])==0) M.esq2<-M.esq2-1
				ptransat<-dunif(pos.qtls[i],min=loc.marc[M.esq2],max=loc.marc[M.dir2],log=TRUE)  
				#
				paceit<-exp(densacum+numer+probacumat+ptransat-densacumat-denom-probacum-ptrans)
				aux3<-runif(1)
				if (aux3<paceit){
					lixo<-pos.qtls[i]
					pos.qtls[i]<-loc.cand
					mat.delinea[,2*i]<-vet.delinea[,2*i]
					mat.delinea[,(2*i)+1]<-vet.delinea[,(2*i)+1]
					if (length(esta_epist)>0){			
						mat.delinea.epist[,esta_epist]<-vet.delinea.epist[,esta_epist]
						epis.ef[epis.ef[,1]==lixo,1]<-loc.cand
						epis.ef[epis.ef[,2]==lixo,2]<-loc.cand}}}}
		#
		#### update QTLs' genotype
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
   				esta_epist<-sort(c(which(pos.qtls[i]%in%epis.ef[,1]),which(pos.qtls[i]%in%epis.ef[,2])))
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
						dom<-0.5-abs(gen[l])
						vet.delinea<-mat.delinea[j,]
						vet.delinea.epist<-mat.delinea.epist[j,]
						vet.delinea[2*i]<-gen[l]
						vet.delinea[(2*i)+1]<-dom
						if (length(esta_epist)>0){
							for(ep in esta_epist){
								if (epis.ef[ep,3]==0) vet.delinea.epist[ep]<-vet.delinea[2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[2*which(pos.qtls==epis.ef[ep,2])]
								if (epis.ef[ep,3]==3) vet.delinea.epist[ep]<-vet.delinea[(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[(2*which(pos.qtls==epis.ef[ep,2]))+1]
								if (epis.ef[ep,3]==1) vet.delinea.epist[ep]<-vet.delinea[2*which(pos.qtls==epis.ef[ep,1])]*vet.delinea[(2*which(pos.qtls==epis.ef[ep,2]))+1]
								if (epis.ef[ep,3]==2) vet.delinea.epist[ep]<-vet.delinea[(2*which(pos.qtls==epis.ef[ep,1]))+1]*vet.delinea[2*which(pos.qtls==epis.ef[ep,2])]}}
						logdens[l]<-dnorm(dados[j,1],c(vet.delinea,vet.delinea.epist)%*%rbind(vet.coef,vet.coef.epist),sqrt(sigma2.vig),log=TRUE)}
					prob<-exp(probQTL+logdens-max(probQTL+logdens))/sum(exp(probQTL+logdens-max(probQTL+logdens)))
					mat.delinea[j,2*i]<-gen[rDiscreta(prob)]
					mat.delinea[j,(2*i)+1]<-0.5-abs(mat.delinea[j,2*i])
					if (length(esta_epist)>0){
						for(ep in esta_epist){
							if (epis.ef[ep,3]==0) mat.delinea.epist[j,ep]<-mat.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*mat.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]
							if (epis.ef[ep,3]==3) mat.delinea.epist[j,ep]<-mat.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*mat.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
							if (epis.ef[ep,3]==1) mat.delinea.epist[j,ep]<-mat.delinea[j,2*which(pos.qtls==epis.ef[ep,1])]*mat.delinea[j,(2*which(pos.qtls==epis.ef[ep,2]))+1]
							if (epis.ef[ep,3]==2) mat.delinea.epist[j,ep]<-mat.delinea[j,(2*which(pos.qtls==epis.ef[ep,1]))+1]*mat.delinea[j,2*which(pos.qtls==epis.ef[ep,2])]}}}}}
		#
		residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))
		#
		#### update the general mean
		#
		vet.coef[1,1]<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[1,1],mat.delinea[,1])[[1]]
		residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))
		#
		#### update QTLs' effects
		#
		if (num.QTLs>0){
			for (i in 1:num.QTLs){
				vet.coef[(2*i),1]<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[(2*i),1],mat.delinea[,(2*i)])[[1]]
				residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))
				vet.coef[(2*i)+1,1]<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef[(2*i)+1,1],mat.delinea[,(2*i)+1])[[1]]
				residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))}}
		#
		if (length(epis.ef)>0){
			for (i in 1:nrow(vet.coef.epist)){
				vet.coef.epist[i,1]<-poster_effect(media.prior,sigma2.prior,sigma2.vig,residuos,vet.coef.epist[i,1],mat.delinea.epist[,i])[[1]]
				residuos<-dados[,1]-(cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist))}}
		#
		#### update the error variance
		#
		sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
		#
		############################## record informations for posterior inference
		#
		if (int>burnin & int%%jumps==0){
			logveros<-round(sum(dnorm(dados[,1],mean=cbind(mat.delinea,mat.delinea.epist)%*%rbind(vet.coef,vet.coef.epist),sd=sqrt(sigma2.vig),log=TRUE)),4)
			cat('',logveros,file=paste(directory,"/log_likeli.txt",sep=""),append=T)	
			cat('',num.QTLs,file=paste(directory,"/number_QTLs.txt",sep=""),append=T)
			cat('',pos.qtls,file=paste(directory,"/QTLs_position.txt",sep=""),append=T)
			cat('',vet.coef,file=paste(directory,"/QTLs_effects.txt",sep=""),append=T)
			cat('',sigma2.vig,file=paste(directory,"/error_variance.txt",sep=""),append=T)
			cat('',vet.coef.epist,file=paste(directory,"/epist_effects.txt",sep=""),append=T)
			cat('',epis.ef,file=paste(directory,"/epist_type.txt",sep=""),append=T)
			cat('',length(epis.ef),file=paste(directory,"/number_epist.txt",sep=""),append=T)}}
}
#
##################
##################
# Running a toy example
##################
##################
#
library(Rcpp)
library(compiler)
#
directory<-"/Users/Daiane/Downloads" # include a directory in your computer where the data files are recorded
dados<-matrix(scan(file=paste(directory,"/dados_epistasia.txt",sep="")),nrow=300)
mapa<-read.table(paste(directory,"/mapa_epistasia.txt",sep=""),h=F)
#
loc.marc<-mapa[,2]
const<-0
for (i in 2:nrow(mapa)){
  if (mapa[i,3]>mapa[i-1,3]) const<-loc.marc[i-1]+100
  loc.marc[i]<-loc.marc[i]+const}
#
loc.marc<-loc.marc/100
gen.marc<-dados[,-1]
sai<-c(6,10,13,20,21,22,27,30,38,40,50,58,65,66,71,75,77,79,80,93,98,99,100,102) ### perfect association with flanker markers
loc.marc<-loc.marc[-sai]
gen.marc<-gen.marc[,-sai]
#
set.seed(100)
QTLs<-sort(sample(1:ncol(gen.marc),10,replace=F))
loc.marc[QTLs]
#
set.seed(10)
mi.verd<-20
sigma2.verd<-0.25
ef.adi.verd<-c(-0.60,0.90,0.25,-0.40,0.40)
ef.dom.verd<-c(0.30,0.05,-0.25,0.15,-0.15)
ef.epis.verd<-c(0.02,-0.19,-0.50,0.70,-1.08)
#
gen.QTL<-gen.marc[,QTLs]
res<-rnorm(nrow(gen.marc),0,sqrt(sigma2.verd))
fenot<-mi.verd+
  ef.adi.verd[1]*gen.QTL[,3]+ef.dom.verd[1]*(0.5-abs(gen.QTL[,3]))+
  ef.adi.verd[2]*gen.QTL[,5]+ef.dom.verd[2]*(0.5-abs(gen.QTL[,5]))+
  ef.adi.verd[3]*gen.QTL[,7]+ef.dom.verd[3]*(0.5-abs(gen.QTL[,7]))+
  ef.adi.verd[4]*gen.QTL[,8]+ef.dom.verd[4]*(0.5-abs(gen.QTL[,8]))+
  ef.adi.verd[5]*gen.QTL[,9]+ef.dom.verd[5]*(0.5-abs(gen.QTL[,9]))+
  ef.epis.verd[1]*(gen.QTL[,2]*(0.5-abs(gen.QTL[,6])))+
  ef.epis.verd[2]*((0.5-abs(gen.QTL[,2]))*(0.5-abs(gen.QTL[,5])))+
  ef.epis.verd[3]*(gen.QTL[,4]*gen.QTL[,10])+
  ef.epis.verd[4]*((0.5-abs(gen.QTL[,1]))*gen.QTL[,5])+
  ef.epis.verd[5]*(gen.QTL[,5]*(0.5-abs(gen.QTL[,9])))+
  res
#
dados<-cbind(fenot,gen.marc[,-QTLs])
loc.marc<-loc.marc[-QTLs]
dados<-dados[,c(1,2,10,18,26,34,42,50,58,66,71)] # just for the example to run faster
loc.marc<-loc.marc[c(1,9,17,25,33,41,49,57,65,70)]
#
sourceCpp(paste(directory,"/cppKruskall.cpp",sep=""))
#set.seed(100)
DDRJ_mainepistQTL_nondensemap(dados=dados,loc.marc=loc.marc,directory=directory)
#
#### Analysing the results
#
log.likel<-scan(file=paste(directory,"/log_likeli.txt",sep=""))
num.qtls<-scan(file=paste(directory,"/number_QTLs.txt",sep=""))
pos.qtls<-scan(file=paste(directory,"/QTLs_position.txt",sep=""))
coef<-scan(file=paste(directory,"/QTLs_effects.txt",sep=""))
sigma2<-scan(file=paste(directory,"/error_variance.txt",sep=""))
num.epis<-scan(file=paste(directory,"/number_epist.txt",sep=""))/3
type.epis<-scan(file=paste(directory,"/epist_type.txt",sep=""))
coef.epis<-scan(file=paste(directory,"/epist_effects.txt",sep=""))
#
finsample<-length(num.qtls)
plot(log.likel,type='l')
library(coda)
seq<-mcmc(log.likel)
geweke.diag(seq)
#
mat.loc<-matrix(0,nrow=finsample,ncol=max(num.qtls))
mat.coef<-matrix(0,nrow=finsample,ncol=(max(num.qtls)*2+1)) # the first column of this matrix has the MCMC sample for the general mean, the second column has simulated values for the additive effect of the first QTL, the third column has simulated values for the dominance effect of the first QTL, the fourth column has simulated values for the additive effect of the second QTL, the fifth column has simulated values for the dominance effect of the second QTL and so on.
mat.coef.epis<-matrix(0,nrow=finsample,ncol=(max(num.epis))) # the first column of this matrix has the MCMC sample for the first epistasis effect, the second column has the effect of second epistais and so on.
mat.tipo.epis<-matrix(0,nrow=finsample,ncol=(3*max(num.epis))) # the first G columns of this matrix has the location of one QTL involved in each of the G epistasis of the current model and the second G columns has the loation of the second QTL involved in each of the G epistasis of the current model. The following G columns specify the type of each of the G epistasis, where 0 means AA epistasis, 1 means AD, 2 means DA e 3 means DD epistasis
obs<-1
obs2<-1
obs3<-1
obs4<-1
for (i in 1:finsample){
	if (num.qtls[i]>0){
		mat.loc[i,1:num.qtls[i]]<-pos.qtls[obs:(obs+num.qtls[i]-1)]
		mat.coef[i,1:(num.qtls[i]*2+1)]<-round(coef[obs2:(obs2+(num.qtls[i]*2+1)-1)],3)}
	if (num.epis[i]>0) mat.coef.epis[i,1:num.epis[i]]<-round(coef.epis[obs3:(obs3+(num.epis[i]-1))],3)
	if (num.epis[i]>0) mat.tipo.epis[i,1:(3*num.epis[i])]<-round(type.epis[obs4:(obs4+(3*num.epis[i])-1)],3)
	obs<-obs+num.qtls[i]
	obs2<-obs2+(num.qtls[i]*2+1)
	obs3<-obs3+num.epis[i]
	obs4<-obs4+(3*num.epis[i])}
#
epistasis<-NULL # in this matrix we have informations (QTLs' position and type) about each identified epistasis
for (i in 1:nrow(mat.tipo.epis)){
	for (j in 1:num.epis[i]) epistasis<-rbind(epistasis,mat.tipo.epis[i,c(j,num.epis[i]+j,2*num.epis[i]+j)])}
epistasis<-epistasis[-which(epistasis[,1]==0 & epistasis[,2]==0),]
#
plot(density(pos.qtls))
plot(pos.qtls)
plot(mat.loc[mat.loc[,1]!=0,1],type='l')
plot(mat.loc[mat.loc[,2]!=0,2],type='l')
plot(mat.loc[mat.loc[,3]!=0,3],type='l')
