

library(parallel)
library(rstan)
library(ordinal)
options(mc.cores = 4)
source("src/stanDMfcts.R")
source("src/writeStanCRModel.R")
load("results/saveInits.RData")




#get posterior values of the univariate models as initial values for the multivariate model
date()
inits<-getInits(fit1)
nam<-names(inits)
for(i in 2:5){
	tmp<-getInits(get(i))
	for(j in nam){
		if(j=="rdmInt"){
			inits[[j]]<-cbind(inits[[j]],tmp[[j]])
		}else{
			inits[[j]]<-c(inits[[j]],tmp[[j]])
		}
	}
}
initBetaFx<-inits[["betaFx"]]
inits$tau<-mean(inits$tau)
inits$r1_global<-mean(inits$r1_global)
inits$r2_global<-1


#Create data for Stan fit
stan_CR_data<-createStanData(data[,nams],data[,dose],data[,cycle],data[,patId],length(initBetaFx))
stan_CR_data$Zero<-rep(0,length(grep("grade",colnames(data))))


#Write code
code<-writeStanCRModel(Y="Y",nY=ncol(stan_CR_data$Y),nlvlY=rep(4,ncol(stan_CR_data$Y)),X="X",T="T",lvlSlpX=list(1:3,1:3,1:3,1:3,1:3),lvlSlpT=list(1:3,1:3,1:3,1:3,1:3),sub=FALSE)
outFile<-file("stanCodes/CRModel.stan")
writeLines(code$stanFile, outFile,sep="")
close(outFile)

set.seed(123456)
fitCR<-stan(file="stanCodes/CRModel.stan",data=stan_CR_data,warmup=1000,iter=5000,chain=4,init=list(inits,inits,inits,inits),control = list(adapt_delta = 0.99))

rm(list=ls()[!ls()%in%c(fitCR,stan_CR_data,inits)])

save.image("results/CRModel.RData")

