

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

#keep all intercepts, all 1st coef for dose and all 1st coef for cycle
keep<-c(1:3,4,7, #1st type of tox
	    10:12,13,16, #2nd type of tox
	    19:21,22,25, #...
	    28:30,31,34,
	    37:39,40,43) 
initBetaFx<-initBetaFx[keep]
inits$betaFx<-inits$betaFx[keep]
inits$z<-inits$z[keep]
inits$r1_local<-inits$r1_local[keep]
inits$r2_local<-inits$r2_local[keep]
inits$lambda<-inits$lambda[keep]

#Create data for Stan fit
stan_PO_data<-createStanData(data[,nams],data[,dose],data[,cycle],data[,patId],length(initBetaFx))
stan_PO_data$Zero<-rep(0,length(grep("grade",colnames(data))))


#Write code
code<-writeStanCRModel(Y="Y",nY=ncol(stan_PO_data$Y),nlvlY=rep(4,ncol(stan_PO_data$Y)),X="X",T="T",lvlSlpX=list(1,1,1,1,1),lvlSlpT=list(1,1,1,1,1),sub=FALSE)
outFile<-file("stanCodes/POModel.stan")
writeLines(code$stanFile, outFile,sep="")
close(outFile)

set.seed(123456)
fitPO<-stan(file="stanCodes/POModel.stan",data=stan_PO_data,warmup=1000,iter=5000,chain=4,init=list(inits,inits,inits,inits),control = list(adapt_delta = 0.99))

rm(list=ls()[!ls()%in%c(fitPO,stan_PO_data,inits)])

save.image("results/POModel.RData")

