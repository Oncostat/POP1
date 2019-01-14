##################################################################################################################################
##################################################################################################################################
#
##################################################################################################################################
##################################################################################################################################





#load your data here

library(parallel)
library(rstan)
library(ordinal)
options(mc.cores = 4)
source("src/stanDMfcts.R")
source("src/writeStanCRModel.R")

#insert the names of the columns of the grade of the 5 types of toxicity
nams<-c("gradeC","gradeD","gradeH","gradeG","gradeO") #for c("Cutaneous","Digestive","General disorder","Hematologic","Others")
#Name of the columns of the individual indices
patId<-"patId" #labels should be standardized with values from 1 to the number of patient
#Name of the columns of the dose
dose<-"doseStd"
#Name of the columns of the cycle
cycle<-"cycle"



#Initialize from univariate fixed effect full continuation ratio model, frequentist version for speed
fm1<-clm(as.factor(data[,nams[1]])~1,nominal=~data[,dose]+data[,cycle])
fm2<-clm(as.factor(data[,nams[2]])~1,nominal=~data[,dose]+data[,cycle])
fm3<-clm(as.factor(data[,nams[3]])~1,nominal=~data[,dose]+data[,cycle])
fm4<-clm(as.factor(data[,nams[4]])~1,nominal=~data[,dose]+data[,cycle])
fm5<-clm(as.factor(data[,nams[5]])~1,nominal=~data[,dose]+data[,cycle])

#Write Stan code for univariate model
code<-writeStanCRModel(Y="Y",nY=1,nlvlY=4,X=c("X","T"),lvlSlpX=list(1:3),lvlSlpT=list(1:3),sub=TRUE)
outFile<-file("stanCodes/subModel.stan")
writeLines(code$stanFile,outFile,sep="")
close(outFile)



for(i in 1:5){
	#get initialization values
	initBetaFx<-stdInitBeta(coef(get(paste0("fm",i))))
	#Create data for Stan fit
	stan_data<-createStanData(data[,nams[i]],data[,dose],data[,cycle],data[,patId],length(initBetaFx))
	#fit
	assign(paste0("fit",i),stan(file="stanCodes/subModel.stan",data=stan_data,warmup=500,iter=1000,chain=1,init=list(list(betaFx=initBetaFx)),control = list(adapt_delta = 0.99)))
}

save.image("results/saveInits.RData")





