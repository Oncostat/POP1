writeStanCRModel<-function(Y,nY,nlvlY,X,T=NULL,lvlSlpX=c(),lvlSlpT=c(),sub=FALSE,FEpriors=c("HS","Cauchy")){
	match.arg(FEpriors)
	if(nY!=length(nlvlY)) stop("The number of level for each Y should be provided")
	if(nY!=length(lvlSlpX)) stop("The number of X slope parameters for each Y should be provided")
	if(!is.null(T)) if(nY!=length(lvlSlpT)) stop("The number of T slope parameters for each Y should be provided")

	####################### 
	#define coef if not provided
	maxCateg<-nlvlY-1
	if(!length(lvlSlpX)){
		lvlSlpX<-list()
		for(i in 1:nY){
			lvlSlpX<-c(lvlSlpX,list(1:maxCateg[i]))
		}
		warnings("Without levels provided for the effects of X (lvlSlpX), the same effect was assigned to all levels (proportional odds assumption)")
	}
	if(!is.null(T)){
		if(!length(lvlSlpT)){
			lvlSlpT<-list()
			for(i in 1:nY){
				lvlSlpT<-c(lvlSlpT,list(1:maxCateg[i]))
			}
			warnings("Without levels provided for the effects of T (lvlSlpT), the same effect was assigned to all levels (proportional odds assumption)")
		}
	}else{
		lvlSlpT<-list()
		for(i in 1:nY){
			lvlSlpT<-c(lvlSlpT,list(rep(0,maxCateg[i])))
		}
	}



	#######################
	#assign default number at parameter to ease writing of likelihood
	#first parameter will be all intercepts (for all Y), then slpX, then slpC
	labLvlInt<-lapply(nlvlY-1,function(l) 1:l)
	labLvlSlpX<-lvlSlpX
	if(length(lvlSlpT)) labLvlSlpT<-lvlSlpT

	lvlInt<-lapply(labLvlInt,function(l) 1:length(l))
	lvlSlpX<-lapply(lvlSlpX,function(l){
		if(sum(l==0)){
			return(NULL)
		}else{
			return(1:length(l))
		}
	})
	if(length(lvlSlpT)){
		lvlSlpT<-lapply(lvlSlpT,function(l){
			if(sum(l==0)){
				return(NULL)
			}else{
				return(1:length(l))
			}
		})
	}


	labs<-list()
	for(i in 1:length(lvlInt)){labs<-c(labs,list(list(lvlInt[[i]],lvlSlpX[[i]],lvlSlpT[[i]])))}

	iter<-1
	for(i in 1:length(labs)){
		for(j in 1:length(labs[[i]])){
			if(!is.null(labs[[i]][[j]])){
				for(k in 1:length(labs[[i]][[j]])){
					labs[[i]][[j]][[k]]<-iter
					iter<-iter+1
				}
			}
		}
	}

	params<-labs
	for(i in 1:length(labs)){
		for(j in 1:length(labs[[i]])){
			tmp<-list()
			if(j==1){
				for(k in 1:maxCateg[i]){
					if(sum(labLvlInt[[i]]!=0)){
						tmp<-c(tmp,list(labs[[i]][[j]][labLvlInt[[i]]<=k]))
					}
				}
			}else{
				if(j==2){
					if(sum(labLvlSlpX[[i]]!=0)){
						for(k in 1:maxCateg[i]){tmp<-c(tmp,list(labs[[i]][[j]][labLvlSlpX[[i]]<=k]))}
					}
				}else{
					if(sum(labLvlSlpT[[i]]!=0)){
						for(k in 1:maxCateg[i]){tmp<-c(tmp,list(labs[[i]][[j]][labLvlSlpT[[i]]<=k]))}
					}
				}
			}
			params[[i]][[j]]<-tmp
		}
	}



	####################### 
	#Begin to write
	func<-paste0("functions{",
		"\n\treal crm_lpdf(matrix Y,",paste("vector",c("X",unlist(ifelse(is.null(T),list(NULL),"T")),"betaFx"),collapse=","),", int N, int J,int P,matrix rdmInt){",
			paste0("\n\t\treal eta",1:max(unlist(maxCateg)),";",collapse=""),
			paste0("\n\t\treal p",0:max(unlist(maxCateg)),";",collapse=""),
			"\n\t\treal out;",
			"\n\t\tvector[N] prob;",
			"\n\t\tfor(i in 1:N){",
			"\n\t\t\tprob[i]=0;",
			paste(sapply(1:nY,function(k){
					paste0(
						paste0("\n\t\t\teta1=betaFx[",params[[k]][[1]][[1]],"]+rdmInt[i,",k,"]",
							ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[1]]),paste0("+(",paste0("betaFx[",params[[k]][[2]][[1]],"]"),")*X[i]"),""),""),
							ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[1]]),paste0("+(",paste0("betaFx[",params[[k]][[3]][[1]],"]"),")*T[i]"),""),""),
						";"),
						paste0("\n\t\t\teta2=",paste(paste0("betaFx[",params[[k]][[1]][[2]],"]"),collapse = "+"),"+rdmInt[i,",k,"]",
							ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[2]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[2]][[2]],"]"),collapse="+"),")*X[i]"),""),""),
							ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[2]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[3]][[2]],"]"),collapse="+"),")*T[i]"),""),""),
						";"),
						paste0("\n\t\t\teta3=",paste(paste0("betaFx[",params[[k]][[1]][[3]],"]"),collapse = "+"),"+rdmInt[i,",k,"]",
							ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[3]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[2]][[3]],"]"),collapse="+"),")*X[i]"),""),""),
							ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[3]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[3]][[3]],"]"),collapse="+"),")*T[i]"),""),""),
						";"),
						"\n\t\t\tp0=1-inv_logit(eta1);\n\t\t\tp1=1-inv_logit(eta2)-p0;\n\t\t\tp2=1-inv_logit(eta3)-p0-p1;\n\t\t\tp3=1-p0-p1-p2;",
						paste0("if(Y[i,",k,"]==0){"),
							"\n\t\t\t\tprob[i]+=log(p0);",
						"\n\t\t\t}else{",
							paste0("\n\t\t\t\tif(Y[i,",k,"]==1){"),
								"\n\t\t\t\t\tprob[i]+=log(p1);",
							"\n\t\t\t\t}else{",
								paste0("\n\t\t\t\t\tif(Y[i,",k,"]==2){"),
									"\n\t\t\t\t\t\tprob[i]+=log(p2);",
								"\n\t\t\t\t\t}else{",
									paste0("\n\t\t\t\t\t\tif(Y[i,",k,"]==3){"),
										"\n\t\t\t\t\t\t\tprob[i]+=log(p3);",
									"\n\t\t\t\t\t\t}",
								"\n\t\t\t\t\t}",
							"\n\t\t\t\t}",
						"\n\t\t\t}"
					)
				}),
			collapse=""),
			"\n\t\t}",      
			"\n\t\tout=sum(prob);",
			"\n\t\treturn(out);",
			"\n\t}\n",
			"\n\tvector vecLik(matrix Y,",paste("vector",c("X",unlist(ifelse(is.null(T),list(NULL),"T")),"betaFx"),collapse=","),", int N, int J,int P,matrix rdmInt){",
				paste0("\n\t\treal eta",1:max(unlist(maxCateg)),";",collapse=""),
				paste0("\n\t\treal p",0:max(unlist(maxCateg)),";",collapse=""),
				"\n\t\treal out;",
				"\n\t\tvector[N] prob;",
				"\n\t\tfor(i in 1:N){",
				"\n\t\t\tprob[i]=0;",
				paste(sapply(1:nY,function(k){
						paste0(
							paste0("\n\t\t\teta1=betaFx[",params[[k]][[1]][[1]],"]+rdmInt[i,",k,"]",
								ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[1]]),paste0("+(",paste0("betaFx[",params[[k]][[2]][[1]],"]"),")*X[i]"),""),""),
								ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[1]]),paste0("+(",paste0("betaFx[",params[[k]][[3]][[1]],"]"),")*T[i]"),""),""),
							";"),
							paste0("\n\t\t\teta2=",paste(paste0("betaFx[",params[[k]][[1]][[2]],"]"),collapse = "+"),"+rdmInt[i,",k,"]",
								ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[2]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[2]][[2]],"]"),collapse="+"),")*X[i]"),""),""),
								ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[2]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[3]][[2]],"]"),collapse="+"),")*T[i]"),""),""),
							";"),
							paste0("\n\t\t\teta3=",paste(paste0("betaFx[",params[[k]][[1]][[3]],"]"),collapse = "+"),"+rdmInt[i,",k,"]",
								ifelse(length(params[[k]][[2]]),ifelse(length(params[[k]][[2]][[3]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[2]][[3]],"]"),collapse="+"),")*X[i]"),""),""),
								ifelse(length(params[[k]][[3]]),ifelse(length(params[[k]][[3]][[3]]),paste0("+(",paste(paste0("betaFx[",params[[k]][[3]][[3]],"]"),collapse="+"),")*T[i]"),""),""),
							";"),
							"\n\t\t\tp0=1-inv_logit(eta1);\n\t\t\tp1=1-inv_logit(eta2)-p0;\n\t\t\tp2=1-inv_logit(eta3)-p0-p1;\n\t\t\tp3=1-p0-p1-p2;",
							paste0("if(Y[i,",k,"]==0){"),
								"\n\t\t\t\tprob[i]+=log(p0);",
							"\n\t\t\t}else{",
								paste0("\n\t\t\t\tif(Y[i,",k,"]==1){"),
									"\n\t\t\t\t\tprob[i]+=log(p1);",
								"\n\t\t\t\t}else{",
									paste0("\n\t\t\t\t\tif(Y[i,",k,"]==2){"),
										"\n\t\t\t\t\t\tprob[i]+=log(p2);",
									"\n\t\t\t\t\t}else{",
										paste0("\n\t\t\t\t\t\tif(Y[i,",k,"]==3){"),
											"\n\t\t\t\t\t\t\tprob[i]+=log(p3);",
										"\n\t\t\t\t\t\t}",
									"\n\t\t\t\t\t}",
								"\n\t\t\t\t}",
							"\n\t\t\t}"
						)
					}),
				collapse=""),
				"\n\t\t}",      
				"\n\t\treturn(prob);",
			"\n\t}",
			"\n}\n\n"
		)

	data<-paste0("data {",
					"\n\tint<lower=0> N;",
					"\n\tint<lower=0> J;",
					"\n\tint<lower=0> P;",
					"\n\tmatrix[N,J] Y;",
					"\n\tvector[N] X;",
					ifelse(is.null(T),"","\n\tvector[N] T;"),
					ifelse(sub,"","\n\tvector[J] Zero;"),
					"\n\tint<lower=0> nRdmLabs;\n\tint rdmLabs[N];",
					"\n}\n\n")

	params<-paste0("parameters {",
		"\n\tvector[P] z;",
		ifelse(FEpriors=="HS",paste0(
			"\n\treal<lower=0> r1_global;",
			"\n\treal<lower=0> r2_global;",
			"\n\tvector<lower=0>[P] r1_local;",
			"\n\tvector<lower=0>[P] r2_local;"),
			"\n\tvector[P] betaFx;"),
		"\n\tcholesky_factor_corr[J] Lcorr;",
		"\n\tvector<lower=0>[J] sigma_rdmInt;",
		"\n\tmatrix[nRdmLabs,J] rdmInt;",
		"\n}\n\n")

	transParams<-paste0("transformed parameters{",
		"\n\tmatrix[N,J] rdmIntId;",
		ifelse(FEpriors=="HS",paste0(
			"\n\tvector[P] betaFx;",
			"\n\treal<lower=0> tau;",
			"\n\tvector<lower=0>[P] lambda;",
			"\n\ttau=r1_global .* sqrt(r2_global);",
			"\n\tlambda=r1_local .* sqrt(r2_local);",
			"\n\tbetaFx=z .* lambda*tau;"
			),""),
		"\n\tfor(i in 1:N){rdmIntId[i,]=rdmInt[rdmLabs[i],];}",
		"\n}\n\n",collapse="")

	model<-paste0("model {",
		"\n\tz~normal(0,1);",
		ifelse(FEpriors=="HS",
			paste0("\n\tr1_local~normal(0,1);",
				"\n\tr2_local~inv_gamma(0.5,0.5);",
				"\n\tr1_global~normal (0,1);",
				"\n\tr2_global~inv_gamma(0.5,0.5);"),
			"\n\tfor(i in 1:P){betaFx[i] ~ cauchy(0,2.5);}"),
		ifelse(sub,"\n\tsigma_rdmInt ~ cauchy(0, 5);\n\tfor(j in 1:nRdmLabs){\n\t\trdmInt[j] ~ normal(0,sigma_rdmInt);\n\t}",
			"\n\tsigma_rdmInt ~ cauchy(0, 5);\n\tLcorr ~ lkj_corr_cholesky(1);\n\tfor(j in 1:nRdmLabs){rdmInt[j,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_rdmInt, Lcorr));}"),
		ifelse(is.null(T),"\n\ttarget+=crm_lpdf(Y|X,betaFx,N,J,P,rdmIntId);","\n\ttarget+=crm_lpdf(Y|X,T,betaFx,N,J,P,rdmIntId);"),
		"\n}\n\n")

	genQ<-paste0("generated quantities {",
		"\n\tvector[N] ll;",
		"\n\tmatrix[J,J] Omega;",
		"\n\tmatrix[J,J] Sigma;",
		ifelse(is.null(T),"\n\tll=vecLik(Y,X,betaFx,N,J,P,rdmIntId);","\n\tll=vecLik(Y,X,T,betaFx,N,J,P,rdmIntId);"),
		"\n\tOmega = multiply_lower_tri_self_transpose(Lcorr);",
		"\n\tSigma = quad_form_diag(Omega, sigma_rdmInt); ",
		"\n}\n\n")

	if(sub){
		return(list(stanFile=c(func,data,params,transParams,model)))
	}else{
		return(list(stanFile=c(func,data,params,transParams,model,genQ)))
	}
}


