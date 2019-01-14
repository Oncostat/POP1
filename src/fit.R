colVars <- function(mat){
	n<-dim(mat)[[1]]
	v<-dim(mat)[[2]]
	return(.colMeans(((mat-matrix(.colMeans(mat,n,v),nrow=n,ncol=v,byrow=TRUE))^2),n,v)*n/(n-1))
}

WAIC<-function(fit){
	log_lik<-extract(fit,"ll")$ll
	dim(log_lik)<-if(length(dim(log_lik))==1) c(length(log_lik),1) else c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
	S<-nrow(log_lik)
	n<-ncol(log_lik)
	lpd<-log(colMeans(exp(log_lik)))
	p_waic<-colVars(log_lik)
	elpd_waic<-lpd-p_waic
	waic<- -2*elpd_waic
	pointwise<-cbind(waic,lpd,p_waic,elpd_waic)
	total<-colSums(pointwise)
	se<-sqrt(n*colVars(pointwise))
	return(list(waic=total["waic"],elpd_waic=total["elpd_waic"],p_waic=total["p_waic"],pointwise=pointwise,total=total,se=se))
}




get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


inv_logit<-function(x) exp(x)/(1+exp(x))




gen_quant_r<-function(Y,X,T,model,nSpls=1000,seed=123456){
	maxGrade<-max(c(Y))
	set.seed(seed)
	postBetas<-apply(beta_post,2,function(x){sample(x,size=nSpls)})
	postRdms<-array(0,dim=c(nrow(Y),nSpls,dim(rdm_post)[3]))
	for(k in 1:dim(rdm_post)[3]){
		for(n in 1:nrow(Y)){
			postRdms[n,,k]<-sample(rdm_post[,n,k],size=nSpls)
		}
	}
	preds<-array(0,dim=c(nrow(Y),maxGrade+1,nSpls,ncol(Y)))
	sim<-array(0,dim=c(nrow(Y),nSpls,ncol(Y)))
	for(i in 1:nrow(Y)){
		for(k in 1:ncol(Y)){
			if(model=="CR"){
				eta1=postBetas[,1+(k-1)*9]+postRdms[i,,k]+(postBetas[,4+(k-1)*9])*X[i]+postBetas[,7+(k-1)*9]*T[i];
				eta2=postBetas[,1+(k-1)*9]+postBetas[,2+(k-1)*9]+postRdms[i,,k]+(postBetas[,4+(k-1)*9]+postBetas[,5+(k-1)*9])*X[i]+(postBetas[,7+(k-1)*9]+postBetas[,8+(k-1)*9])*T[i];
				eta3=postBetas[,1+(k-1)*9]+postBetas[,2+(k-1)*9]+postBetas[,3+(k-1)*9]+postRdms[i,,k]+(postBetas[,4+(k-1)*9]+postBetas[,5+(k-1)*9]+postBetas[,6+(k-1)*9])*X[i]+(postBetas[,7+(k-1)*9]+postBetas[,8+(k-1)*9]+postBetas[,9+(k-1)*9])*T[i];
			}else{
				eta1<-postBetas[,1+(k-1)*5]+postRdms[i,,k]+(postBetas[,4+(k-1)*5])*X[i]+postBetas[,5+(k-1)*5]*T[i];
				eta2<-postBetas[,1+(k-1)*5]+postBetas[,2+(k-1)*5]+postRdms[i,,k]+(postBetas[,4+(k-1)*5])*X[i]+(postBetas[,5+(k-1)*5])*T[i];
				eta3<-postBetas[,1+(k-1)*5]+postBetas[,2+(k-1)*5]+postBetas[,3+(k-1)*5]+postRdms[i,,k]+(postBetas[,4+(k-1)*5])*X[i]+(postBetas[,5+(k-1)*5])*T[i];
			}
			p0<-1-inv_logit(eta1);
			p1<-1-inv_logit(eta2)-p0;
			p2<-1-inv_logit(eta3)-p0-p1;
			p3<-1-p0-p1-p2;
			preds[i,1,,k]<-p0
			preds[i,2,,k]<-p1
			preds[i,3,,k]<-p2
			preds[i,4,,k]<-p3
		}
	}
	for(k in 1:ncol(Y)){
		for(i in 1:nrow(Y)){
			sim[i,,k]<-do.call(rbind,lapply(1:nSpls,function(n){
					if(sum(preds[i,,n,k]<0)){
						return(NA)
					}else{
						which(rmultinom(1,1,preds[i,,n,k])==1)-1
					}
				}))
		}
	}
	out<-list(sim=sim,preds=preds)
	return(out)
}





predDist<-function(fit,stan_data,model=c("CR","PO")){
	match.arg(model)
	if(length(model)>1){
		model<-"CR"
		warnings("Too many arg for 'model', 'CR' was used as default choice")
	}
	ext_fit0<-extract(fit0)
	beta_post<-ext_fit0$betaFx
	rdm_post<-ext_fit0$rdmIntId
	yPred<-gen_quant_r(stan_data$Y,stan_data$X,stan_data$T,model,1000) #Simulate y based on new x
	gc()
	Y<-stan_data$Y
	X<-stan_data$X
	T<-stan_data$T
	obsProbs<-array(0,dim=c(max(T),max(c(Y))+1,ncol(Y)))
	expProbs<-array(0,dim=c(max(T),max(c(Y))+1,1000,ncol(Y)))
	for(k in 1:dim(obsProbs)[3]){
		for(grade in 0:(dim(obsProbs)[2]-1)){
			for(cyc in 1:dim(obsProbs)[1]){
				obsProbs[cyc,grade+1,k]<-sum(Y[T==cyc,k]>=grade)/sum(T==cyc)
				expProbs[cyc,grade+1,,k]<-sapply(1:1000,function(j){
					tmp<-yPred$sim[T==cyc,j,k]
					tmp<-tmp[!is.na(tmp)]
					sum(tmp>=grade)/sum(T[!is.na(tmp)]==cyc)
				})
			}
		}
	}
	return(tmp)
}






