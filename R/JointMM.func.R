JointMM.func <-
function(data,cov.name.long,cov.name.surv,threshold.zero.prop = 0.1,is.longi.model.only = FALSE, quad.n=10){
	if(!is.data.frame(data)){
		print("please organize `data` in a data frame")
		return (1)
	}
    ind.match.arguments <- sum(is.na(match(c('id','Y','obstime','time','disease','start','stop','event'),colnames(data))))
	if(ind.match.arguments >=1){
	print("please make sure that 1) the subject ID in the data frame is named as `id`, 2) the relative abundance of the taxon from the microbiome sample is named as `Y`, 3) the observed time of the corresponding microbiome sample is named as `obstime`, 4) the time to disease onset is named as `time`, the event onset indicator is named as `disease`, 5) the starting time point of the microbiome sample represented is named as `start` which is equal to `obstime`, 6) the stopping time point of the microbiome sample represented is named as `stop` which is the end of the time period, 7) the indicator of whether the event happened during this time period is named as `event`, and re-run the function")
	return (0)
	}
	ind.match.cov.names <- sum(is.na(match(c(cov.name.long,cov.name.surv),colnames(data))))
    
	if(ind.match.cov.names >=1){
	print("please make sure that the names of the covariates in the longitudinal and survival sub-models are in the data frame")
	return (1)
	}
	dat = data
  dat = dat[order(dat$id,dat$obstime),]
  dat.id = dat[!duplicated(dat$id),] 
  X = dat[,cov.name.long]
  Y = dat[,'Y']
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if(!is.longi.model.only){
  dat.surv = dat.id[,c('time','disease')]
  Z = as.matrix(dat.id[,cov.name.surv])
  }
  
  X.aug <- cbind(intersept = 1, X)
  nsample = length(unique(dat$id))
  N.obs <- sapply(as.character(unique(dat$id)),function(x)sum(dat$id==x)) 
  N = sum(N.obs)
  tmp = unlist(sapply(N.obs,function(x){c(rep(1,x),rep(0,N))}))
  prod.mat = matrix(tmp[1:(nsample*N)],byrow = TRUE,nrow = nsample,ncol = N)
  gherm <- generate_gaussian_quad_points(quad.n = quad.n)
  gh.weights <- matrix(rep(gherm$weights,nsample),nrow = nsample,byrow = TRUE) 
  gh.nodes <- matrix(rep(gherm$nodes,N),nrow = N,byrow = TRUE)
  bivar.gherm = generate_gaussian_quad_points_surv(quad.n = quad.n)
  X.test.coeff.index <- rep(FALSE,ncol(X.aug))
  if(!is.longi.model.only){
  Z.test.coeff.index <-  rep(FALSE,2) 
  }
  if(!is.longi.model.only){
	d = list(dat = dat,dat.surv = dat.surv, X = X,X.aug = X.aug, Y= Y, Z = Z, nsample =nsample,N.obs = N.obs, N = N,prod.mat = prod.mat,gherm = gherm,gh.weights = gh.weights,gh.nodes = gh.nodes,bivar.gherm = bivar.gherm,X.test.coeff.index = X.test.coeff.index,Z.test.coeff.index =Z.test.coeff.index) 
  }else{
	d = list(dat = dat,X = X,X.aug = X.aug, Y= Y,nsample =nsample,N.obs = N.obs, N = N,prod.mat = prod.mat,gherm = gherm,gh.weights = gh.weights,gh.nodes = gh.nodes,bivar.gherm = bivar.gherm,X.test.coeff.index = X.test.coeff.index)
  }
	zero.prop = mean(data$Y==0)
	reduced.model.indicator = zero.prop < threshold.zero.prop
   if(!reduced.model.indicator){
		if(!is.longi.model.only){
			d = fit.JointMM.fullmodel (d,quad.n=quad.n)
		}else{
			d = fit.JointMM.fullmodel.longitudinal.part(d,quad.n=quad.n)
		}
   }else{
		n.zero.prop = sum(d$Y==0)
		pseudo.prop = sort(d$Y[d$Y>0])[1:5]
		d$Y[d$Y==0] = sample(pseudo.prop,size = n.zero.prop,replace = TRUE)
		if(!is.longi.model.only){
			d = fit.JointMM.reducedmodel (d,quad.n=quad.n)
		}else{
			d = fit.JointMM.reducedmodel.longitudinal.part (d,quad.n=quad.n)
		}

  }
  if(!is.longi.model.only){
	res = d[c("par.est.Wald","SEs","est.hessian","Wald.Ts","pvals.Wald","status.Wald")]  
  }else{
	res = d[c("par.est.Wald","SEs","est.hessian","Wald.Ts.Longonly","pvals.Wald.Longonly","status.Wald")]
  }
  
  return(res)
}
