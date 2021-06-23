fit.JointMM.fullmodel.longitudinal.part.cross.part.corr <-
function(d,quad.n){ 
  dat = d$dat
  dat.surv = d$dat.surv
  X = d$X
  X.aug = d$X.aug
  Y= d$Y
  nsample =d$nsample
  N.obs = d$N.obs
  N = d$N
  prod.mat = d$prod.mat
  gherm = d$gherm
  gh.weights = d$gh.weights
  gh.nodes = d$gh.nodes
  bivar.gherm = d$bivar.gherm
  X.test.coeff.index = d$X.test.coeff.index
##
    start2 = c(.5,rep(0,sum(!X.test.coeff.index)),.5,20,1-1e-3,rep(0,sum(!X.test.coeff.index)),0) 
 names(start2) = c(c('se1',paste0('alpha',1:ncol(X.aug))),
                      'se2','phi','q','beta1',c(paste0('beta',2:ncol(X.aug))),'delta0')
  lower = c(c(1e-5,rep(-Inf,sum(!X.test.coeff.index))),
            rep(1e-3,2),max(Y)+1e-3,c(rep(-Inf,sum(!X.test.coeff.index)-1),-Inf,-Inf)
			)
  upper = c(rep(Inf,sum(!X.test.coeff.index)+1),
            rep(Inf,2),1,c(rep(Inf,sum(!X.test.coeff.index)-1),Inf,Inf)
			)
  opt.H1 <- optim(
    par=start2, 
    fn=joint.loglike.cross.part.corr,
    gr = grad.joint.loglike.cross.part.corr,
    method = 'L-BFGS-B',
    control= list(maxit = 1e3),
	lower = lower,
    upper = upper,
    X.aug=X.aug,Y=Y,prod.mat=prod.mat,
    X.test.coeff.index=X.test.coeff.index,
    Z=Z,dat.surv=dat.surv, 
    Z.test.coeff.index=Z.test.coeff.index,
    nsample = nsample,N= N,
    quad.n=quad.n,
    bivar.gherm=bivar.gherm,
    hessian = TRUE
  )
  est.H1 = opt.H1$par
  est.hessian = opt.H1$hessian
  status = opt.H1$convergence
  ##
  	  est.H1 = c(opt.H1$par)
	  
	  est.hessian = opt.H1$hessian
	  status = opt.H1$convergence
	  I.inverse = solve(est.hessian +1e-8)
	vars = diag(I.inverse)
	vars[which(vars<0)] = -vars[which(vars<0)]
	SEs = sqrt(vars)     
	Wald.all.Longonly = Wald.presence.Longonly = Wald.proportion.Longonly = pval.all.Longonly = pval.presence.Longonly = pval.proportion.Longonly = rep(NA,ncol(X))
	for (j in 1:ncol(X)){ 
		ind = c(2+j,ncol(X)+6+j)
		Wald.all.Longonly[j] = ifelse(sum(is.na(I.inverse[ind,ind])),NA,t(est.H1[ind]) %*% solve(I.inverse[ind,ind]) %*% est.H1[ind])
	   
		Wald.presence.Longonly[j] = ifelse(is.na(SEs[ind[1]]),NA,(est.H1[ind[1]]/SEs[ind[1]])^2)
		Wald.proportion.Longonly[j] = ifelse(is.na(SEs[ind[2]]),NA,(est.H1[ind[2]]/SEs[ind[2]])^2)
		
		pval.all.Longonly[j]= pchisq(Wald.all.Longonly[j],df =2,lower.tail = FALSE)
		pval.presence.Longonly[j] = pchisq(Wald.presence.Longonly[j],df =1,lower.tail = FALSE)
		pval.proportion.Longonly[j]= pchisq(Wald.proportion.Longonly[j],df =1,lower.tail = FALSE)
    
  }
  Wald.Ts.Longonly = rbind(Wald.all.Longonly,Wald.presence.Longonly,Wald.proportion.Longonly)
  
  pvals.Longonly = rbind(pval.all.Longonly,pval.presence.Longonly,pval.proportion.Longonly)
  colnames(Wald.Ts.Longonly) = colnames(pvals.Longonly) = colnames(X)

  d$par.est.Longonly = est.H1
  d$SEs = SEs
  names(d$par.est.Longonly) = names(d$SEs)= names(start2)
  d$est.hessian = est.hessian
  d$I.inverse = I.inverse
  d$Wald.Ts.Longonly = Wald.Ts.Longonly
  d$pvals.Wald.Longonly = pvals.Longonly
  d$status.Wald.Longonly = status
	#wald test of the cross part correlation
  Wald.cross.part.corr = ifelse(is.na(SEs[length(SEs)]),NA,(est.H1[length(SEs)]/SEs[length(SEs)])^2)
  pval.cross.part.corr = pchisq(Wald.cross.part.corr,df =1,lower.tail = FALSE)
  d$Wald.cross.part.corr =Wald.cross.part.corr
  d$pval.cross.part.corr = pval.cross.part.corr
  return(d)
}
