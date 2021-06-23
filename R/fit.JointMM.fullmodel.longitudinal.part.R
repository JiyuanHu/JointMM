fit.JointMM.fullmodel.longitudinal.part <-
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
  
  opt.H1.logit <- optim(par=c(.5,rep(0,sum(!X.test.coeff.index))), 
                        fn=loglik.A,
                        method = 'L-BFGS-B',
                        lower = c(0.00001,
                                  rep(-Inf,sum(!X.test.coeff.index))),
                        upper = Inf,
                        X.aug=X.aug,Y=Y,prod.mat=prod.mat,
                        X.test.coeff.index = X.test.coeff.index,
                        gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                        control= list(maxit = 1e3),
                        hessian = TRUE
  )
  opt.H1.Sbeta <- optim(par = c(rep(0,sum(!X.test.coeff.index)-1),0.5,20,max(Y)+1e-4,-0.5),
                        fn=loglik.B.change.para.order,
                        gr = gradloglik.B.change.para.order,
                        method = 'L-BFGS-B',
                        lower = c(rep(-Inf,sum(!X.test.coeff.index)-1),rep(1e-3,2),max(Y)+1e-3,-Inf),
                        upper = c(rep(Inf,sum(!X.test.coeff.index)-1),Inf,Inf,1,Inf),
                        X.aug=X.aug,Y=Y,prod.mat=prod.mat,
                        X.test.coeff.index = X.test.coeff.index,
                        gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                        control= list(maxit = 1e3),
                        hessian = TRUE
  )
	  est.H1 = c(opt.H1.logit$par,opt.H1.Sbeta$par)
	  est.hessian = adiag(opt.H1.logit$hessian,opt.H1.Sbeta$hessian)
	  status = unique(c(opt.H1.logit$convergence,opt.H1.Sbeta$convergence))
	  I.inverse = solve(est.hessian +1e-8)
	vars = diag(I.inverse)
	vars[which(vars<0)] = -vars[which(vars<0)]
	SEs = sqrt(vars)     
	Wald.all.Longonly = Wald.presence.Longonly = Wald.proportion.Longonly = pval.all.Longonly = pval.presence.Longonly = pval.proportion.Longonly = rep(NA,ncol(X))
	for (j in 1:ncol(X)){ 
		ind = c(2+j,ncol(X)+2+j)
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
  names(d$par.est.Longonly) = names(d$SEs)= c(c('se1',paste0('alpha',1:(ncol(X)+1))),
                                          c(paste0('beta',2:(ncol(X)+1)),'se2','phi','q','beta1'))
  d$est.hessian = est.hessian
  d$I.inverse = I.inverse
  d$Wald.Ts.Longonly = Wald.Ts.Longonly
  d$pvals.Wald.Longonly = pvals.Longonly
  d$status.Wald.Longonly = status
  
  return(d)
}
