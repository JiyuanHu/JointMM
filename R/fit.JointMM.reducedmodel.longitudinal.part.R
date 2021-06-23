fit.JointMM.reducedmodel.longitudinal.part <-
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
  eps = function(maxY){
    ifelse(maxY>0.01,0.01,ifelse(maxY>1e-3,1e-3,1e-4))
  }
  q <- max(Y) +eps (max(Y))
  opt.H1.Sbeta <- optim(
	par =c(.5,5,rep(0,sum(!X.test.coeff.index))), ## s2,phi,q,beta
    fn=loglik.B.present.part.only,
	method = 'L-BFGS-B',
    lower = c(rep(1e-10,2),
              rep(-Inf,sum(!X.test.coeff.index))),
    upper = c(rep(Inf,2),rep(Inf,sum(!X.test.coeff.index))),
    control= list(maxit = 1e3),
	hessian = TRUE,
	X.aug=X.aug,Y=Y,prod.mat=prod.mat,
    X.test.coeff.index = X.test.coeff.index,
    gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n)
  est.H1 = opt.H1.Sbeta$par
  names(est.H1) = c('se2','phi',paste0('beta',1:ncol(X.aug)-1))
  hes2 = opt.H1.Sbeta$hessian
  
	status = c(opt.H1.Sbeta$convergence)
	I.inverse = solve(hes2 +1e-8)
	vars = diag(I.inverse)
	vars[which(vars<0)] = -vars[which(vars<0)]
	SEs = sqrt(vars)    
 
 	Wald.proportion.Longonly = pval.proportion.Longonly = rep(NA,ncol(X))
	for (j in 1:ncol(X)){
		ind = 3+j
		Wald.proportion.Longonly[j] = ifelse(is.na(SEs[ind]),NA,(est.H1[ind]/SEs[ind])^2)
		pval.proportion.Longonly[j]= pchisq(Wald.proportion.Longonly[j],df =1,lower.tail = FALSE)
  }
 
 names(Wald.proportion.Longonly) = names(pval.proportion.Longonly) = colnames(X)

  d$par.est.Longonly = est.H1
  d$SEs = SEs
  d$est.hessian = hes2
  d$I.inverse = I.inverse

  d$Wald.Ts.Longonly = Wald.proportion.Longonly
  d$pvals.Wald.Longonly = pval.proportion.Longonly
  d$status.Wald.Longonly = status
  return(d)
}
