fit.JointMM.reducedmodel <-
function(d,quad.n){
  dat = d$dat
  dat.surv = d$dat.surv
  X = d$X
  X.aug = d$X.aug
  Y= d$Y
  Z = d$Z
  nsample =d$nsample
  N.obs = d$N.obs
  N = d$N
  prod.mat = d$prod.mat
  gherm = d$gherm
  gh.weights = d$gh.weights
  gh.nodes = d$gh.nodes
  bivar.gherm = d$bivar.gherm
  X.test.coeff.index = d$X.test.coeff.index

  d$Z.test.coeff.index = FALSE
  Z.test.coeff.index = d$Z.test.coeff.index;
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
  par2 = opt.H1.Sbeta$par
  names(par2) = c('se2','phi',paste0('beta',1:ncol(X.aug)-1))
  hes2 = opt.H1.Sbeta$hessian
  se2.hat <- par2[1]

   opt.H1 <- optim(
    par=c(1+1e-3,rep(0,ncol(Z)+2)),
    fn=joint.loglike.one.part.model,
    gr = grad.joint.loglike.one.part.model,
    method = 'L-BFGS-B',
    lower = c(1+1e-8,rep(-Inf,ncol(Z)+1),-Inf),
    upper = c(Inf,rep(Inf,ncol(Z)+2)),
    control= list(maxit = 1e3),
    para.nuisance = c(par2[-c(1:3)],par2[1:2],q,par2[3]),
    X.aug=X.aug,Y=Y,prod.mat=prod.mat,
    X.test.coeff.index=X.test.coeff.index,
    Z=Z,dat.surv=dat.surv, 
    Z.test.coeff.index=Z.test.coeff.index,
    nsample = nsample,N= N,
    gh.weights = gh.weights,
    gh.nodes=gh.nodes,
    quad.n=quad.n,
    hessian = TRUE
  )
  est.H1 = c(par2,opt.H1$par)
  names(est.H1) = c(names(par2),'scaleK',
                          paste0('gamma',0:ncol(Z)),'delta2'
  )
	est.hessian = adiag(hes2,opt.H1$hessian)  
	status = unique(c(opt.H1.Sbeta$convergence,opt.H1$convergence))
	I.inverse = solve(est.hessian +1e-8)
	vars = diag(I.inverse)
	vars[which(vars<0)] = -vars[which(vars<0)]
	SEs = sqrt(vars)    
  Wald.abundance = ifelse(is.na(SEs[length(SEs)]),NA,(est.H1[length(SEs)]/SEs[length(SEs)])^2)
  pval.Wald.abundance = pchisq(Wald.abundance,df =1,lower.tail = FALSE)
  
  d$par.est.Wald = est.H1
  d$SEs = SEs
  d$est.hessian = est.hessian
  d$I.inverse = I.inverse
  d$Wald.Ts = c(Wald.abundance = Wald.abundance)
  d$pvals.Wald = c(pval.Wald.abundance=pval.Wald.abundance)
  d$status.Wald = status
  return(d)
}
