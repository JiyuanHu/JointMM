fit.JointMM.fullmodel <-
function(d,quad.n=quad.n){
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
  Z.test.coeff.index =d$Z.test.coeff.index
  
  verbose=FALSE
  opt.H1.logit <- nlminb(start=c(.5,rep(0,sum(!X.test.coeff.index))), ## s1,alpha
                           objective=loglik.A,
                           lower = c(0.00001,
                                     rep(-Inf,sum(!X.test.coeff.index))),
                           upper = Inf,
                           X.aug=X.aug,Y=Y,prod.mat=prod.mat,
                           X.test.coeff.index = X.test.coeff.index,
                           gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                           control=list(trace=ifelse(verbose,2,0))
    )
    par1 = opt.H1.logit$par
    se1.hat <- par1[1]
    opt.H1.Sbeta <- nlminb(start=c(.5,20,1-1e-3,rep(0,sum(!X.test.coeff.index))), 
							objective=loglik.B,
                           lower = c(rep(1e-3,2),max(Y)+1e-2,
                                     rep(-Inf,sum(!X.test.coeff.index))),
                           upper = c(rep(Inf,2),1,rep(Inf,sum(!X.test.coeff.index))),
                           X.aug=X.aug,Y=Y,prod.mat=prod.mat,
                           X.test.coeff.index = X.test.coeff.index,
                           gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                           control=list(trace=ifelse(verbose,2,0))
    )
    opt.H1.Sbeta
    par2 = opt.H1.Sbeta$par
    se2.hat <- par2[1]
    start = c(1.1,.5,rep(0,ncol(Z)),rep(0,sum(!Z.test.coeff.index)))
    opt.H1.Surv <- nlminb(start = start,
      objective=loglik.C,
      lower = c(1+1e-5,1e-3,rep(-Inf,ncol(Z)+sum(!Z.test.coeff.index))),
      upper = c(Inf, 1,rep(Inf,ncol(Z)+sum(!Z.test.coeff.index))),
      se1.hat= se1.hat,se2.hat = se2.hat,
      Z=Z,dat.surv=dat.surv, 
      Z.test.coeff.index=rep(TRUE,2),
      nsample = nsample,
      bivar.gherm=bivar.gherm,
      control=list(trace=ifelse(verbose,2,0))
    )
    par3 = opt.H1.Surv$par[1:3]
    start2 = c(par1,par2[-c(1:4)],par2[1:4],par3,1e-2,1e-2) #try2     
    start2.H0 = c(par1,par2[-c(1:4)],par2[1:4],par3,0,0)
    names(start2) = c(c('se1',paste0('alpha',1:ncol(X.aug))),
                      c(paste0('beta',2:ncol(X.aug)),'se2','phi','q','beta1'),
                      c('shapeK','scaleB','gamma','delta1','delta2'))
    d$par.est.marginal = start2
  lower = c(c(1e-5,rep(-Inf,sum(!X.test.coeff.index))),
            c(rep(-Inf,sum(!X.test.coeff.index)-1),rep(1e-3,2),max(Y)+1e-3,-Inf),
            c(1+1e-5,1e-3,rep(-Inf,ncol(Z)),#rep(-Inf,sum(!Z.test.coeff.index)))
              rep(-Inf,sum(!Z.test.coeff.index)))
  )
  upper = c(rep(Inf,sum(!X.test.coeff.index)+1),
            c(rep(Inf,sum(!X.test.coeff.index)-1),rep(Inf,2),1,Inf),
            c(Inf, 1,rep(Inf,ncol(Z)+sum(!Z.test.coeff.index))))
  opt.H1 <- optim(
    par=start2, 
    fn=joint.loglike,
    gr = grad.joint.loglike,
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
  Tmax=1
  if((sum(is.na(est.hessian))| status !=0) & Tmax <5){ #use another start value.
    print('NaN found in the hessian matrix')
    if(Tmax ==1){
      start2 = c(par1,par2[-c(1:4)],par2[1],20,1-1e-3,par2[4],par3,1e-2,1e-2) #try1
    }else{
      start2 = c(par1,par2[-c(1:4)],par2[1],20,1-1e-3,par2[4],par3,1e-2,1e-2) - runif(15,0,1e-2)
    }
    opt.H1 <- optim(#start=start1,
      par=start2, 
      fn=joint.loglike,
      gr = grad.joint.loglike,
      method = 'L-BFGS-B',
      lower = lower,
      upper = upper,
	  control= list(maxit = 1e3),
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
    Tmax = Tmax+1
  }
  print(c('Number of optimization to converge:',Tmax))
  I.inverse = solve(opt.H1$hessian)
  SEs = sqrt(diag(I.inverse)) #standard errors
  ind = (ncol(I.inverse)-1):ncol(I.inverse)
  Wald.all = ifelse(sum(is.na(I.inverse[ind,ind])),NA,t(est.H1[ind]) %*% solve(I.inverse[ind,ind]) %*% est.H1[ind])
  
  Wald.presence.absence = ifelse(is.na(SEs[ind[1]]),NA,(est.H1[ind[1]]/SEs[ind[1]])^2)
  Wald.abundance = ifelse(is.na(SEs[ind[2]]),NA,(est.H1[ind[2]]/SEs[ind[2]])^2)

  pval.Wald.all = pchisq(Wald.all,df =2,lower.tail = FALSE)
  pval.Wald.presence.absence = pchisq(Wald.presence.absence,df =1,lower.tail = FALSE)
  pval.Wald.abundance = pchisq(Wald.abundance,df =1,lower.tail = FALSE)
  status.H1 = opt.H1$convergence
  
  d$par.est.Wald = est.H1
  d$SEs = SEs
  d$est.hessian = est.hessian
  d$I.inverse = I.inverse
  names(d$par.est.Wald) = names(d$SEs) = c('se1',paste0('alpha',1:ncol(X.aug)),
                          paste0('beta',2:ncol(X.aug)),'se2','phi','q','beta1',
                          'shapeK','scaleB','gamma','delta1','delta2')
  d$Wald.Ts = c(Wald.all = Wald.all,Wald.delta1 = Wald.presence.absence,Wald.delta2 = Wald.abundance)
  d$pvals.Wald = c(pval.Wald.all = pval.Wald.all,pval.Wald.delta1 = pval.Wald.presence.absence,
                   pval.Wald.delta2 = pval.Wald.abundance)
  
  d$status.Wald = status.H1
  return(d)
}
