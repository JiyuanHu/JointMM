grad.joint.loglike <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                              Z, dat.surv, Z.test.coeff.index,
                              nsample,N,
                              quad.n,bivar.gherm){
  eps = 1e-4
  eps.mat = diag(eps,length(para))
  
  grad = NULL;
  for(j in 1:length(para)){
  grad[j] = 1/eps*(joint.loglike(para+eps.mat[,j],X.aug,Y,prod.mat,X.test.coeff.index,
                                 Z, dat.surv, Z.test.coeff.index,
                                 nsample,N,
                                 quad.n,bivar.gherm)-
                  joint.loglike(para,X.aug,Y,prod.mat,X.test.coeff.index,
                                   Z, dat.surv, Z.test.coeff.index,
                                   nsample,N,
                                   quad.n,bivar.gherm)
                  )
  }
  return(grad)
}
