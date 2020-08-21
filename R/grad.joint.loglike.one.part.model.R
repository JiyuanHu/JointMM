grad.joint.loglike.one.part.model <-
function(para,para.nuisance,X.aug,Y,prod.mat,X.test.coeff.index,
                                                Z, dat.surv, Z.test.coeff.index,
                                                nsample,N,
                                                gh.weights,gh.nodes,quad.n){
  eps = 1e-4
  eps.mat = diag(eps,length(para))
  
  grad = NULL;
  for(j in 1:length(para)){
    grad[j] = 1/eps*(joint.loglike.one.part.model(para+eps.mat[,j],para.nuisance,X.aug,Y,prod.mat,X.test.coeff.index,Z, dat.surv, Z.test.coeff.index,nsample,N,gh.weights,gh.nodes,quad.n)-joint.loglike.one.part.model(para,para.nuisance,X.aug,Y,prod.mat,X.test.coeff.index,Z, dat.surv, Z.test.coeff.index, nsample,N, gh.weights,gh.nodes,quad.n)
    )
  }
  return(grad)
}
