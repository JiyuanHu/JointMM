joint.loglike.one.part.model <-
function(para,para.nuisance,X.aug,Y,prod.mat,X.test.coeff.index,
                                                          Z, dat.surv, Z.test.coeff.index,
                                                          nsample,N,
                                                          gh.weights,gh.nodes,quad.n){
  beta =  para.nuisance[1:(sum(!X.test.coeff.index)-1)]
  se2 = para.nuisance[sum(!X.test.coeff.index)]
  phi = para.nuisance[sum(!X.test.coeff.index)+1]
  q = para.nuisance[sum(!X.test.coeff.index)+2]
  beta = c(para.nuisance[sum(!X.test.coeff.index)+3],beta) 
  paraB = c(se2,phi,q,beta)
  paraC = c(se2,para) 
  logLike.partB <- loglik.B_i.one.part.model(para=paraB,X.aug,Y,prod.mat,X.test.coeff.index,
                                             gh.nodes,N) 
  logLike.partC <- loglik.C_i.one.part.model(para = paraC,Z, dat.surv, Z.test.coeff.index,gh.nodes,nsample) 
  joint.logL = sum(log(rowSums(gh.weights/sqrt(pi)*exp(logLike.partB+logLike.partC))))
  return(-joint.logL)
}
