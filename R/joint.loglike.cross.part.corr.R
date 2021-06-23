joint.loglike.cross.part.corr <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                         Z, dat.surv, Z.test.coeff.index,
                         nsample,N,
                         quad.n,bivar.gherm){
  se1 = para[1]
  alpha = para[2:(sum(!X.test.coeff.index)+1)]
  paraA = c(se1,alpha)
  
  se2 = para[sum(!X.test.coeff.index)+2]
  phi = para[sum(!X.test.coeff.index)+3]
  q = para[sum(!X.test.coeff.index)+4]
  beta =  para[(sum(!X.test.coeff.index)+5):(2*sum(!X.test.coeff.index)+4)]
  #beta = c(para[2*sum(!X.test.coeff.index)+4],beta) 
  delta0 = para[2*sum(!X.test.coeff.index)+5]
  
  paraB = c(delta0,se1,se2,phi,q,beta)
  
  ##paraBi = c(delta0,se1,se2,phi,q,beta)
  bivar.gh.weights = matrix(rep(bivar.gherm$weights,nsample),nrow = nsample,byrow = TRUE)

  logLike.partA <- loglik.A_i(paraA,X.aug,Y,prod.mat,X.test.coeff.index,
                              bivar.gherm,N)
  logLike.partB <- loglik.B_i.cross.part.corr(para=paraB,X.aug,Y,prod.mat,X.test.coeff.index,
                              bivar.gherm,N)
  joint.logL = sum(log(rowSums(bivar.gh.weights/pi*exp(logLike.partA+logLike.partB))))
  return(-joint.logL)
}
