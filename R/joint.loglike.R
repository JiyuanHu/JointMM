joint.loglike <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                         Z, dat.surv, Z.test.coeff.index,
                         nsample,N,
                         quad.n,bivar.gherm){
  se1 = para[1]
  alpha = para[2:(sum(!X.test.coeff.index)+1)]
  paraA = c(se1,alpha)
  beta =  para[(sum(!X.test.coeff.index)+2):(2*sum(!X.test.coeff.index))]
  se2 = para[2*sum(!X.test.coeff.index)+1]
  phi = para[2*sum(!X.test.coeff.index)+2]
  q = para[2*sum(!X.test.coeff.index)+3]
  beta = c(para[2*sum(!X.test.coeff.index)+4],beta) 
  paraB = c(se2,phi,q,beta)
  paraC = c(se1,se2,para[-c(1:(2*sum(!X.test.coeff.index)+4))])
  bivar.gh.weights = matrix(rep(bivar.gherm$weights,nsample),nrow = nsample,byrow = TRUE)

  logLike.partA <- loglik.A_i(paraA,X.aug,Y,prod.mat,X.test.coeff.index,
                              bivar.gherm,N)
  logLike.partB <- loglik.B_i(para=paraB,X.aug,Y,prod.mat,X.test.coeff.index,
                              bivar.gherm,N)
  logLike.partC <- loglik.C_i(para = paraC,Z, dat.surv, Z.test.coeff.index,bivar.gherm,nsample)
  joint.logL = sum(log(rowSums(bivar.gh.weights/pi*exp(logLike.partA+logLike.partB+logLike.partC))))
  return(-joint.logL)
}
