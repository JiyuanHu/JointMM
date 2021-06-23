loglik.B_i.cross.part.corr <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                      bivar.gherm,N){
  Y <- as.vector(Y)
  delta0 <- para[1]
  se1 <- para[2]
  se2 <- para[3] 
  phi  <- para[4] 
  q <- para[5] 
  beta <- rep(NA,length(X.test.coeff.index))
  beta[X.test.coeff.index]  <- 0 
  beta[!X.test.coeff.index] <- para[-(1:5)][1:sum(!X.test.coeff.index)]
  beta <- as.matrix(beta)
#paraBi = c(delta0,se1,se2,phi,q,beta)
  bivar.nodesa = matrix(rep(bivar.gherm$nodes[,1],N),nrow = N,byrow = TRUE)
  bivar.nodesb = matrix(rep(bivar.gherm$nodes[,2],N),nrow = N,byrow = TRUE)
  quad.n = ncol(bivar.nodesb)
  u <- 1 / (1 + exp(-(X.aug %*% beta[,rep(1,quad.n)] +delta0*bivar.nodesa * se1* sqrt(2)+ bivar.nodesb * se2* sqrt(2)))) 
  suppressWarnings(
    log.i <- -(phi-1)*log(q)+ lgamma(phi)-lgamma(u*phi)-lgamma((1-u)*phi) + (u*phi-1)*log(Y)+((1-u)*phi-1)*log(q-Y)
  )
  log.i[is.infinite(log.i)] <- 0
  logL_Bi = prod.mat %*% log.i
  return(logL_Bi)
}
