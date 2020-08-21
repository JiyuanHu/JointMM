loglik.B_i <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                      bivar.gherm,N){
  Y <- as.vector(Y)
  se2 <- para[1] 
  phi  <- para[2] 
  q <- para[3] 
  beta <- rep(NA,length(X.test.coeff.index))
  beta[X.test.coeff.index]  <- 0 
  beta[!X.test.coeff.index] <- para[-(1:3)][1:sum(!X.test.coeff.index)]
  beta <- as.matrix(beta)

  bivar.nodesb = matrix(rep(bivar.gherm$nodes[,2],N),nrow = N,byrow = TRUE)
  quad.n = ncol(bivar.nodesb)
  u <- 1 / (1 + exp(-(X.aug %*% beta[,rep(1,quad.n)] + bivar.nodesb * se2* sqrt(2)))) 
  suppressWarnings(
    log.i <- -(phi-1)*log(q)+ lgamma(phi)-lgamma(u*phi)-lgamma((1-u)*phi) + (u*phi-1)*log(Y)+((1-u)*phi-1)*log(q-Y)
  )
  log.i[is.infinite(log.i)] <- 0
  logL_Bi = prod.mat %*% log.i
  return(logL_Bi)
}
