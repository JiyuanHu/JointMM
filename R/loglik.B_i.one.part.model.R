loglik.B_i.one.part.model <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                                      gh.nodes,N){
  Y <- as.vector(Y)
  se2 <- para[1] 
  phi  <- para[2]
  q <- para[3] 
  beta <- rep(NA,length(X.test.coeff.index))
  beta[X.test.coeff.index]  <- 0 
  beta[!X.test.coeff.index] <- para[-(1:3)][1:sum(!X.test.coeff.index)]
  beta <- as.matrix(beta)
  
  quad.n = ncol(gh.nodes)
  u <- 1 / (1 + exp(-(X.aug %*% beta[,rep(1,quad.n)] + gh.nodes * se2* sqrt(2)))) 
  suppressWarnings(
    log.i <- -(phi-1)*log(q)+ lgamma(phi)-lgamma(u*phi)-lgamma((1-u)*phi) + (u*phi-1)*log(Y)+((1-u)*phi-1)*log(q-Y)
  )
  log.i[is.infinite(log.i)] <- 0
  logL_Bi = prod.mat %*% log.i
  return(logL_Bi)
}
