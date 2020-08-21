loglik.B.change.para.order <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                    gh.weights,gh.nodes,
                    quad.n){
  Y <- as.vector(Y)
  ncov = ncol(X.aug)-1
  beta1 = para[1:ncov]
  
  se2 <- para[ncov+1] 
  phi  <- para[ncov+2]
  q <- para[ncov+3]
  beta0 <- para[ncov+4]
  beta <- as.matrix(c(beta0,beta1))
  u <- 1 / (1 + exp(-(X.aug %*% beta[,rep(1,quad.n)] + gh.nodes * se2 * sqrt(2))))
  suppressWarnings(
    log.i <- -(phi-1)*log(q)+ lgamma(phi)-lgamma(u*phi)-lgamma((1-u)*phi) + (u*phi-1)*log(Y)+((1-u)*phi-1)*log(q-Y)
  )
  log.i[is.infinite(log.i)] <- 0
  logL <- sum(log(rowSums(gh.weights / sqrt(pi) * exp(prod.mat %*% log.i),
                          na.rm=TRUE)),
              na.rm=TRUE)
  return(-logL)
}
