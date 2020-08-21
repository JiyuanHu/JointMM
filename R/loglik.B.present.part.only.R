loglik.B.present.part.only <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                    gh.weights,gh.nodes,
                    quad.n){
  Y <- as.vector(Y)
  se2 <- para[1] 
  phi  <- para[2]
  eps = function(maxY){
    ifelse(maxY>0.01,0.01,ifelse(maxY>1e-3,1e-3,1e-4))
  }
  q <- max(Y) +eps (max(Y))
  beta <- rep(NA,length(X.test.coeff.index))
  beta[X.test.coeff.index]  <- 0 ## initialized as 0
  beta[!X.test.coeff.index] <- para[-(1:2)][1:sum(!X.test.coeff.index)]
  beta <- as.matrix(beta)
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
