loglik.A <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                    gh.weights,gh.nodes,
                    quad.n){
  se1 <- para[1] 
  alpha <- rep(NA,length(X.test.coeff.index))
  alpha[X.test.coeff.index]  <- 0 
  alpha[!X.test.coeff.index] <- para[-1][1:sum(!X.test.coeff.index)]
  alpha <- as.matrix(alpha)
  p  <- 1 / (1 + exp(-(X.aug%*%alpha[,rep(1,quad.n)] + gh.nodes*se1*sqrt(2))) )
  p[p<10^(-7)] <- 10^(-7)
  p[p>(1-10^(-7))] <- (1-10^(-7))
  p[Y == 0,] <- 1 - p[Y == 0,]
  tmp = gh.weights / sqrt(pi) * exp(prod.mat %*% log(p))
  logL <- sum(log(rowSums(gh.weights / sqrt(pi) * exp(prod.mat %*% log(p)))))
  return(-logL)
}
