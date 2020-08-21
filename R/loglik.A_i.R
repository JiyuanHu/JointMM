loglik.A_i <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,bivar.gherm,N){
  se1 <- para[1] 
  alpha <- rep(NA,length(X.test.coeff.index))
  alpha[X.test.coeff.index]  <- 0 
  alpha[!X.test.coeff.index] <- para[-1][1:sum(!X.test.coeff.index)]
  alpha <- as.matrix(alpha)
  bivar.nodesa = matrix(rep(bivar.gherm$nodes[,1],N),nrow = N,byrow = TRUE)
  quad.n = ncol(bivar.nodesa)
  p  <- 1 / (1 + exp(-(X.aug%*%alpha[,rep(1,quad.n)] + bivar.nodesa*se1*sqrt(2))) )
  p[p<10^(-7)] <- 10^(-7)
  p[p>(1-10^(-7))] <- (1-10^(-7))
  p[Y == 0,] <- 1 - p[Y == 0,]
  logL_Ai = prod.mat %*% log(p)
  return(logL_Ai)
}
