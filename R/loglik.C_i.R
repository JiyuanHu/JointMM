loglik.C_i <-
function(para,Z,dat.surv, Z.test.coeff.index,bivar.gherm,nsample){
  time <- as.vector(dat.surv$time)
  disease<- as.vector(dat.surv$disease)
  se1 <- para[1]
  se2 <- para[2]
  k <- para[3] 
  b <- para[4] 
  gamma <- para[5:(4+ncol(Z))]
  gamma = as.matrix(gamma)
  
  delta <- rep(NA,length(Z.test.coeff.index))
  delta[Z.test.coeff.index]  <- 0 
  delta[!Z.test.coeff.index] <- para[-(1:(4+ncol(Z)))][1:sum(!Z.test.coeff.index)]
  delta1 = delta[1]
  delta2 = delta[2]
  bivar.nodesa = matrix(rep(bivar.gherm$nodes[,1],nsample),nrow = nsample,byrow = TRUE)
  bivar.nodesb = matrix(rep(bivar.gherm$nodes[,2],nsample),nrow = nsample,byrow = TRUE)
  quad.n = ncol(bivar.nodesa)
  time.long = matrix(rep(time,quad.n),ncol=quad.n, byrow = FALSE)
  disease.long = matrix(rep(disease,quad.n),ncol=quad.n, byrow = FALSE)
  log.i <- disease.long* (log(b)+log(k)+(k-1)*log(time.long)+ Z %*% gamma[,rep(1,quad.n)]+
                            delta1*bivar.nodesa*se1*sqrt(2)+delta2*bivar.nodesb*se2*sqrt(2)) - 
    b *time.long^k * exp(Z %*% gamma[,rep(1,quad.n)]+
                             delta1*bivar.nodesa*se1*sqrt(2)+
                             delta2*bivar.nodesb*se2*sqrt(2))
  log.i[is.infinite(log.i)] <- 0
  logL_Ci = log.i
  return(logL_Ci)
}
