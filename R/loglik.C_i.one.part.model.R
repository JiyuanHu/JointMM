loglik.C_i.one.part.model <-
function(para,Z,dat.surv, Z.test.coeff.index,gh.nodes,nsample){
  Z.aug = cbind(intersept = 1,Z)
  time <- as.vector(dat.surv$time)
  disease<- as.vector(dat.surv$disease)
  se2 <- para[1]
  k <- para[2] 
  gamma <- para[3:(2+ncol(Z.aug))]
  gamma = as.matrix(gamma)
  
  delta <- para[-(1:(2+ncol(Z.aug)))]
  delta2 = delta
  gh.nodes.survival = gh.nodes[1:nsample,]
  quad.n = ncol(gh.nodes)
  time.long = matrix(rep(time,quad.n),ncol=quad.n, byrow = FALSE)
  disease.long = matrix(rep(disease,quad.n),ncol=quad.n, byrow = FALSE)
  log.i <- disease.long* (log(k)+(k-1)*log(time.long)+ Z.aug %*% gamma[,rep(1,quad.n)]+
                            delta2*gh.nodes.survival*se2*sqrt(2)) - 
    time.long^k * exp(Z.aug %*% gamma[,rep(1,quad.n)]+
                        delta2*gh.nodes.survival*se2*sqrt(2))
  log.i[is.infinite(log.i)] <- 0
  logL_Ci = log.i
  return(logL_Ci)
}
