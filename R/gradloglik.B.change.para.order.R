gradloglik.B.change.para.order <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                        gh.weights,gh.nodes,
                        quad.n){
  eps = 1e-4
  eps.mat = diag(eps,length(para))
  
  grad = NULL;
  for(j in 1:length(para)){
    grad[j] = 1/eps*(loglik.B.change.para.order(para+eps.mat[,j],X.aug,Y,prod.mat,X.test.coeff.index,
                              gh.weights,gh.nodes,quad.n)-
                       loglik.B.change.para.order(para,X.aug,Y,prod.mat,X.test.coeff.index,
                                gh.weights,gh.nodes,
                                quad.n)
    )
  }
  return(grad)
  
}
