gradloglik.B.present.part.only <-
function(para,X.aug,Y,prod.mat,X.test.coeff.index,
                        gh.weights,gh.nodes,
                        quad.n){
  eps = 1e-4
  eps.mat = diag(eps,length(para))
  
  grad = NULL;
  for(j in 1:length(para)){
    grad[j] = 1/eps*(loglik.B.present.part.only(para+eps.mat[,j],X.aug,Y,prod.mat,X.test.coeff.index,
                              gh.weights,gh.nodes,quad.n)-
                       loglik.B.present.part.only(para,X.aug,Y,prod.mat,X.test.coeff.index,
                                gh.weights,gh.nodes,
                                quad.n)
    )
  }
  return(grad)
  
}
