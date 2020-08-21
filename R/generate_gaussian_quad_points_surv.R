generate_gaussian_quad_points_surv <-
function(quad.n){
  gherm <- gauss.quad(quad.n, kind = "hermite")
  nodes = gherm$nodes
  weights = gherm$weights
  weights = as.numeric(outer(weights,weights))
  nodes = cbind(rep(nodes,quad.n),rep(nodes,each = quad.n))
  bivar.gherm = list(nodes = nodes, weights = weights)
  return(bivar.gherm)
}
