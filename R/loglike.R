# compute log-likelihood for blossom tree
loglike.blossom.tree <- function (xtrain, xheld, graph, log_blossoms) {
  d <- ncol(xtrain)
  n <- nrow(xtrain)
  m <- nrow(xheld)
  h1 <- rep(0, d)
  h2 <- rep(0, d)
  for (j in c(1:d)) {
    h1[j] <- bandwidth.nrd2d(xtrain[, j], is1D = T)
    h2[j] <- bandwidth.nrd2d(xtrain[, j], is1D = F)
  }
  K1 <- matrix(0, nrow = m, ncol = d)
  for (j in 1:d) {
    K1[, j] <- kde1d.new(x = xtrain[, j], h = h1[j], x.new = xheld[, j])$y
  }
  edges <- ends(graph, 1:ecount(graph))
  log_pair_score <- array(0, m)
  for (j in 1:ecount(graph)) {
    cur_edge_from <- edges[j, 1]
    cur_edge_to <- edges[j, 2]
    K2 <- kde2d.new(x = xtrain[, cur_edge_from], y = xtrain[, cur_edge_to], h = c(h2[cur_edge_from], h2[cur_edge_to]), x.new = xheld[, cur_edge_from], y.new = xheld[, cur_edge_to])
    cur_pair_score <- log(diag(K2$z)) - log(K1[, cur_edge_from]) - log(K1[, cur_edge_to])
    log_pair_score <- log_pair_score + cur_pair_score
  }
  return(mean(log_blossoms + log_pair_score))
}
