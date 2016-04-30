#' @import bisoreg glasso huge igraph MASS mvtnorm
#' @export

# fit blossom tree graphical model
blossomTree <- function (data, lambda = NULL, refit = FALSE, verbose = TRUE) {
  set.seed(100)
  npn.data <- huge.npn(data, verbose = verbose)
  index <- sample(1:dim(npn.data)[1])[1:round(dim(npn.data)[1] / 2)]
  xtrain <- npn.data[index, ]
  xheld <- npn.data[-index, ]
  num_kde_pt <- 128
  d <- ncol(xtrain)
  if (verbose) {
    cat("Estimating the pairwise negentropy....")
  }
  ne <- negentropy(xtrain, m = num_kde_pt)$ne
  if (verbose) {
    cat("done.\n")
  }
  wt <- t(ne)[lower.tri(t(ne))]
  g <- graph.full(d)
  E(g)$weight <- wt
  if (verbose) {
    cat("Constructing a family of blossom tree structures and computing held-out log-likelihood....\n")
  }
  fit <- pure.blossom.tree(xtrain, xheld, g, wt, bound = 8, prune = FALSE, lambda, refit, verbose)
  if (verbose) {
    cat("done.\n")
  }
  best.loglike <- max(fit$loglike_btde)
  best.adj <- graph.adjacency(fit$seq_adj[[which.max(fit$loglike_btde)]], mode = "undirected")
  return(list(loglike = fit$loglike_btde, adj = fit$seq_adj, best.loglike = best.loglike, best.adj = best.adj))
}
