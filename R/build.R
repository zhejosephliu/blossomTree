# construct blossom tree
pure.blossom.tree <- function (xtrain, xheld, g, wt, bound = 8, prune = FALSE, lambda = NULL, refit = FALSE, verbose = TRUE) {
  max_edges <- Inf
  d <- ncol(xtrain)
  seq_adj <- list()
  forestnode <- rep(0, d)
  node <- rep(1, d)
  blossom <- rep(0, d)
  parcor <- matrix(Inf, d, d)
  new_wts <- c()
  loglike_btde <- c()
  sub_g <- graph.empty(vcount(g))
  sub_g <- as.undirected(sub_g)
  wt_sort_ids <- sort(wt, decreasing = T, index.return = T)$ix
  res <- loglike.glasso(xtrain, xheld, pool.gl = lambda, refit)
  loglike_btde <- c(loglike_btde, max(res$loglike.normal))
  adj <- matrix(0, d, d)
  bedges <- ends(res$g.gl, 1:ecount(res$g.gl))
  for (k in 1:dim(bedges)[1]) {
    adj[bedges[k, 1], bedges[k, 2]] <- 1
    adj[bedges[k, 2], bedges[k, 1]] <- 1
  }
  seq_adj[[length(seq_adj) + 1]] <- adj
  cluster_id <- 1:vcount(g)
  cluster_size <- array(1, dim = vcount(g))
  if (max_edges >= ecount(g)) {
    max_edges <- ecount(g)
  }
  if (ecount(g) > 0) {
    for (i in 1:max_edges) {
      cur_vs <- ends(g, wt_sort_ids[i])
      v1 <- cur_vs[1]
      v2 <- cur_vs[2]
      v1i <- v1
      v2i <- v2
      cur_wt <- wt[wt_sort_ids[i]]
      if (cluster_id[v1i] != cluster_id[v2i] && (!prune || cluster_size[cluster_id[v1i]] + cluster_size[cluster_id[v2i]] <= bound)) {
        cluster_size[cluster_id[v1i]] <- cluster_size[cluster_id[v1i]] + cluster_size[cluster_id[v2i]]
        old_id <- cluster_id[v2i]
        new_id <- cluster_id[v1i]
        for (j in 1:length(cluster_id)) {
          if (cluster_id[j] == old_id) {
            cluster_id[j] <- new_id
          }
        }
        new_wts <- c(new_wts, cur_wt)
        sub_g <- add.edges(sub_g, c(v1, v2))
        edges <- ends(sub_g, 1:ecount(sub_g))
        edgeset1 <- c(edges[which(edges[, 1] == v1), 2], edges[which(edges[,2] == v1), 1])
        edgeset2 <- c(edges[which(edges[, 1] == v2), 2], edges[which(edges[,2] == v2), 1])
        forestnode[v1] <- 1
        forestnode[v2] <- 1
        node[v1] <- 0
        node[v2] <- 0
        blossom[v1] <- 0
        blossom[v2] <- 0
        for (t in 1:d) {
          if (node[t] == 0) {
            next
          }
          for (v in edgeset1) {
            parcor[t, v1] <- min(parcor[t, v1], abs(np.partial.cor(xtrain[, t], xtrain[, v1], xtrain[, v])))
          }
          for (u in edgeset2) {
            parcor[t, v2] <- min(parcor[t, v2], abs(np.partial.cor(xtrain[, t], xtrain[, v2], xtrain[, u])))
          }
          blossom[t] <- (which(forestnode == 1))[which.max(parcor[t, which(forestnode == 1)])]
        }
        adj <- matrix(0, d, d)
        for (k in 1:dim(edges)[1]) {
          adj[edges[k, 1], edges[k, 2]] <- 1
          adj[edges[k, 2], edges[k, 1]] <- 1
        }
        log_blossoms <- 0
        for (j in 1:d) {
          if (forestnode[j] == 0) {
            next
          }
          index <- c(j, which(blossom == j))
          if (length(index) == 1) {
            log_blossoms <- log_blossoms + mean(log(dnorm(x = xheld[, index], mean = mean(xtrain[, index]), sd = sd(xtrain[, index]))))
            next
          }
          res <- loglike.glasso(xtrain[, index], xheld[, index], pool.gl = lambda, refit)
          if (refit) {
            log_blossoms <- log_blossoms + max(res$loglike.gl)
          } else {
            log_blossoms <- log_blossoms + max(res$loglike.normal)
          }
          if (ecount(res$g.gl) == 0) {
            next
          }
          bedges <- ends(res$g.gl, 1:ecount(res$g.gl))
          for (k in 1:dim(bedges)[1]) {
            adj[index[bedges[k, 1]], index[bedges[k, 2]]] <- 1
            adj[index[bedges[k, 2]], index[bedges[k, 1]]] <- 1
          }
        }
        seq_adj[[length(seq_adj) + 1]] <- adj
        loglike_btde_add <- loglike.blossom.tree(xtrain, xheld, sub_g, log_blossoms)
        if (verbose) {
          cat(paste('    add edge (', v1, ', ', v2, '),', '\t', 'loglike: ', round(loglike_btde_add, 3), '\n', sep = ""))
        }
        loglike_btde <- c(loglike_btde, loglike_btde_add)
      }
    }
  }
  if (ecount(sub_g) > 0) {
    E(sub_g)$weight <- new_wts
  }
  return(list(loglike_btde = loglike_btde, seq_adj = seq_adj))
}
