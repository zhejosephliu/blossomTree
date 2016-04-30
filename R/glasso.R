# path regularization of graphical lasso
glasso.path <- function (x.trn, pool.gl) {
  n <- nrow(x.trn)
  d <- ncol(x.trn)
  if (is.null(pool.gl)) {
    pool.gl <- seq(0.005, 0.2, length = 50)
  }
  glPath <- list()
  out <- glasso(var(x.trn), rho = pool.gl[1], penalize.diagonal = FALSE, start = "cold")
  glPath[[1]] <- out$wi
  for (j in c(2:length(pool.gl))) {
    out <- glasso(var(x.trn), rho = pool.gl[j], penalize.diagonal = FALSE, start = "warm", w.init = out$w, wi.init = out$wi)
    glPath[[j]] <- out$wi
  }
  return(list(glPath = glPath))
}

# fit graphical lasso
loglike.glasso <- function (x.st1, x.st2, pool.gl, refit = FALSE) {
  if (is.null(pool.gl)) {
    pool.gl <- seq(0.005, 0.2, length = 50)
  }
  d <- ncol(x.st1)
  gl.path <- glasso.path(x.trn = x.st1, pool.gl)$glPath
  loglike.gl <- rep(0, length(gl.path))
  loglike.normal <- rep(0, length(gl.path))
  sparsity.gl <- rep(0, length(gl.path))
  mu.hat <- apply(x.st1, 2, mean)
  for (j in c(1:length(gl.path))) {
    sigma.hat <- solve(gl.path[[j]])
    sigma.hat <- (sigma.hat + t(sigma.hat)) / 2
    tmp <- log(dmvnorm(x = x.st2, mean = mu.hat, sigma = sigma.hat))
    loglike.gl[j] <- mean(tmp)
    sparsity.gl[j] <- (sum(abs(gl.path[[j]]) > 1e-6) - d) / 2
    if (sum(gl.path[[j]] == 0) > 0) {
      out <- glasso(var(x.st1), rho = 1e-10, zero = which(gl.path[[j]] == 0, arr.ind = TRUE), penalize.diagonal = FALSE, start = "cold")
    } else {
      out <- glasso(var(x.st1), rho = 1e-10, penalize.diagonal = FALSE, start = "cold")
    }
    sigma.normal <- solve(out$wi)
    sigma.normal <- (sigma.normal + t(sigma.normal)) / 2
    loglike.normal[j] <- mean(log(dmvnorm(x = x.st2, mean = mu.hat, sigma = sigma.normal)))
  }
  if (refit) {
    j <- which(loglike.normal == max(loglike.normal))[1]
  } else {
    j <- which(loglike.gl == max(loglike.gl))[1]
  }
  g.gl <- graph.adjacency((gl.path[[j]] + t(gl.path[[j]])), weighted = TRUE, mode = "undirected", diag = FALSE)
  return(list(loglike.gl = loglike.gl, loglike.normal = loglike.normal, sparsity.gl = sparsity.gl, gl.path = gl.path, g.gl = g.gl))
}
