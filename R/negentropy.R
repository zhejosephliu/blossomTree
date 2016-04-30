# negentropy
negentropy <- function (x, m) {
  d <- ncol(x)
  ne <- matrix(0, ncol = d, nrow = d)
  h2 <- rep(0, d)
  for (j in c(1:d)) {
    h2[j] <- bandwidth.nrd2d(x[, j], is1D = F)
    if (h2[j] == 0) {
      stop(paste('Error: not enough distinct values to compute bandwidth in feature', j))
      return(-1)
    }
  }
  for (j in c(1:(d - 1))) {
    for (k in c((j + 1):d)){
      cur_K2 <- kde2d(x[, j], x[, k], n = m, h = c(h2[j], h2[k]) * 4)
      cur_K2_val <- cur_K2$z + 1e-200
      x_range <- max(cur_K2$x) - min(cur_K2$x)
      y_range <- max(cur_K2$y) - min(cur_K2$y)
      first_term <- x_range * y_range * sum(cur_K2_val * log(cur_K2_val)) / (m - 1) ^ 2
      ne[j, k] <- first_term + log(var(x[, j]) * var(x[, k]) - cov(x[, j], x[, k]) ^ 2) / 2 + log(2 * pi) + 1
    }
  }
  return(list(ne = ne))
}
