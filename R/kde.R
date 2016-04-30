# bandwidth selectors
bandwidth.nrd2d <- function (x, is1D = TRUE) {
  if (!is1D) {
    1.06 * min(sqrt(var(x)), (quantile(x, c(0.25, 0.75))[2L] - quantile(x, c(0.25, 0.75))[1L]) / 1.34) * length(x) ^ (-1/6)
  } else {
    1.06 * min(sqrt(var(x)), (quantile(x, c(0.25, 0.75))[2L] - quantile(x, c(0.25, 0.75))[1L]) / 1.34) * length(x) ^ (-1/5)
  }
}

# 1-d kernel density estimate at a point x.new
kde1d.new <- function (x, h, x.new) {
  nx <- length(x)
  if (any(!is.finite(x))) {
    stop("Error: missing or infinite values in the data are not allowed")
  }
  gx <- x.new
  if (missing(h)) {
    h <- bandwidth.nrd2d(x)
  }
  ax <- outer(gx, x, "-") / h
  z <- apply(dnorm(ax), 1, mean) / h
  return(list(x = gx, y = z))
}

# 2-d kernel density estimate at a point (x.new, y.new)
kde2d.new <- function (x, y, h, x.new, y.new, lims = c(range(x), range(y))) {
  nx <- length(x)
  n.new <- length(x.new)
  if (length(y) != nx) {
    stop("Error: data vectors must be the same length")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Error: missing or infinite values in the data are not allowed")
  }
  if (any(!is.finite(lims))) {
    stop("Error: only finite values are allowed in 'lims'")
  }
  gx <- x.new
  gy <- y.new
  if (missing(h)){
    h <- c(bandwidth.nrd2d(x, FALSE), bandwidth.nrd2d(y, FALSE))
  }
  ax <- outer(gx, x, "-") / h[1L]
  ay <- outer(gy, y, "-") / h[2L]
  z <- matrix(dnorm(ax), n.new, nx) %*% t(matrix(dnorm(ay), n.new, nx)) / (nx * h[1L] * h[2L])
  return(list(x = gx, y = gy, z = z))
}
