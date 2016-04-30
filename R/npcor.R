# nonparametric partial correlation
np.partial.cor <- function (y, z, x) {
  return(cor(y - predict(loess.wrapper(x, y)), z - predict(loess.wrapper(x, z))))
}
