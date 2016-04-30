\name{blossomTree}
\alias{blossomTree}
\title{Bossom tree graph estimation and density estimation}
\description{
The main function for bossom tree graph estimation and density estimation.
}
\usage{
blossomTree(data, lambda = NULL, refit = FALSE, verbose = TRUE)
}
\arguments{
\item{data}{
An \code{n} by \code{d} data matrix, where \code{n} is the sample size and \code{d} is the dimension.
}
\item{lambda}{
A sequence of positive numbers to control the regularization of the graphical lasso.
}
\item{refit}{
If \code{refit = TRUE}, the refit graphical lasso is applied. The default value is \code{FALSE}.
}
\item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
Partition the data matrix into two subsamples with equal size by random allocation. The first subsample is used to construct a family of blossom tree graphs and the corresponding density estimators, while the second subsample is then used to determine the optimal blossom tree by maximizing the held-out log-likelihood.

The refit graphical lasso is a two-step procedure: in the first step, a sparse inverse covariance matrix is obtained by the graphical lasso; in the second step, a Gaussian model is refit without regularization, but enforcing the sparsity pattern obtained in the first step.
}
\value{
\item{loglike}{
Held-out log-likelihood for a family of blossom tree density estimators with \code{m=1,...,d-1} forest edges.
}
\item{adj}{
Adjacency matrices of a family of blossom tree graphs with \code{m=1,...,d-1} forest edges.
}
\item{best.loglike}{
Maximum of the held-out log-likelihood.
}
\item{best.adj}{
Adjacency matrix of the optimal blossom tree graph corresponding to the maximum of the held-out log-likelihood.
}
}
\references{
Zhe Liu and John Lafferty. Blossom tree graphical models. \emph{Advances in Neural Information Processing Systems}, 27:1458--1465, 2014.
}
\examples{
library(igraph)

fit <- blossomTree(data)
bt <- fit$best.adj
plot(bt)
}
