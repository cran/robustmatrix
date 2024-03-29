% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmatnorm.R
\name{rmatnorm}
\alias{rmatnorm}
\title{Simulate from a Matrix Normal Distribution}
\usage{
rmatnorm(n, mu = NULL, cov_row, cov_col)
}
\arguments{
\item{n}{the number of samples required.}

\item{mu}{a \eqn{p \times q} matrix containing the means.}

\item{cov_row}{a \eqn{p \times p} positive-definite symmetric matrix specifying the rowwise covariance matrix}

\item{cov_col}{a \eqn{q \times q} positive-definite symmetric matrix specifying the columnwise covariance matrix}
}
\value{
If \eqn{n = 1} a matrix with \eqn{p} rows and \eqn{q} columns, o
otherwise a 3d array of dimensions \eqn{(p,q,n)} with a sample in each slice.
}
\description{
Simulate from a Matrix Normal Distribution
}
\examples{
n = 1000; p = 2; q = 3
mu = matrix(rep(0, p*q), nrow = p, ncol = q)
cov_row = matrix(c(5,2,2,4), nrow = p, ncol = p)
cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
X[,,9] #printing the 9th sample.
}
