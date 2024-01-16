#' Simulate from a Matrix Normal Distribution
#'
#' @param n the number of samples required.
#' @param mu a \eqn{p \times q} matrix containing the means.
#' @param cov_row a \eqn{p \times p} positive-definite symmetric matrix specifying the rowwise covariance matrix
#' @param cov_col a \eqn{q \times q} positive-definite symmetric matrix specifying the columnwise covariance matrix
#' @return If \eqn{n = 1} a matrix with \eqn{p} rows and \eqn{q} columns, o
#' otherwise a 3d array of dimensions \eqn{(p,q,n)} with a sample in each slice.
#' @export
#'
#' @examples
#' n = 1000; p = 2; q = 3
#' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
#' cov_row = matrix(c(5,2,2,4), nrow = p, ncol = p)
#' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
#' X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
#' X[,,9] #printing the 9th sample.
rmatnorm <- function(n, mu = NULL, cov_row, cov_col){
  p = nrow(cov_row)
  q = nrow(cov_col)
  if(is.null(mu)|all(mu == 0)){
    mu <- matrix(0, nrow = p, ncol = q)
  }

  cov_row_eigen <- eigen(cov_row)
  cov_col_eigen <- eigen(cov_col)

  cov_row_sqrt <- cov_row_eigen$vectors %*% diag(sqrt(cov_row_eigen$values)) %*% t(cov_row_eigen$vectors)
  cov_col_sqrt <- cov_col_eigen$vectors %*% diag(sqrt(cov_col_eigen$values)) %*% t(cov_col_eigen$vectors)

  X_sample <- array(dim = c(p,q,n))
  for(i in 1:n){
    X_sample[,,i] <- mu + cov_row_sqrt%*%matrix(rnorm(n = p*q), nrow = p, ncol = q)%*%t(cov_col_sqrt)
  }
  if(n == 1){
    X_sample <- matrix(X_sample, nrow = p, ncol = q)
  }
  X_sample
}
