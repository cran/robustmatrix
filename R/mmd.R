#' Matrix Mahalanobis distance
#'
#' @param inverted Logical. FALSE by default.
#' If TRUE \code{cov_row} and \code{cov_col} are supposed to contain the inverted rowwise and columnwise covariance matrices, respectively.
#' @inheritParams rmatnorm
#' @inheritParams mmle
#'
#' @return Squared Mahalanobis distance(s) of observation(s) in \code{X}.
#' @export
#'
#' @examples
#' n = 1000; p = 2; q = 3
#' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
#' cov_row = matrix(c(1,0.5,0.5,1), nrow = p, ncol = p)
#' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
#' X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
#' ind <- sample(1:n, 0.3*n)
#' X[,,ind] <- rmatnorm(n = length(ind), matrix(rep(10, p*q), nrow = p, ncol = q), cov_row, cov_col)
#' distances <- mmd(X, mu, cov_row, cov_col)
#' plot(distances)
#' abline(h = qchisq(0.99, p*q), lty = 2, col = "red")
mmd <- function(X, mu, cov_row, cov_col, inverted = FALSE){
  if(length(dim(X)) == 2){
    MMD(X = X, mu = mu, cov_row = cov_row, cov_col = cov_col, inverted = inverted)
  } else if(length(dim(X)) == 3){
    as.vector(TensorMMD(X = X, mu = mu, cov_row = cov_row, cov_col = cov_col, inverted = inverted))
  } else{
    stop("X must be either a p times q matrix or a 3d array of dimensions (p,q,n)")
  }
}
