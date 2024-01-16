#' Outlier explanation based on Shapley values for matrix-variate data
#'
#' \code{matrixShapley} decomposes the squared matrix Mahalanobis distance (\code{\link{mmd}}) into additive outlyingness contributions of
#' the rows, columns, or cell of a matrix \insertCite{mayrhofer2023multivariate,mayrhofer2024}{robustmatrix}.
#'
#' @param type Character. Either "row", "col", or "cell" (default) to compute rowwise, columnwise, or cellwise Shapley values.
#' @inheritParams mmd
#'
#' @return Rowwise, columnwise, or cellwise Shapley value(s).
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{mmd}}.
#'
#' @export
#'
#' @examples
#' n = 1000; p = 2; q = 3
#' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
#' cov_row = matrix(c(5,2,2,4), nrow = p, ncol = p)
#' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
#' X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
#' distances <- mmd(X, mu, cov_row, cov_col)
matrixShapley <- function(X, mu = NULL, cov_row, cov_col, inverted = FALSE, type = "cell"){
  if(is.null(mu)){
    mu <- array(0, dim = dim(X))
  }
  if(inverted){
    cov_row_inv <- cov_row
    cov_col_inv <- cov_col
  } else{
    cov_row_inv <- chol2inv(chol(cov_row))
    cov_col_inv <- chol2inv(chol(cov_col))
  }

  if(length(dim(X)) == 2){
    if(type == "row"){
      diag(cov_row_inv%*%(X-mu)%*%cov_col_inv%*%t(X-mu))
    } else if(type == "col"){
      diag(t(X-mu)%*%cov_row_inv%*%(X-mu)%*%cov_col_inv)
    } else{
      (X-mu)*cov_row_inv%*%(X-mu)%*%cov_col_inv
    }
  } else if(length(dim(X)) == 3){
    if(type == "row"){
      apply(X, 3, function(x) diag(cov_row_inv%*%(x-mu)%*%cov_col_inv%*%t(x-mu)))
    } else if(type == "col"){
      apply(X, 3, function(x) diag(t(x-mu)%*%cov_row_inv%*%(x-mu)%*%cov_col_inv))
    } else{
      shv <- array(NA, dim = dim(X))
      shv[] <- apply(X, 3, function(x) (x-mu)*cov_row_inv%*%(x-mu)%*%cov_col_inv)
      shv
    }
  } else{
    stop("X must be either a p times q matrix or a 3d array of dimensions (p,q,n)")
  }
}
