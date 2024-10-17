#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <thread>
#include <iostream>

using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

///////////////////////////////////////////////////////////////////////////////////////
// Creating classes for matrix-variate estimators, MLE, and C-step
///////////////////////////////////////////////////////////////////////////////////////
class Matrix_Est {
  public:
    arma::mat mu;
    arma::mat cov_row;
    arma::mat cov_col;
    arma::mat cov_row_inv;
    arma::mat cov_col_inv;

    void init(const int n, const int p, const int q) {
      mu.zeros(p,q);
      cov_row.eye(p,p);
      cov_col.eye(q,q);
      cov_row_inv.eye(p,p);
      cov_col_inv.eye(q,q);
    }
};


class MLE_res {
  public:
    Matrix_Est est;
    double norm;
    int iterations;
};

class Cstep_res {
  public:
    Matrix_Est est;
    double det;
    arma::uvec h_subset;
    arma::vec dets;
    int iterations;
};

///////////////////////////////////////////////////////////////////////////////////////
// Computation of rowwise and columnwise covariance matrix using MLE
///////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat MLErow(const arma::cube& X_std, const arma::mat& cov_col_inv) {
  // const int nthreads = std::thread::hardware_concurrency();
  const int n = X_std.n_slices;
  const int p = X_std.n_rows;
  const int q = X_std.n_cols;

  arma::mat cov_row = arma::zeros(p,p);
  for (int i = 0; i < n; i++){
    cov_row += X_std.slice(i) * cov_col_inv * X_std.slice(i).t();
  }
  // cov_row /= (n * q);
  cov_row = symmatu(cov_row/(n * q));
  return(cov_row);
}

// [[Rcpp::export]]
arma::mat MLEcol(const arma::cube& X_std, const arma::mat& cov_row_inv) {
  const int n = X_std.n_slices;
  const int p = X_std.n_rows;
  const int q = X_std.n_cols;

  arma::mat cov_col = arma::zeros(q,q);
  for (int i = 0; i < n; i++){
    cov_col += X_std.slice(i).t() * cov_row_inv * X_std.slice(i);
  }
  // cov_col /= (n * p);
  cov_col = symmatu(cov_col/(n * p));
  return(cov_col);
}


// [[Rcpp::export]]
double KroneckerNorm(arma::mat A, arma::mat B, arma::mat C, arma::mat D){
  double norm = trace(A * A) * trace(B * B) - 2*trace(A * C)*trace(B * D) + trace(C * C) * trace(D * D);
  return(norm);
}

///////////////////////////////////////////////////////////////////////////////////////
// Computation iterative matrix normal MLE algorithm (flip-flop)
///////////////////////////////////////////////////////////////////////////////////////

MLE_res mmleCpp(const arma::cube& X,
                const int max_iter = 100,
                const double lambda = 0,
                const bool silent = false){
  const int n = X.n_slices;
  const int p = X.n_rows;
  const int q = X.n_cols;

  const float pDq = (float)p/(float)q;
  const float qDp = (float)q/(float)p;

  int minimal_n = ceil(std::max(pDq,qDp)) + 1;
  int smallest_n = ceil(pDq + qDp) + 2;

  if(smallest_n > n && ! silent){
    warning({"At least n >= ceil(p/q + q/p) + 2 = " + std::to_string((int)smallest_n) + " samples are needed to ensure convergence of the mmcd algorithm and uniqueness of the estimators.\n" +
      "If at least n >= ceil(max(p/q + q/p)) + 1 = " + std::to_string((int)minimal_n) + " samples are provided the algorithm should return positive definite estimates without those guarantees."});
    if(minimal_n > n){
      stop("Only n = " + std::to_string((int)n) + " samples provided while at least "+ std::to_string((int)minimal_n) + " = ceil(max(p/q + q/p)) + 1) are requiered");
    }
  }
  // calculate mean
  arma::mat mu = mean(X, 2);

  // center X
  const arma::cube X_std = X.each_slice() - mu;

  // initialize rowwise and columnwise covariance
  arma::mat cov_row(p, p, arma::fill::eye);
  arma::mat cov_col(q, q, arma::fill::eye);
  // initialize rowwise and columnwise covariance inverse
  arma::mat cov_row_inv;
  arma::mat cov_col_inv;

  bool conv_crit = true;
  double norm;
  double scale_factor;

  int i = 0;
  while(conv_crit){
    arma::mat cov_row_old = cov_row;
    arma::mat cov_col_old = cov_col;

    bool flag_cov_col_inv = inv_sympd(cov_col_inv, cov_col);
    if(!flag_cov_col_inv){
      stop("Matrix 'cov_col' is singular");
    }
    cov_row = MLErow(X_std, cov_col_inv);

    if(lambda != 0){
      cov_row = eye(p, p)*lambda + (1 - lambda)*cov_row;
    }

    bool flag_cov_row_inv = inv_sympd(cov_row_inv, cov_row);
    if(!flag_cov_row_inv){
      stop("Matrix 'cov_row' is singular");
    }
    cov_col = MLEcol(X_std, cov_row_inv);

    if(lambda != 0){
      cov_col = eye(q, q)*lambda + (1 - lambda)*cov_col;
    }

    scale_factor = cov_col(0,0);
    cov_row *= scale_factor;
    cov_col /= scale_factor;

    norm = KroneckerNorm(cov_row_old, cov_col_old, cov_row, cov_col);
    conv_crit = norm > 0.001 && i < max_iter;
    i++;
  }
  bool flag_cov_row_inv = inv_sympd(cov_row_inv, cov_row);
  if(!flag_cov_row_inv){
    stop("Matrix 'cov_row' is singular");
  }
  bool flag_cov_col_inv = inv_sympd(cov_col_inv, cov_col);
  if(!flag_cov_col_inv){
    stop("Matrix 'cov_col' is singular");
  }

  // bool flag_cov_row_inv = true;
  // if(cov_row.is_symmetric()){
  //   flag_cov_row_inv = inv_sympd(cov_row_inv, cov_row);
  // } else{
  //   flag_cov_row_inv = inv(cov_row_inv, cov_row);
  // }
  //
  // if(!flag_cov_row_inv){
  //   stop("Matrix 'cov_row' is singular");
  // }
  //
  // bool flag_cov_col_inv = true;
  // if(cov_col.is_symmetric()){
  //   flag_cov_col_inv = inv_sympd(cov_col_inv, cov_col);
  // } else{
  //   flag_cov_col_inv = inv_sympd(cov_col_inv, cov_col);
  // }
  // if(!flag_cov_col_inv){
  //   stop("Matrix 'cov_col' is singular");
  // }
  //cov_row_inv /= scale_factor;
  //cov_col_inv *= scale_factor;

  MLE_res res;
  res.est.mu = mu;
  res.est.cov_row = cov_row;
  res.est.cov_col = cov_col;
  res.est.cov_row_inv = cov_row_inv;
  res.est.cov_col_inv = cov_col_inv;
  res.norm = norm;
  res.iterations = i;


  return(res);
}

//' Maximum Likelihood Estimation for Matrix Normal Distribtuion
//'
//' \code{mmle} computes the Maximum Likelihood Estimators (MLEs) for the matrix normal distribution
//' using the iterative flip-flop algorithm \insertCite{Dutilleul1999}{robustmatrix}.
//'
//' @param X a 3d array of dimension \eqn{(p,q,n)}, containing \eqn{n} matrix-variate samples of \eqn{p} rows and \eqn{q} columns in each slice.
//' @param max_iter upper limit of iterations.
//' @param lambda a smooting parameter for the rowwise and columnwise covariance matrices.
//' @param silent Logical. If FALSE (default) warnings and errors are printed.
//'
//' @return A list containing the following:
//' \item{\code{mu}}{Estimated \eqn{p \times q} mean matrix.}
//' \item{\code{cov_row}}{Estimated \eqn{p} times \eqn{p} rowwise covariance matrix.}
//' \item{\code{cov_col}}{Estimated \eqn{q} times \eqn{q} columnwise covariance matrix.}
//' \item{\code{cov_row_inv}}{Inverse of \code{cov_row}.}
//' \item{\code{cov_col_inv}}{Inverse of \code{cov_col}.}
//' \item{\code{norm}}{Forbenius norm of squared differences between covariance matrices in final iteration.}
//' \item{\code{iterations}}{Number of iterations of the mmle procedure.}
//'
//' @references
//' \insertAllCited{}
//'
//' @seealso For robust parameter estimation use \code{\link{mmcd}}.
//'
//' @export
//'
//' @examples
//' n = 1000; p = 2; q = 3
//' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
//' cov_row = matrix(c(1,0.5,0.5,1), nrow = p, ncol = p)
//' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
//' X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
//' par_mmle <- mmle(X)
// [[Rcpp::export]]
Rcpp::List mmle(const arma::cube& X,
                     const int max_iter = 100,
                     const double lambda = 0,
                     const bool silent = false){
  MLE_res MLE = mmleCpp(X, max_iter, lambda, silent);
  List res = List::create(Named("mu") = MLE.est.mu,
                          Named("cov_row") = MLE.est.cov_row,
                          Named("cov_col") = MLE.est.cov_col,
                          Named("cov_row_inv") = MLE.est.cov_row_inv,
                          Named("cov_col_inv") = MLE.est.cov_col_inv,
                          Named("norm") = MLE.norm,
                          Named("iterations") = MLE.iterations);
  return(res);
}


// [[Rcpp::export]]
double MMD(const arma::mat X,
           arma::mat mu,
           arma::mat cov_row,
           arma::mat cov_col,
           const bool inverted = false) {
  if(!inverted){
    bool flag_cov_row = inv_sympd(cov_row, cov_row);
    if(!flag_cov_row){
      stop("Matrix 'cov_row' is singular");
    }
    bool flag_cov_col = inv_sympd(cov_col, cov_col);
    if(!flag_cov_col){
      stop("Matrix 'cov_row' is singular");
    }
  }
  const arma::mat X_std = X - mu;
  double MD = trace(cov_col * X_std.t() * cov_row * X_std);
  return(MD);
}

// [[Rcpp::export]]
arma::vec TensorMMD(const arma::cube X,
                    arma::mat mu,
                    arma::mat cov_row,
                    arma::mat cov_col,
                    const bool inverted = false) {
  int n = X.n_slices;
  arma::vec MDs;
  MDs.zeros(n);

  if(!inverted){
    bool flag_cov_row = inv_sympd(cov_row, cov_row);
    if(!flag_cov_row){
      stop("Matrix 'cov_row' is singular");
    }
    bool flag_cov_col = inv_sympd(cov_col, cov_col);
    if(!flag_cov_col){
      stop("Matrix 'cov_row' is singular");
    }
  }

  const arma::cube X_std = X.each_slice() - mu;

  int i;
  for (i = 0; i < n; i++){
    MDs(i) =  trace(cov_col * X_std.slice(i).t() * cov_row * X_std.slice(i));
  }

  return(MDs);
}

Cstep_res cstepCpp(const arma::cube& X,
                   double alpha = 0.5,
                   int h_init = -1,
                   bool init = true,
                   const int max_iter = 100,
                   const int max_iter_MLE = 100,
                   const double lambda = 0,
                   const bool adapt_alpha = true,
                   Matrix_Est est = Matrix_Est(),
                   arma::uvec h_subset = arma::uvec()
                   ) {
  const int n = X.n_slices;
  const int p = X.n_rows;
  const int q = X.n_cols;

  const float pDq = (float)p/(float)q;
  const float qDp = (float)q/(float)p;

  int h = floor(alpha*(float)n);
  int h_minimal = ceil(std::max(pDq,qDp)) + 1;
  int h_smallest = ceil(pDq + qDp) + 2;



  if(alpha < 0.5){
    stop("'alpha' must be between 0.5 and 1");
  }

  // Adapting alpha to obtain maximum breakdown point for alpha = 0.5
  if(adapt_alpha){
    int h_max_bp = floor((n + h_smallest)/2);
    h = floor(2*(float)h_max_bp - (float)n + 2*((float)n-(float)h_max_bp)*alpha);
    alpha = (float)h/(float)n;
  }

  if(h_smallest > h){
    warning({"At least h >= ceil(p/q + q/p) + 2 = " + std::to_string((int)h_smallest) + " samples are needed in h-subset to ensure convergence of the mmcd algorithm and uniqueness of the estimators.\n" +
      "If at least h >= ceil(max(p/q + q/p)) + 1 = " + std::to_string((int)h_minimal) + " samples are provided the algorithm should return positive definite estimates without those guarantees."});
    if(h_minimal > h){
      stop("Only h = " + std::to_string((int)h) + " samples provided while at least "+ std::to_string((int)h_minimal) + " = ceil(max(p/q + q/p)) + 1) are requiered");
    }
  }

  // Create a vector containing elements from 1 to n
  arma::uvec population = linspace<uvec>(0, n-1, n);

  // C-Step counter
  int i = 0;
  // convergence criterion for difference of determinants
  double det_diff = 0;
  // vector to store determinants
  arma::vec dets;
  dets.zeros(100);

  if(alpha == 1){
    MLE_res MLE = mmleCpp(X, max_iter_MLE, lambda, true);
    est.mu = MLE.est.mu;
    est.cov_row = MLE.est.cov_row;
    est.cov_col = MLE.est.cov_col;
    est.cov_row_inv = MLE.est.cov_row_inv;
    est.cov_col_inv = MLE.est.cov_col_inv;
    i = 1;
  } else {
    ////////////////////////
    // Initialization //////
    ////////////////////////
    bool check_dim = ((int)est.mu.size() != p*q ||
                      (int)est.cov_row.size() != p*p ||
                      (int)est.cov_col.size() != q*q ||
                      (int)est.cov_row_inv.size() != p*p ||
                      (int)est.cov_col_inv.size() != q*q);
    if(!init && check_dim){
      warning("'init' was set to FALSE and the dimensions of the initial estimators 'est' do not match the dimension of the data 'X'.\nTo avoid an error the 'init' is set to TRUE. Estimators in 'est' will be discarded and properly initialized");
      init = TRUE;
    }
    if(init){
      // Initialize h-subset
      if(h_init < h_smallest){
        if(h >= h_smallest){
          h_init = h_smallest; // the smallest sample size needed to ensure convergence and uniqueness
        } else{
          h_init = h_minimal; // the minimal possible h that might still work ---- h will not be smaller since we checked for it above
        }
      }
      // Sample h elements without replacement
      h_subset = arma::randperm(n,h);


      // Use first h_init observations of h_subset to construct initial h-subset
      arma::uvec h_init_subset = h_subset.head(h_init);

      bool conv_crit = true;
      while(conv_crit && h_init < h){
        h_init += 1;
        try{
          MLE_res MLE = mmleCpp(X.slices(h_init_subset), max_iter_MLE, lambda, true);
          conv_crit = false;
          est = MLE.est;
          arma::vec MDs = TensorMMD(X, est.mu, est.cov_row_inv, est.cov_col_inv, true);
          arma::uvec MD_order = sort_index(MDs);
          h_subset = MD_order.head(h);
        } catch(...){
          h_init_subset = h_subset.head(h_init);
        }
      }
    }

    ////////////////////////
    ////////// Main   //////
    ////////////////////////
    bool det_diff_crit = true;
    while(det_diff_crit && i < max_iter){
      MLE_res MLE = mmleCpp(X.slices(h_subset), max_iter_MLE, lambda, true);
      est = MLE.est;
      arma::vec MDs = TensorMMD(X, est.mu, est.cov_row_inv, est.cov_col_inv, true);
      arma::uvec MD_order = sort_index(MDs);
      h_subset = MD_order.head(h);

      dets(i) = q*log_det_sympd(est.cov_row) + p*log_det_sympd(est.cov_col);
      if(i == 0){
        det_diff_crit = true;
      } else {
        det_diff = dets(i -1) - dets(i); // this is always larger or equal to zero
        det_diff_crit = det_diff  > 0.001; // while difference between determinants is larger then thresh, continue
      }
      i++;
    }
  }

  Cstep_res res;
  res.est = est;
  res.det = dets(i-1);
  res.h_subset = h_subset;
  res.dets = dets.head(i);
  res.iterations = i;

  return(res);
}



//' C-step of Matrix Minimum Covariance Determinant (MMCD) Estimator
//'
//' This function is part of the FastMMCD algorithm \insertCite{mayrhofer2024}{robustmatrix}.
//'
//' @param h_init Integer. Size of initial h-subset. If smaller than 0 (default) size is chosen automatically.
//' @param init Logical. If TRUE (default) elemental subsets are used to initialize the procedure.
//' @param max_iter upper limit of C-step iterations (default is 100)
//' @inheritParams mmcd
//'
//' @return A list containing the following:
//' \item{\code{mu}}{Estimated \eqn{p \times q} mean matrix.}
//' \item{\code{cov_row}}{Estimated \eqn{p} times \eqn{p} rowwise covariance matrix.}
//' \item{\code{cov_col}}{Estimated \eqn{q} times \eqn{q} columnwise covariance matrix.}
//' \item{\code{cov_row_inv}}{Inverse of \code{cov_row}.}
//' \item{\code{cov_col_inv}}{Inverse of \code{cov_col}.}
//' \item{\code{md}}{Squared Mahalanobis distances.}
//' \item{\code{md_raw}}{Squared Mahalanobis distances based on \emph{raw} MMCD estimators.}
//' \item{\code{det}}{Value of objective function (determinant of Kronecker product of rowwise and columnwise covariane).}
//' \item{\code{dets}}{Objective values for the final h-subsets.}
//' \item{\code{h_subset}}{Final h-subset of \emph{raw} MMCD estimators.}
//' \item{\code{iterations}}{Number of C-steps.}
//'
//' @seealso \code{\link{mmcd}}
//'
//' @export
//'
//' @examples
//' n = 1000; p = 2; q = 3
//' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
//' cov_row = matrix(c(1,0.5,0.5,1), nrow = p, ncol = p)
//' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
//' X <- rmatnorm(n = 1000, mu, cov_row, cov_col)
//' ind <- sample(1:n, 0.3*n)
//' X[,,ind] <- rmatnorm(n = length(ind), matrix(rep(10, p*q), nrow = p, ncol = q), cov_row, cov_col)
//' par_mmle <- mmle(X)
//' par_cstep <- cstep(X)
//' distances_mmle <- mmd(X, par_mmle$mu, par_mmle$cov_row, par_mmle$cov_col)
//' distances_cstep <- mmd(X, par_cstep$mu, par_cstep$cov_row, par_cstep$cov_col)
//' plot(distances_mmle, distances_cstep)
//' abline(h = qchisq(0.99, p*q), lty = 2, col = "red")
//' abline(v = qchisq(0.99, p*q), lty = 2, col = "red")
// [[Rcpp::export]]
Rcpp::List cstep(const arma::cube& X,
                       const double alpha = 0.5,
                       int h_init = -1,
                       bool init = true,
                       const int max_iter = 100,
                       const int max_iter_MLE = 100,
                       const double lambda = 0,
                       const bool adapt_alpha = true
){
  Cstep_res cstep = cstepCpp(X, alpha, h_init, init, max_iter, max_iter_MLE, lambda, adapt_alpha);
  List res = List::create(Named("mu") = cstep.est.mu,
                          Named("cov_row") = cstep.est.cov_row,
                          Named("cov_col") = cstep.est.cov_col,
                          Named("cov_row_inv") = cstep.est.cov_row_inv,
                          Named("cov_col_inv") = cstep.est.cov_col_inv,
                          Named("det") = cstep.det,
                          Named("h_subset") = cstep.h_subset + 1,
                          Named("dets") = cstep.dets,
                          Named("iterations") = cstep.iterations);
  return(res);
}

//' The Matrix Minimum Covariance Determinant (MMCD) Estimator
//'
//' \code{mmcd} computes the robust MMCD estimators of location and covariance for matrix-variate data
//' using the FastMMCD algorithm \insertCite{mayrhofer2024}{robustmatrix}.
//'
//' @param nsamp number of initial h-subsets (default is 500).
//' @param alpha numeric parameter between 0.5 (default) and 1. Controls the size \eqn{h \approx alpha * n} of the h-subset over which the determinant is minimized.
//' @param max_iter_cstep upper limit of C-step iterations (default is 100)
//' @param max_iter_MLE upper limit of MLE iterations (default is 100)
//' @param max_iter_cstep_init upper limit of C-step iterations for initial h-subsets (default is 2)
//' @param max_iter_MLE_init upper limit of MLE iterations for initial h-subsets (default is 2)
//' @param adapt_alpha Logical. If TRUE (default) alpha is adapted to take the dimension of the data into account.
//' @param reweight Logical. If TRUE (default) the reweighted MMCD estimators are computed.
//' @param scale_consistency Character. Either "quant" (default) or "mmd_med". If "quant", the consistency factor is chosen to achieve consistency under the matrix normal distribution.
//' If "mmd_med", the consistency factor is chosen based on the Mahalanobis distances of the observations.
//' @param outlier_quant numeric parameter between 0 and 1. Chi-square quantile used in the reweighting step.
//' @param nthreads Integer. If 1 (default), all computations are carried out sequentially.
//' If larger then 1, C-steps are carried out in parallel using \code{nthreads} threads.
//' If < 0, all possible threads are used.
//' @inheritParams mmle
//'
//' @return A list containing the following:
//' \item{\code{mu}}{Estimated \eqn{p \times q} mean matrix.}
//' \item{\code{cov_row}}{Estimated \eqn{p} times \eqn{p} rowwise covariance matrix.}
//' \item{\code{cov_col}}{Estimated \eqn{q} times \eqn{q} columnwise covariance matrix.}
//' \item{\code{cov_row_inv}}{Inverse of \code{cov_row}.}
//' \item{\code{cov_col_inv}}{Inverse of \code{cov_col}.}
//' \item{\code{md}}{Squared Mahalanobis distances.}
//' \item{\code{md_raw}}{Squared Mahalanobis distances based on \emph{raw} MMCD estimators.}
//' \item{\code{det}}{Value of objective function (determinant of Kronecker product of rowwise and columnwise covariane).}
//' \item{\code{alpha}}{The (adjusted) value of alpha used to determine the size of the h-subset.}
//' \item{\code{consistency_factors}}{Consistency factors for raw and reweighted MMCD estimators.}
//' \item{\code{dets}}{Objective values for the final h-subsets.}
//' \item{\code{best_i}}{ID of subset with best objective.}
//' \item{\code{h_subset}}{Final h-subset of \emph{raw} MMCD estimators.}
//' \item{\code{h_subset_reweighted}}{Final h-subset of \emph{reweighted} MMCD estimators.}
//' \item{\code{iterations}}{Number of C-steps.}
//' \item{\code{dets_init_first}}{Objective values for the \code{nsamp} initial h-subsets after \code{max_iter_cstep_init} C-steps.}
//' \item{\code{subsets_first}}{Subsets created in subsampling procedure for large \code{n}.}
//' \item{\code{dets_init_second}}{Objective values of the 10 best initial subsets after executing C-steps until convergence.}
//'
//' @details The MMCD estimators generalize the well-known Minimum Covariance Determinant (MCD)
//' \insertCite{Rousseeuw1985,Rousseeuw1999}{robustmatrix} to the matrix-variate setting.
//' It looks for the \eqn{h} observations, \eqn{h = \alpha * n}, whose covariance matrix has the smallest determinant.
//' The FastMMCD algorithm is used for computation and is described in detail in \insertCite{mayrhofer2024}{robustmatrix}.
//' NOTE: The procedure depends on \emph{random} initial subsets. Currently setting a seed is only possible if \code{nthreads = 1}.
//'
//' @references
//' \insertAllCited{}
//'
//' @seealso The \code{mmcd} algorithm uses the \code{\link{cstep}} and \code{\link{mmle}} functions.
//'
//' @export
//'
//' @examples
//' n = 1000; p = 2; q = 3
//' mu = matrix(rep(0, p*q), nrow = p, ncol = q)
//' cov_row = matrix(c(1,0.5,0.5,1), nrow = p, ncol = p)
//' cov_col = matrix(c(3,2,1,2,3,2,1,2,3), nrow = q, ncol = q)
//' X <- rmatnorm(n = n, mu, cov_row, cov_col)
//' ind <- sample(1:n, 0.3*n)
//' X[,,ind] <- rmatnorm(n = length(ind), matrix(rep(10, p*q), nrow = p, ncol = q), cov_row, cov_col)
//' par_mmle <- mmle(X)
//' par_mmcd <- mmcd(X)
//' distances_mmle <- mmd(X, par_mmle$mu, par_mmle$cov_row, par_mmle$cov_col)
//' distances_mmcd <- mmd(X, par_mmcd$mu, par_mmcd$cov_row, par_mmcd$cov_col)
//' plot(distances_mmle, distances_mmcd)
//' abline(h = qchisq(0.99, p*q), lty = 2, col = "red")
//' abline(v = qchisq(0.99, p*q), lty = 2, col = "red")
// [[Rcpp::export]]
Rcpp::List mmcd(const arma::cube& X,
                const int nsamp = 500,
                double alpha = 0.5,
                const double lambda = 0,
                const int max_iter_cstep = 100,
                const int max_iter_MLE = 100,
                const int max_iter_cstep_init = 2,
                const int max_iter_MLE_init = 2,
                const bool adapt_alpha = true,
                const bool reweight = true,
                const std::string scale_consistency = "quant",
                const double outlier_quant = 0.975,
                int nthreads = 1){

  if(alpha < 0.5){
    stop("'alpha' must be between 0.5 and 1");
  }

  const int n = X.n_slices;
  const int p = X.n_rows;
  const int q = X.n_cols;

  const float pDq = (float)p/(float)q;
  const float qDp = (float)q/(float)p;

  // int h_minimal = ceil(std::max(pDq,qDp)) + 1;
  int h_smallest = floor(pDq + qDp) + 2;
  int h = floor(alpha*(float)n);

  // Adapting alpha to obtain maximum breakdown point for alpha = 0.5
  if(adapt_alpha){
    int h_max_bp = floor((n + h_smallest)/2);
    h = floor(2*(float)h_max_bp - (float)n + 2*((float)n-(float)h_max_bp)*alpha);
    alpha = (float)h/(float)n;
  }

  if(h_smallest > h){
    stop({"At least h >= ceil(p/q + q/p) + 2 = " + std::to_string((int)h_smallest) +
      " samples are needed in h-subset to ensure convergence of the mmcd algorithm and uniqueness of the estimators.\n"});
  }


  if(nthreads < 0){
    nthreads = std::thread::hardware_concurrency();
  }


  arma::uvec population = linspace<uvec>(0, n-1, n); // Create a vector containing elements from 1 to n
  int nn; // size of subpopulation
  arma::uvec sub_population; // vector for subpopulation

  // variables to store the results of first initialization
  vector<arma::uvec> subsets_first;
  vector<Cstep_res> csteps_init_first;
  arma::uvec best_iter_init_first;
  arma::vec dets_init_first;

  // variables to store the results of second initialization
  vector<Cstep_res> csteps_init_second;
  arma::uvec best_iter_init_second;
  arma::vec dets_init_second;

  // variables to store the final results
  vector<Cstep_res> csteps;
  arma::vec dets;
  arma::uvec iterations;

  // Rcout << "Vector size after init = " << csteps_init_first.size();
  int n_subsets = 1;
  if(n > 600){
    nn = min(1500,n);
    sub_population = arma::randperm(n,nn);
    n_subsets = floor(nn/300);
    arma::uvec size_sub_populations = linspace<uvec>(0, nn-1, n_subsets + 1);
    for(int i = 0; i < n_subsets; i++){
      int first = size_sub_populations(i);
      int last = size_sub_populations(i+1)-1;
      if(i == n_subsets - 1){
        last ++;
      }
      subsets_first.push_back(sub_population.subvec(first, last));
    }
  } else{
    nn = n;
    sub_population = population;
    subsets_first.push_back(population);
  }

  ////////////////////////////////////////////////////////////////////////
  // First initialization on subsets (up to 5 subsets of size 300 to 600)
  ////////////////////////////////////////////////////////////////////////
  int n_samp_init = ceil((float)nsamp/(float)n_subsets);
  csteps_init_first.resize(n_samp_init*n_subsets);
  dets_init_first.zeros(csteps_init_first.size());

  try{
    std::vector<unsigned int> seeds(nthreads);
    #if defined(_OPENMP)
      #pragma omp parallel num_threads(nthreads)
      #pragma omp for collapse(2)
    #endif
    for (int i = 0; i < n_subsets; i++){
      for (int j = 0; j < n_samp_init; j++){
        csteps_init_first[i * n_samp_init + j] = cstepCpp(X.slices(subsets_first[i]),       // X (data)
                                                          alpha,                            // alpha (percentage of samples in h_subset)
                                                          -1,                               // h_init (set to -1 to use elemental subsets)
                                                          true,                             // init (use initialization)
                                                          max_iter_cstep_init,              // max_iter
                                                          max_iter_MLE_init,                // max_iter_init
                                                          lambda,                           // penalty for MLE
                                                          false                             // adapting alpha
                                                          );
        dets_init_first(i * n_samp_init + j) = csteps_init_first[i * n_samp_init + j].det;
      }
    }
    for (int i = 0; i < n_subsets; i++){
      arma::uvec dets_init_first_ordered = i * n_samp_init + sort_index(dets_init_first.subvec(i * n_samp_init, (i+1) * n_samp_init - 1));
      best_iter_init_first = join_cols(best_iter_init_first,dets_init_first_ordered.head(min(10,n_samp_init)));
    }
  } catch(...) {
    stop("mmcd error in first initialization");
  }

  ////////////////////////////////////////////////////////////////////////
  // Second initialization on subpopulation of a size of up to 1500
  ////////////////////////////////////////////////////////////////////////
  csteps_init_second.resize(best_iter_init_first.size());
  dets_init_second.zeros(best_iter_init_first.size());

  try{
    #if defined(_OPENMP)
      #pragma omp parallel num_threads(nthreads)
      #pragma omp for
    #endif
    for (int i = 0; i < (int)best_iter_init_first.size(); i++){
      arma::uvec h_subset_init_second;
      if(n > 600){
        int h_init_second = floor(alpha*(float)nn);
        arma::vec MDs_init_second = TensorMMD(X.slices(sub_population), csteps_init_first[best_iter_init_first(i)].est.mu, csteps_init_first[best_iter_init_first(i)].est.cov_row_inv, csteps_init_first[best_iter_init_first(i)].est.cov_col_inv, true);
        arma::uvec MD_order_init_second = sort_index(MDs_init_second);
        h_subset_init_second = MD_order_init_second.head(h_init_second);
      } else {
        h_subset_init_second = csteps_init_first[best_iter_init_first(i)].h_subset;
      }
      csteps_init_second[i] = cstepCpp(X.slices(sub_population),                            // X (data)
                                       alpha,                                               // alpha (percentage of samples in h_subset)
                                       -1,                                                  // h_init (will be ignored since we provide initial estimators in 'est')
                                       false,                                               // init (use initialization)
                                       max_iter_cstep_init,                                 // max_iter
                                       max_iter_MLE_init,                                   // max_iter_init
                                       lambda,                                              // penalty for MLE
                                       false,                                               // adapting alpha
                                       csteps_init_first[best_iter_init_first(i)].est,      // initial estimators
                                       h_subset_init_second                                 // initial h_subset
                                       );
      dets_init_second(i) = csteps_init_second[i].det;
    }

    arma::uvec dets_init_second_ordered = sort_index(dets_init_second);
    best_iter_init_second = dets_init_second_ordered.head(min(10,(int)dets_init_second.size()));
  } catch(...) {
    stop("mmcd error in second initialization");
  }

  ////////////////////////////////////////////////////////////////////////
  // Final c-steps on best subsets
  ////////////////////////////////////////////////////////////////////////
  csteps.resize(best_iter_init_second.size());
  dets.zeros(best_iter_init_second.size());
  iterations.zeros(best_iter_init_second.size());

  try{
    #if defined(_OPENMP)
      #pragma omp parallel num_threads(nthreads)
      #pragma omp for
    #endif
    for (int i = 0; i < (int)best_iter_init_second.size(); i++){
      arma::uvec h_subset_init_final;
      if(n > 600){
        arma::vec MDs_init_final = TensorMMD(X, csteps_init_second[best_iter_init_second(i)].est.mu, csteps_init_second[best_iter_init_second(i)].est.cov_row_inv, csteps_init_second[best_iter_init_second(i)].est.cov_col_inv, true);
        arma::uvec MD_order_init_final = sort_index(MDs_init_final);
        h_subset_init_final = MD_order_init_final.head(h);
      } else {
        h_subset_init_final = csteps_init_second[best_iter_init_second(i)].h_subset;
      }
      csteps[i] = cstepCpp(X,                                                           // X (data)
                           alpha,                                                       // alpha (percentage of samples in h_subset)
                           -1,                                                          // h_init (will be ignored since we provide initial estimators in 'est')
                           false,                                                       // init (use initialization)
                           max_iter_cstep,                                              // max_iter
                           max_iter_MLE,                                                // max_iter_init
                           lambda,                                                      // penalty for MLE
                           false,                                                       // adapting alpha
                           csteps_init_second[best_iter_init_second(i)].est,            // initial estimators
                           h_subset_init_final                                          // initial h_subset);
                           );
      dets(i) = csteps[i].det;
      iterations(i) = csteps[i].iterations;
    }
  } catch(...) {
    stop("mmcd error in final c-steps");
  }

  int best_i = dets.index_min();
  Matrix_Est est = csteps[best_i].est;

  // scale for consistency in case of matrix normal distribution
  arma::vec MDs = TensorMMD(X, est.mu, est.cov_row_inv, est.cov_col_inv, true);
  double scale_factor = 1;
  if(scale_consistency.compare("quant") == 0){
    scale_factor = alpha/R::pgamma(R::qchisq(alpha,p*q, true, false)/2, p*q/2 + 1, 1, true, false); //gamma distribution is the same as chi-square in this case
  } else if(scale_consistency.compare("mmd_med") == 0){
    scale_factor = median(MDs)/R::qchisq(0.5, p*q, true, false);
  }
  est.cov_row *= scale_factor; //just multiply cov_row (not cov_col) since we need to scale the Kronecker product [Omega (x) Sigma]
  est.cov_row_inv /= scale_factor;
  MDs /= scale_factor;

  arma::vec MDs_reweighted = MDs;
  double scale_factor_reweighted = 1;
  arma::uvec h_subset_reweighted;
  if(reweight){
    h_subset_reweighted = find(MDs < R::qchisq(outlier_quant, p*q, true, false));

    // ensure that the reweighted estimator has positive weights for at least as many observations as are contained in the h-subsets.
    if(h_subset_reweighted.size() >= csteps[best_i].h_subset.size()){
      MLE_res MLE = mmleCpp(X.slices(h_subset_reweighted), max_iter_MLE, lambda, true);
      est = MLE.est;

      // scale for consistency AGAIN in case of matrix normal distribution
      MDs_reweighted = TensorMMD(X, est.mu, est.cov_row_inv, est.cov_col_inv, true);
      if(scale_consistency.compare("quant") == 0){
        double alpha1 = (float)h_subset_reweighted.size()/(float)n;
        scale_factor_reweighted = alpha1/R::pgamma(R::qchisq(alpha1,p*q, true, false)/2, p*q/2 + 1, 1, true, false); //gamma distribution is the same as chi-square in this case
      } else if(scale_consistency.compare("mmd_med") == 0){
        scale_factor_reweighted = median(MDs)/R::qchisq(0.5, p*q, true, false);
      }
      est.cov_row *= scale_factor_reweighted; //just multiply cov_row (not cov_col) since we need to scale the Kronecker product [Omega (x) Sigma]
      est.cov_row_inv /= scale_factor_reweighted;
      MDs_reweighted /= scale_factor_reweighted;
    } else{
      h_subset_reweighted = csteps[best_i].h_subset;
    }
  }

  double scale_factor_diag = est.cov_col(0,0);
  est.cov_row *= scale_factor_diag;
  est.cov_col /= scale_factor_diag;
  est.cov_row_inv /= scale_factor_diag;
  est.cov_col_inv *= scale_factor_diag;

  arma::vec scale_factors = {scale_factor,scale_factor_reweighted};

  // arma::uvec MD_order = sort_index(MDs);
  List res = List::create(Named("mu") = est.mu,
                          Named("cov_row") = est.cov_row,
                          Named("cov_col") = est.cov_col,
                          Named("cov_row_inv") = est.cov_row_inv,
                          Named("cov_col_inv") = est.cov_col_inv,
                          Named("md") = MDs_reweighted,
                          Named("md_raw") = MDs,
                          Named("det") = csteps[best_i].det,
                          Named("alpha") = alpha,
                          Named("consistency_factors") = scale_factors,
                          Named("dets") = dets,
                          Named("best_i") = best_i,
                          Named("h_subset") = csteps[best_i].h_subset + 1,
                          Named("h_subset_reweighted") = h_subset_reweighted + 1,
                          Named("iterations") = iterations,
                          Named("dets_init_first") = dets_init_first,
                          Named("subsets_first") = subsets_first,
                          Named("dets_init_second") = dets_init_second
  );
  return(res);
}
