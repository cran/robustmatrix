#' Probability of obtaining at least one clean h-subset in the \code{\link{mmcd}} function.
#'
#' @param p number of rows.
#' @param q number of columns.
#' @param n_subsets number of elemental h-substs (default is 500).
#' @param contamination level of contamination (default is 0.5).
#'
#' @return Probability of obtaining at least one clean h-subset in the \code{\link{mmcd}} function.
#'
#' @export
clean_prob_mmcd <- function(p, q, n_subsets = 500, contamination = 0.5){
  h_subset_size = floor(p/q + q/p) + 2
  1 - (1-(1-contamination)^(h_subset_size))^n_subsets
}

#' Number of subsets that are required to obtain at least one clean h-subset in the \code{\link{mmcd}} function with probability \code{prob}.
#'
#' @param p number of rows.
#' @param q number of columns.
#' @param prob probability (default is 0.99).
#' @param contamination level of contamination (default is 0.5).
#'
#' @return Number of subsets that are required to obtain at least one clean h-subset in the \code{\link{mmcd}} function with probability \code{prob}.
#'
#' @export
n_subsets_mmcd <- function(p, q, prob = 0.99, contamination = 0.5){
  h_subset_size = floor(p/q + q/p) + 2
  c("n_subsets" = ceiling(log(1-prob)/log(1-(1-contamination)^(h_subset_size))), "d" = floor(p/q + q/p))
}


