#' Extracts a term with one column per fixed or random-level into a matrix
#'
#' Extracts a 3D array for a term with "[" sample, linear-model,
#' random/fixed-effect "]" order and returns a matrix with samples as rows
#' and columns in order 1) all random/fixed effects for lm1, 2) all
#' random/fixed effects for lm2, etc.
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param term String, term name
#'
#' @return Matrix
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' a <- extract_term_to_matrix(br_fit, "a")
#' }
extract_term_to_matrix <- function(object, term) {
  term_as_array <- rstan::extract(object$stanfit, pars = term)[[1]]
  do.call(cbind, purrr::map(1:object$data$lmN, ~term_as_array[, ., ]))
}


#' Extract a term and replicates it randomN times for each linear model
#'
#' Extract a term and replicates it randomN times for each linear model.
#' Used for population mean or variance terms.
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param term String, term name
#'
#' @return Matrix
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' bH_mu <- extract_replicate_term_to_matrix(br_fit, "bH_mu")
#' }
extract_replicate_term_to_matrix <- function(object, term) {
  term_as_array <- rstan::extract(object$stanfit, pars = term)[[1]]

  do.call(cbind,
          purrr::map(1:object$data$lmN,
                     ~matrix(rep(term_as_array[, .], object$data$randomN),
                             ncol = object$data$randomN)))
}
