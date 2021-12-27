#' Evaluates values for a fixed history parameter
#'
#' Expects either a single value within a valid range or
#' \code{randomN} values.
#'
#' @param param_name Name of the parameter.
#' @param param_value A single value or \code{randomN} numeric values.
#' @param randomN Number of levels for the random variable.
#' @param upperLimit Upper limit for a valid \code{param_value}.
#'
#' @return A numeric vector \code{randomN} long.
#'
#' @export
#' @keywords internal
#'
#' @examples
#' check_fixed_history_parameter("tau", 1, 10, Inf)
check_fixed_history_parameter <- function(param_name, param_value, randomN, upperLimit){
  if (!is.numeric(param_value)) stop(sprintf("%s must be numeric", param_name))
  if (!is.atomic(param_value)) stop(sprintf("%s must be an atomic vector", param_name))

  if (length(randomN) != 1) stop(sprintf("Number of random level for %s must be a scalar", param_name))
  if (as.integer(randomN) != randomN) stop(sprintf("Number of random level for %s must be an integer", param_name))
  if (randomN < 1) stop(sprintf("Number of random level for %s must be positive", param_name))

  if (length(param_value) != 1 && length(param_value) != randomN) {
    stop(sprintf("Number ofvalues for %s must be either 1 or %d (number of random clusters)", param_name, randomN))
  }
  if (any(param_value < 0)) stop(sprintf("%s must be non-negative", param_name))
  if (any(param_value > upperLimit)) stop(sprintf("%s must be less or equal to %g", param_name, upperLimit))

  if (length(param_value) != randomN) param_value <- rep(param_value, randomN)
  param_value
}




