#' Checks for validity of values for use as normal distribution parameters.
#'
#' Should a pair of numeric values, second value should be non-zero.
#' Stops execution with an error.
#'
#' @param values Parameters for normal distribution.
#' @param parameter Name of the parameter for which the prior is defined.
#'
#' @return Logical TRUE, if none of the tests fail
#' @export
#' @keywords internal
#'
#' @examples
#' check_normal_prior(c(0, 1), "tau")
check_normal_prior <- function(values, parameter){
  if (!is.numeric(values)) stop(sprintf("prior parameters for %s is not a numeric", parameter))
  if (length(values) != 2) stop(sprintf("prior parameters for %s should be a pair of values", parameter))
  if (values[2] <= 0) stop(sprintf("variance for %s should be strictly positive", parameter))
  TRUE
}
