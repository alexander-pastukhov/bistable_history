#' Computes predicted cumulative history using posterior predictive distribution.
#'
#' Computes predicted cumulative history using fitted model. This is just a wrapper
#' for \code{predict(object, summary, probs, full_length, predict_history=history_type)}.
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param history_type  \code{"1"} or \code{"2"} for cumulative history for the first or second perceptual
#' states (with indexes 1 and 2, respectively),
#' \code{"dominant"} or \code{"suppressed"} for cumulative history for states that either dominant
#' or suppressed during the following phase, \code{"difference"} for difference between suppressed
#' and dominant. See cumulative history vignette for details.
#' @param summary Whether summary statistics should be returned instead of
#' raw sample values. Defaults to \code{TRUE}
#' @param probs The percentiles used to compute summary, defaults to NULL (no CI).
#' @param full_length Only for \code{summary = TRUE}, whether the summary table should
#' include rows with no predictions. I.e., rows with mixed phases, first/last dominance
#' phase in the run, etc. See [preprocess_data()]. Defaults to \code{TRUE}.
#' @param ... Unused
#'
#' @return If \code{summary=FALSE}, a numeric matrix iterationsN x clearN.
#' If \code{summary=TRUE} but \code{probs=NULL} a vector of requested cumulative history values.
#' If \code{summary=TRUE} and \code{probs} is not \code{NULL}, a data.frame
#' with a column _"Predicted"_ (mean) and a column for each specified quantile.

#' @export
#'
#' @seealso \code{\link{fit_cumhist}}, \code{\link{predict.cumhist}}
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' history_difference_summary <- predict_history(br_fit, "difference")
#'
#' # full posterior prediction samples
#' history_difference <- predict_history(br_fit,
#'                                       "difference",
#'                                       summary = FALSE,
#'                                       full_length = TRUE)
#' }
predict_history <- function(object,
                            history_type,
                            summary = TRUE,
                            probs = NULL,
                            full_length = TRUE,
                            ...){
  predict(object, summary, probs, full_length, predict_history = history_type)
}
