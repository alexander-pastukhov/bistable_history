#' Computes R-squared using Bayesian R-squared approach.
#'
#' For detail refer to:
#' Andrew Gelman, Ben Goodrich, Jonah Gabry, and Aki Vehtari (2018).
#' R-squared for Bayesian regression models. The American Statistician
#' \doi{10.1080/00031305.2018.1549100} and
#' [https://avehtari.github.io/bayes_R2/bayes_R2.html](https://avehtari.github.io/bayes_R2/bayes_R2.html)
#'
#' @name bayes_R2
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param summary Whether summary statistics should be returned instead of
#' raw sample values. Defaults to \code{TRUE}
#' @param probs The percentiles used to compute summary, defaults to 89% credible interval.
#' @param ... Unused.
#'
#' @return vector of values or a data.frame with summary
#'
#' @importFrom dplyr %>% bind_cols
#' @importFrom rlang .data
#' @importFrom rstantools bayes_R2
#' @importFrom stats quantile
#' @importFrom tibble tibble as_tibble
#'
#' @method bayes_R2 cumhist
#' @export bayes_R2
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' bayes_R2(br_fit)
#' }
bayes_R2.cumhist <- function(object, summary=TRUE, probs=c(0.055, 0.945), ...){
  if (is.null(object$stanfit)) stop("The object has no fitted stan model")

  # compute bayes R2 via rstantools
  predicted <- predict(object, summary=FALSE)
  r2s <- rstantools::bayes_R2(predicted, object$data$clear_duration)

  # returning
  if (!summary) return(r2s)

  r2_mean <- apply(as.matrix(r2s), MARGIN=2, FUN=mean)
  if (is.null(probs)) return(r2_mean)

  # full summary
  tibble::tibble(R2 = r2_mean) %>%
    dplyr::bind_cols(tibble::as_tibble(t(apply(as.matrix(r2s), MARGIN=2, FUN=quantile, probs=probs))))
}
