#' Extract Model Coefficients
#'
#' Extracts models population-level coefficients history-specific terms and
#' fixed-effect terms for every modeled distribution parameter.
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param summary Whether summary statistics should be returned instead of
#' raw sample values. Defaults to \code{TRUE}
#' @param probs The percentiles used to compute summary, defaults to 89% credible interval.
#' @param ... Unused.
#'
#' @return data.frame with values or summary
#'
#' @importFrom dplyr %>% mutate relocate bind_rows arrange
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock,
#'                       state = "State",
#'                       duration = "Duration",
#'                       fixed_effects = "Time")
#' coef(br_fit)
#' }
coef.cumhist <- function(object, summary=TRUE, probs=c(0.055, 0.945), ...){
  if (!is(object, "cumhist")) stop('The object must be of class "cumhist"')
  if (is.null(object$stanfit)) stop("The object has no fitted stan model")

    coef <-
      bistablehistory::historyef(object, summary, probs) %>%
      dplyr::mutate(Term = "History") %>%
      dplyr::relocate(.data$Term, .after= .data$DistributionParameter)

  if (object$data$fixedN > 0) {
    coef <-
      dplyr::bind_rows(coef,
                       bistablehistory::fixef(object, summary, probs)) %>%
      dplyr::arrange(.data$DistributionParameter)

  }
  coef
}


