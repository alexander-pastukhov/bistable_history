#' Computes predicted dominance phase durations using posterior predictive distribution.
#'
#' Computes predicted dominance phase durations using fitted model. Returns predicted
#' values only for the dominance phases that were marked for use. I.e., excluding first
#' and last dominance phases, mixed phases, etc. See [preprocess_data()].
#'
#' @param object An object of class [cumhist][cumhist-class()]
#' @param summary Whether summary statistics should be returned instead of
#' raw sample values. Defaults to \code{TRUE}
#' @param probs The percentiles used to compute summary, defaults to NULL (no CI).
#' @param full_length Only for \code{summary = TRUE}, whether the summary table should
#' include rows with no predictions. I.e., rows with mixed phases, first/last dominance
#' phase in the run, etc. See [preprocess_data()]. Defaults to \code{TRUE}.
#' @param ... Unused
#'
#' @return If \code{summary=FALSE}, a numeric matrix iterationsN x clearN.
#' If \code{summary=TRUE} but \code{probs=NULL} a vector of mean predicted durations.
#' If \code{summary=TRUE} and \code{probs} is not \code{NULL}, a data.frame
#' with a column _"Predicted"_ (mean) and a column for each specified quantile.
#'
#' @importFrom dplyr bind_cols
#' @importFrom rlang .data
#' @importFrom rstan extract
#' @importFrom stats quantile predict
#' @importFrom tibble tibble as_tibble
#'
#' @method predict cumhist
#' @export
#'
#' @seealso \code{\link{fit_cumhist}}
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' predict(br_fit)
#'
#' # full posterior prediction samples
#' predictions_samples <- predict(br_fit, summary=FALSE)
#' }
predict.cumhist <-  function(object, summary = TRUE, probs = NULL, full_length = TRUE, ...) {
  if (is.null(object$stanfit)) stop("The object has no fitted stan model")

  # extracting parameters
  lm_params <- rstan::extract(object$stanfit, pars="lm_param")$lm_param

  if (object$family == "gamma") {
    predictions <- exp(lm_params[, 1, ]) * exp(lm_params[, 2, ])
  } else if (object$family == "lognormal") {
    sigma <- rstan::extract(object$stanfit, pars="sigma")$sigma
    predictions <- exp(exp(lm_params[, 1, ]) + sigma / 2)
  } else if (object$family == "normal") {
    predictions <- lm_params[, 1, ]
  }

  # raw samples
  if (!summary) return(predictions)

  # means
  predictions_summary <- tibble::tibble(Predicted = apply(as.matrix(predictions), MARGIN=2, FUN=mean))

  # full summary
  if (!is.null(probs)) {
    predictions_summary <-
      dplyr::bind_cols(predictions_summary,
                       tibble::as_tibble(t(apply(as.matrix(predictions),
                                         MARGIN = 2,
                                         FUN = quantile,
                                         probs = probs))))
  }

  # to we need the full length?
  if (!full_length) {
    if (is.null(probs)) {
      return(predictions_summary$Predicted)
    } else {
      return(predictions_summary)
    }
  }

  full_length_predictions <-
    tibble(is_used = object$data$is_used) %>%
    group_by(.data$is_used) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    mutate(id = ifelse(.data$is_used, .data$id, NA)) %>%
    left_join(predictions_summary %>%
                mutate(id = row_number()),
              by = "id") %>%
    select(-c("is_used", "id"))

  if (is.null(probs)) {
    return(full_length_predictions$Predicted)
  } else {
    return(full_length_predictions)
  }
}

