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
#' @importFrom stats quantile predict na.omit
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

  # extracting history effects
  tau_norm <- extract_history_parameter(object, "tau", link_function = exp)
  mixed_state <- extract_history_parameter(object, "mixed_state", link_function = boot::inv.logit)

  # extracting intercepts, note that the matrix order is
  # lm 1 for all participants, then lm 2 for all participants
  a <- extract_term_to_matrix(object, "a")

  # extracting history effect, same order as for "a"
  if (object$data$randomN == 1) {
    bH <- extract_replicate_term_to_matrix(object, "bH_mu")
  } else {
    bH <- extract_replicate_term_to_matrix(object, "bH_mu") +
      extract_replicate_term_to_matrix(object, "bH_sigma") *
      extract_term_to_matrix(object, "bH_rnd")
  }

  # extracting variance (if necessary)
  if (object$data$varianceN > 0) {
    sigma <-  rstan::extract(object$stanfit, pars = "sigma")[[1]]
  } else {
    # placeholder to pass to C code
    sigma <- rep(0, nrow(a))
  }

  # extracting fixed effects (if necessary)
  if (object$data$fixedN > 0) {
    bF <- extract_term_to_matrix(object, "bF")
  } else {
    bF <- matrix(rep(0, nrow(a)), ncol=1)
  }

  # predicting all samples
  predictions <- predict_samples(object$data$family,     # dimensions
                                 object$data$fixedN,
                                 object$data$randomN,
                                 object$data$lmN,
                                 object$data$istate - 1, # data
                                 object$data$duration,
                                 object$data$is_used,
                                 object$data$run_start,
                                 object$data$session_tmean,
                                 object$data$irandom - 1,
                                 object$data$fixed,
                                 tau_norm,               # history parameters
                                 mixed_state,
                                 object$data$history_starting_values,
                                 a,                      # lm parameters
                                 bH,
                                 bF,
                                 sigma)

  # raw samples
  if (!summary) {
    if (full_length) {
      return(predictions)
    } else {
      return(predictions[, !is.na(predictions[1, ])])
    }
  }

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
  if (full_length) {
    if (is.null(probs)) {
      return(predictions_summary$Predicted)
    } else {
      return(predictions_summary)
    }
  }

  # dropping rows with NA
  if (is.null(probs)) {
    return(na.omit(predictions_summary$Predicted))
  } else {
    return(na.omit(predictions_summary))
  }
}

