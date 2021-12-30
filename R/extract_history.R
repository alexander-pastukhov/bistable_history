#' Computes history for a fitted model
#'
#' Computes history for a fitted model, uses only mean values
#' for each history parameter. Uses values for each random cluster,
#' if \code{"random"} or \code{"1|random"} parametrisation was used.
#'
#' @param object An object of class [cumhist][cumhist-class()]
#'
#' @return A matrix of cumulative history values for each state
#'
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state = "State", duration = "Duration")
#' extract_history(br_fit)
#' }
extract_history <- function(object)
{
  if (!is(object, "cumhist")) stop('The object must be of class "cumhist"')
  if (is.null(object$stanfit)) stop("The object has no fitted stan model")

  # data as table
  data <- tibble::tibble(istate = object$data$istate,
                         duration = object$data$duration,
                         irandom= object$data$irandom,
                         run_start = object$data$run_start,
                         session_tmean = object$data$session_tmean)

  # history parameters
  tau <-  bistablehistory::check_fixed_history_parameter("tau",
                                                         history_tau(object, probs=NULL, includePopulationLevel=FALSE)$Estimate,
                                                         object$data$randomN, Inf)
  mixed_state <- bistablehistory::check_fixed_history_parameter("mixed_state",
                                                                history_mixed_state(object, probs=NULL, includePopulationLevel=FALSE)$Estimate,
                                                                object$data$randomN, Inf)
  history_init <- bistablehistory::evaluate_history_init(object$data$history_starting_values)

  # computing the history
  h <- fast_history_compute(data, tau, mixed_state, history_init)


  # adding state names
  colnames(h) <- levels(factor(object$data$state))[1:2]

  h
}
