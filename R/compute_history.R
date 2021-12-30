#' Computes cumulative history for the time-series
#'
#' Computes cumulative history for each state in the time-series.
#'
#' @param data A table with time-series.
#' @param state String, the name of the column that specifies
#' perceptual state. The column type should be a factor with
#' two or three levels (the third level is assumed to correspond to a
#' transition/mixed phase) or should be convertible to a two level
#' factor (as it would be impossible to infer the identity of transition/
#' mixed phase).
#' @param duration String, name of the column with duration of individual
#' perceptual dominance phases. Optional, you can specify \code{onset}
#' instead.
#' @param onset String, name of the column with onsets of the perceptual
#' dominance states. Optional, used to compute duration of the dominance
#' phases, if these are not provided explicitly via \code{duration}
#' parameter.
#' @param random_effect String, name of the column that identifies random effect,
#' e.g. individual participants, stimuli for a single participant, etc.
#' If omitted, no random effect is assumed. If specified and
#' there is more than one level (participant, stimulus, etc.), it is used
#' in a hierarchical model.
#' @param session String, name of the column that identifies unique
#' experimental session for which a mean dominance phase duration will
#' be computed (see \code{norm_tau} parameter). Code assumes that session
#' IDs are different within a participant but can be the same between them.
#' If omitted, a single mean dominance duration based on the entire time series
#' is used.
#' @param run String, name of the column that identifies unique runs/blocks.
#' If omitted, the data is assumed to belong to a single time series. Code
#' assumes that run IDs are different within an experimental session but
#' can be the same between the session. E.g. session A, runs 1, 2, 3.. and
#' session B, runs 1, 2, 3 but not session A, runs 1, 2, 1.
#' @param tau Time constant of exponential growth/decay
#' normalized to the mean duration of clear percepts within each \code{session}.
#' Can be 1) a single positive number (>0) that is used for all participants and runs,
#' 2) \code{NULL} (default) -  a _single_ value will be fitted for all participants and runs,
#' 3) \code{"random"} - an independent tau is fitted for each random cluster,
#' 4) \code{"1|random"}- a tau for a random cluster
#' is sampled from a population distribution, i.e., pooled parameter values via
#' a multilevel model.
#' @param mixed_state Specifies an activation level during
#' transition/mixed phases (state #3, see \code{state}). Either a single
#' number (range 0..1) that will be used as a fixed level or a vector
#' of two numbers \code{c(mu, kappa)} that specifies, correspondingly, mean
#' (range 0..1) and precision (>0) of beta proportion distribution, it
#' should be sampled from. Defaults to a fixed value of \code{0.5}.
#' @param history_init Initial value for cumulative history computation. Either
#' a numeric scalar in 0..1 range or a vector of two numbers in 0..1 range.
#' In the latter case, two histories will start at different levels.
#'
#' @return A matrix \code{nrow(data)} × 2 with computed history values
#' @export
#'
#' @examples
#' df <- compute_history(br_singleblock, state = "State",
#'                       duration = "Duration", tau = 1,
#'                       mixed_state = 0.5, history_init = 0)
compute_history <- function(data,
                            state,
                            duration=NULL,
                            onset=NULL,
                            random_effect=NULL,
                            session=NULL,
                            run=NULL,
                            tau=1,
                            mixed_state=0.5,
                            history_init=0){

  # preprocess data for history computation
  data <- bistablehistory::preprocess_data(data, state, duration, onset, random_effect, session, run)
  randomN <- length(unique(data$random))

  # checking history computation parameters
  tau <- bistablehistory::check_fixed_history_parameter("tau", tau, randomN, Inf)
  mixed_state <- bistablehistory::check_fixed_history_parameter("mixed_state", mixed_state, randomN, 1)
  history_init <- bistablehistory::evaluate_history_init(history_init)

  # computing the history
  h <- fast_history_compute(data, tau, mixed_state, history_init)

  # adding state names
  colnames(h) <- levels(factor(data$State))[1:2]

  h
}

