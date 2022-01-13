#' Extracts a history parameter as a matrix
#'
#' Extracts a history parameter as a matrix with
#' \code{samplesN} rows and randomN (found in \code{object$data$randomN})
#' columns.
#'
#' @param object A [cumhist][cumhist-class()] object
#' @param param_name String, a name of the parameter
#' @param samplesN Number of samples, if NULL is computed from rstan (but it is cheaper to do this once).
#' @param link_function A link function to use (exp or inv.logit) or `NULL` for identity.
#'
#' @return Matrix with \code{samplesN} rows and randomN
#' (found in \code{object$data$randomN}) columns
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration")
#' extract_history_parameter(br_fit, "tau", link_function = exp)
#' }
extract_history_parameter <- function(object, param_name, samplesN=NULL, link_function=NULL) {
  # figuring out number of samples
  if (is.null(samplesN)) {
    sampler_params <- rstan::get_sampler_params(object$stanfit, inc_warmup = FALSE)
    # samplesN = number of samples in one chain * number of chains
    samplesN <- dim(sampler_params[[1]])[1] * length(sampler_params)
  }

  # extracting parameter based on the used option
  option <- object$data[[paste0(param_name, "_option")]]
  if (option == 1) {
    # fixed value
    values <- matrix(rep(object$data[[paste0('fixed_', param_name)]],
                         samplesN * object$data$randomN),
                     ncol = object$data$randomN)
  } else if (option == 2) {
    # single fitted value
    values <- matrix(rep(rstan::extract(object$stanfit, pars=paste0(param_name, "_mu"))[[1]],
                         object$data$randomN),
                     ncol = object$data$randomN)
  } else if (option == 3) {
    # independent random-effect values
    values <- rstan::extract(object$stanfit, pars=paste0(param_name, "_rnd"))[[1]]
  } else if (option == 4) {
    # pooled random-effect values
    param_mu <- matrix(rep(rstan::extract(object$stanfit, pars=paste0(param_name, "_mu"))[[1]],
                           object$data$randomN),
                       ncol = object$data$randomN)
    param_sigma <- matrix(rep(rstan::extract(object$stanfit, pars=paste0(param_name, "_sigma"))[[1]],
                              object$data$randomN),
                          ncol = object$data$randomN)
    param_rnd <- rstan::extract(object$stanfit, pars=paste0(param_name, "_rnd"))[[1]]
    values <- param_mu + param_sigma * param_rnd
  }

  # applying a link function (if necessary)
  if (!is.null(link_function) & option != 1) values <- link_function(values)

  values
}
