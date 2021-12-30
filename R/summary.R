#' Summary for a cumhist object
#'
#' @param object A [cumhist][cumhist-class()] object
#' @param ... Unused
#' @return Nothing, console output only.
#'
#' @importFrom glue glue
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration", fixed_effects="Time")
#' summary(br_fit)
#' }

summary.cumhist <- function(object, ...){
  cat('Call: ')
  print(object$Call)
  cat('\n')

  cat(glue::glue('Family: {object$family}\n\n\n'))

  cat("History parameters:\n")
  cat(glue::glue("    tau = {round(history_tau(object)$Estimate, 2)}\n\n"))
  cat(glue::glue("    mixed state = {round(history_mixed_state(object)$Estimate, 2)}\n\n"))
  cat(glue::glue("    history mix = {round(history_mix(object)$Estimate, 2)}\n\n"))

  cat("\nLinear model:\n")
  print(coef(object))
}


