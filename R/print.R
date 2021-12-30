#' Prints out cumhit object
#'
#' @param x A [cumhist][cumhist-class()] object
#' @param ... Unused
#' @return Nothing, console output only.
#'
#' @importFrom glue glue
#' @export
#'
#' @examples
#' \donttest{
#' br_fit <- fit_cumhist(br_singleblock, state="State", duration="Duration", fixed_effects="Time")
#' br_fit
#' }

print.cumhist <- function(x, ...){
  cat('Call: ')
  print(x$Call)
  cat('\n')

  cat(glue::glue('Family: {x$family}\n\n'))
}

