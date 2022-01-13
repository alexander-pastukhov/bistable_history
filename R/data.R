#' Single run for binocular rivalry stimulus
#'
#' A single subject / single run dataset for binocular rivalry.
#'
#' @format A data frame with 76 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID, all rows contain _"ap"_}
#'   \item{Group}{Display, all rows contain _"BR"_}
#'   \item{Block}{Run / block index, all rows contain _1_}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Index of a perceptually dominant state, _1_, _2_ - clear perceptual state, _3_ mixed / transition phase}
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"br_singleblock"

#' Single experimental session for binocular rivalry stimulus
#'
#' A single subject / multiple runs dataset for binocular rivalry.
#'
#' @format A data frame with 76 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID, all rows contain _"ap"_}
#'   \item{Display}{Display, all rows contain _"BR"_}
#'   \item{Block}{Run / block index}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Index of a perceptually dominant state, _1_, _2_ - clear perceptual state, _3_ mixed / transition phase}
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"br_single_subject"


#' Multirun data for two participants, kinetic-depth effect display
#'
#' @format A data frame with 1186 rows and 5 variables:
#' \describe{
#'   \item{Observer}{Participant ID}
#'   \item{Block}{Run / block index}
#'   \item{State}{Factor variable for state with levels \code{-1} and \code{1} coding two clear perceptual states and \code{-2} the mixed / transition phase}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"kde_two_observers"

#' Binocular rivalry data
#'
#' Dataset on binocular rivalry for eight participants.
#'
#' @format A data frame with 3769 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID.}
#'   \item{Display}{Display, all rows contain \code{"BR"}}
#'   \item{Block}{Run / block index.}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"br"

#' Kinetic-depth effect data
#'
#' Dataset on kinetic-depth effect for eleven participants.
#'
#' @format A data frame with 38698 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID.}
#'   \item{Display}{Display, all rows contain \code{"KD"}}
#'   \item{Block}{Run / block index.}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"kde"

#' Necker cube data
#'
#' Dataset on Necker cube for five participants.
#'
#' @format A data frame with 3464 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID.}
#'   \item{Display}{Display, all rows contain \code{"NC"}}
#'   \item{Block}{Run / block index.}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
#' @source  \doi{10.1167/11.10.12}
"nc"

#' Binocular rivalry, variable contrast
#'
#' Dataset on binocular rivalry with variable but equal
#' contrast for six participants.
#'
#' @format A data frame with 4616 rows and 6 variables:
#' \describe{
#'   \item{Observer}{Participant ID.}
#'   \item{Block}{Run / block index.}
#'   \item{Contrast}{Contrast on scale from 0 to 1.}
#'   \item{Time}{Time relative to the run onset in _seconds_}
#'   \item{State}{Factor with levels \code{"Left"}, \code{"Right"} (clear states), and \code{"Mixed"}}.
#'   \item{Duration}{Duration of a dominance phase in _seconds_. Note that the duration for the last dominance phase is curtailed and, therefore, set to zero.}
#' }
"br_contrast"
