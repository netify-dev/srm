#' East Asian consultation agreements slice from ATOP
#'
#' Consultation alliance agreements involving East Asian countries from 2010 to 2018.
#' Sliced from the ATOP dyad-year dataset (atop5_1dy). Includes agreements where 
#' at least one partner is from Japan, South Korea, Philippines, Thailand, Singapore, 
#' Malaysia, Indonesia, India, China, Vietnam, and North Korea, demonstrating network 
#' dynamics in the region.
#' 
#' @name atop_EA
#' @docType data
#' @usage data(atop_EA)
#' @keywords datasets
#' 
#' @format A data frame with X rows and 3 columns:
#' \describe{
#'   \item{year}{Year of consultation agreement}
#'   \item{country1, country2}{3-letter ISO country codes for agreement partners}
#' }
#'
#' @references Leeds, Brett Ashley, Jeffrey M. Ritter, Sara McLaughlin Mitchell, 
#' and Andrew G. Long. 2002. Alliance Treaty Obligations and Provisions, 1815-1944. 
#' International Interactions 28: 237-260.
#' (\href{http://www.atopdata.org}{ATOP})
#'
#' @source \href{http://www.atopdata.org}{Alliance Treaty Obligations and Provisions (ATOP) dyad-year dataset}
#'
#' @examples
#' data(atop_EA)
#' head(atop_EA)
NULL