#' East Asian consultation agreements from ATOP
#'
#' Consultation alliance agreements involving East and Central Asian
#' countries from 2010 to 2018. Sliced from the ATOP dyad-year dataset
#' (atop5_1dy). The dataset includes 18 countries: an East Asian core
#' (CHN, IND, JPN, KOR, PRK, PHL, THA, VNM) and their consultation
#' treaty partners (ARM, KAZ, KGZ, PAK, RUS, TJK, TKM, UKR, USA, UZB).
#'
#' @name atop_EA
#' @docType data
#' @usage data(atop_EA)
#' @keywords datasets
#'
#' @format A data frame with 156 rows and 3 columns:
#' \describe{
#'   \item{year}{Year of consultation agreement (2010-2018)}
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
