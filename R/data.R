#' Simulated classroom friendship ratings
#'
#' A 12x12 matrix of simulated friendship ratings between students at
#' North Shore High, inspired by *Mean Girls*. Generated with actor,
#' partner, and unique effect components to serve as a simple
#' cross-sectional example with directed (asymmetric) ties.
#'
#' @name classroom
#' @docType data
#' @usage data(classroom)
#' @format A 12x12 numeric matrix with row and column names
#'   corresponding to character names: Cady, Aaron, Gretchen, Karen,
#'   Janis, Damian, Kevin, Glen, Shane, Trang, Regina, Norbury.
#'   Diagonal is zero.
#' @keywords datasets
#' @examples
#' data(classroom)
#' fit = srm(classroom)
#' summary(fit)
NULL

#' Simulated trade network (longitudinal)
#'
#' A named list of three 10x10 matrices representing simulated bilateral
#' trade intensity between 10 countries across three time periods
#' (2015, 2017, 2019). Useful for demonstrating longitudinal SRM analysis.
#'
#' @name trade_net
#' @docType data
#' @usage data(trade_net)
#' @format A named list of three 10x10 numeric matrices. Countries:
#'   USA, CHN, DEU, JPN, GBR, FRA, KOR, IND, BRA, CAN.
#' @keywords datasets
#' @examples
#' data(trade_net)
#' fit = srm(trade_net)
#' summary(fit)
NULL

#' Simulated Small Council pursuit data (bipartite)
#'
#' A 10x7 matrix of simulated "pursuit intensity" scores representing
#' how aggressively each Great House from *Game of Thrones* seeks
#' each Small Council position. Generated with deliberate actor
#' effects (some houses are more ambitious overall), partner effects
#' (some positions are more sought after), and unique dyadic effects
#' (house-specific affinities for particular roles). Serves as a
#' bipartite (two-mode) example where senders and receivers come
#' from different populations.
#'
#' @name small_council
#' @docType data
#' @usage data(small_council)
#' @format A 10x7 numeric matrix. Rows are Great Houses: Stark,
#'   Lannister, Targaryen, Baratheon, Tyrell, Martell, Greyjoy,
#'   Arryn, Tully, Bolton. Columns are Small Council positions:
#'   Hand, Coin, Whispers, Ships, War, Law, Faith.
#' @keywords datasets
#' @examples
#' data(small_council)
#' fit = srm(small_council)
#' summary(fit)
NULL
