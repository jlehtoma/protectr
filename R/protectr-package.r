#' protectr: Systematic Protected Area Design in R
#'
#' Solve Marxan-like systematic reserve design problems using integer
#' programming techniques as implemented by the Gurobi Optimizer.
#'
#' @name protectr
#' @docType package
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(".")
