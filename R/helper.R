#' Expidited search for phenyl fragment sigma values
#'
#' This function is for returning values for the phenyl fragment without running htdesc.
#' @usage helper(type, sigma.select)
#'
#' @param type Matches the HT.type specified in htbatch function value
#' @param sigma.select The type of sigma to be returned; valid inputs include "A", "B", "C", "D", "E", "F", "G", "H", and "U"
#'
#' @return values from library that is inserted directly into qsardataframe.
#'
helper <- function(type, sigma.select){
  hfilter <- dplyr::filter(htdescHelper, HT.type == type, sigma.selection == Sigma.select)
  return(hfilter)
  }

