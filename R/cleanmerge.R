#' Cleaning the dataframe created by htbatch
#'
#' Replaces dataframe NA values with zero and condenses htbatch returned dataframe to only relevant descriptor value columns for modeling.
#'
#'
#' @param qsardataframe returned dataframe from htbatch
#' @param ... inherit arguments
#'
#' @importFrom tidyr replace_na
#'
#' @examples ## Not run: cleanmerge(qsardataframe)
#' @export
cleanmerge <- function(qsardataframe, ...) {
  #Replace na values
  qsardataframe <- tidyr::replace_na(qsardataframe, replace = list(r1.meta1.value = 0, r1.meta2.value = 0, r1.ortho1.value = 0, r1.ortho2.value = 0, r1.para1.value = 0, r2.meta1.value = 0, r2.meta2.value = 0, r2.ortho1.value = 0, r2.ortho2.value = 0, r2.para1.value = 0))

  #Initialize iterator
  i = 1
  n <- nrow (qsardataframe)

  for (i in 1:n) {
    if (is.na(qsardataframe$r1.meta1.smiles[i]) | is.na(qsardataframe$r1.ortho1.smiles[i]) | is.na(qsardataframe$r1.para1.smiles[i]) == FALSE) {
      #Sum up meta/para/ortho values into hammett value if available
      qsardataframe$r1.hammett.value[i] <- sum(qsardataframe$r1.meta1.value[i], qsardataframe$r1.meta2.value[i], qsardataframe$r1.ortho1.value[i], qsardataframe$r1.ortho2.value[i], qsardataframe$r1.para1.value[i])

    } else {
      qsardataframe$r1.hammett.value[i] = 0
    }

    if (is.na(qsardataframe$r2.meta1.smiles[i]) | is.na(qsardataframe$r2.ortho1.smiles[i]) | is.na(qsardataframe$r2.para1.smiles[i]) == FALSE) {
      qsardataframe$r2.hammett.value[i] <- sum(qsardataframe$r2.meta1.value[i], qsardataframe$r2.meta2.value[i], qsardataframe$r2.ortho1.value[i], qsardataframe$r2.ortho2.value[i], qsardataframe$r2.para1.value[i])

    } else {
      qsardataframe$r2.hammett.value[i] = 0
    }
  }
  qsardataframe <- cbind(qsardataframe[5],qsardataframe[55], qsardataframe[56], qsardataframe[63], qsardataframe[64], qsardataframe[65], qsardataframe[72])
  as.data.frame(qsardataframe)
  return(qsardataframe)
}
