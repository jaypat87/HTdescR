#' Cleaning the dataframe created by htbatch
#'
#' Cleaning the dataframe created by htbatch
#'
#'
library(tidyr)
cleanmerge <- function (qsardataframe, ...) {
  qsardataframe <- replace_na(qsardataframe, replace = list(r1.meta1.value = 0, r1.meta2.value = 0, r1.ortho1.value = 0, r1.ortho2.value = 0, r1.para1.value = 0, r2.meta1.value = 0, r2.meta2.value = 0, r2.ortho1.value = 0, r2.ortho2.value = 0, r2.para1.value = 0))

  return(qsardataframe)
}
