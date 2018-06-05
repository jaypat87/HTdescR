#' Cleaning the dataframe created by htbatch
#'
#' Cleaning the dataframe created by htbatch
cleanmerge <- function (qsardataframe, ...) {
  i = 1
  n <- nrow (qsardataframe)

  for (i in 1:n) {
  if (is.na(qsardataframe$r2.meta1.value[i]) == TRUE) {
    qsardataframe$r2.meta1.value[i] = 0
      }
  }
  return (qsardataframe)
}

