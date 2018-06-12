library(testthat)
library(HTdescR)

source("./R/htdesc.R")

test_results <- test_dir("./tests/testthat", reporter = "summary")
