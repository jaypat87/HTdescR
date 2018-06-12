context("HTdescR")
library(HTdescR)

unit_test <- test_that("htdesc function", {
  t <- htdesc(smile = '*Cl', HT.type = 'taft', sigma.selection = 'A')
  expect_identical(t[[1]],1)
  expect_equal(t[[2]],3)
  expect_equal(t[[3]],'Cl*')
  expect_identical(t[[4]],2.36)
})

unit_test
