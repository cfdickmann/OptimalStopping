
library(testthat)
library(OptimalStopping)

test_that("1DPut",
{
  actual<-BSOption1D(51, 0.02, 0.2, 0.06, 0, 40, 36, 1e5, payoff="put")
  expect_equal(actual, 4.47, tolerance=0.01)
})


test_that("2DCall",
{
  actual<-BSOption2D(10,0.3333333,0.2,0.05,0.1,100,90,1e5,"call")
  expect_equal(actual, 8.08, tolerance=0.06)
})
