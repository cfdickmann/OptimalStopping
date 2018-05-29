library(testthat)
library(OptimalStopping)

test_that("1DPut",
          {
            actual<-BSOption1D(51, 0.02, 0.2, 0.06, 0, 40, 36, 1e5, option="put")
            print(actual[[1]])
            print(actual[[2]])
            expect_equal(actual[[1]], 4.47, tolerance=0.02)
            expect_equal(actual[[2]], 4.47, tolerance=0.02)
})

#test_that("1DPut2ex",
#          {
#            actual<-BSOption1D(51, 0.04, 0.2, 0.06, 0, 40, 36, 1e5, option="put")
#            print(actual[[1]])
#            print(actual[[2]])
#            expect_equal(actual[[1]], 4.82, tolerance=0.02)
#            expect_equal(actual[[2]], 4.82, tolerance=0.02)
#})

#test_that("1DPut3exex",
#{
#  actual<-BSOption1D(51, 0.02, 0.4, 0.06, 0, 40, 36, 1e5, option="put")
#  print(actual[[1]])
#  expect_equal(actual[[1]], 8.48, tolerance=0.01)
#  print(actual[[2]])
#  expect_equal(actual[[2]], 8.48, tolerance=0.01)
#})

test_that("2DCall",
{
  actual<-BSOption2D(10,0.3333333,0.2,0.05,0.1,100,90,1e5, option="call")
  print(actual)
  expect_equal(actual, 8.08, tolerance=0.2)
})
