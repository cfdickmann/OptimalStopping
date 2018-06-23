library(testthat)
library(OptimalStopping)

test_that("LSexample",
{
  print("LSexample")
  set.seed(1)

  paths<-list()
  payoffs<-list()
  M<-10000
  N<-51
  dt<-0.02
  r<-0.06
  div<-0
  sigma<-0.2
  strike<-40
  S0<-36

  #Simulating paths
  paths[[1]]<- matrix(rep(c(S0), M), ncol=1, nrow=M)

  for(i in 2:N)
  {
    #Antithetics variates
    mm<-rnorm(M/2, mean=0, sd=1)
    mm<-rbind(mm, -mm)

    #Construct increments
    paths[[i]]<- paths[[i-1]] * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm,ncol=1))
  }

  #Payoffs of a put option
  for(i in 1:N)
  {
    payoffs[[i]]<- exp(-(i-1) * dt * r)*pmax(strike-paths[[i]], 0)
  }

  actual<-LongstaffSchwartz(paths, payoffs)
  print(actual)
  expect_equal(actual, 4.47, tolerance=0.02)
  print("")
})

# test_that("AndersenBroadie",
#           {
#             print("AndersenBroadie")
#             set.seed(1)
#
#             paths<-list()
#             payoffs<-list()
#             M<-10000
#             N<-51
#             dt<-0.02
#             r<-0.06
#             div<-0
#             sigma<-0.2
#             strike<-40
#             S0<-36
#
#             #Simulating paths
#             paths[[1]]<- matrix(rep(c(S0), M), ncol=1, nrow=M)
#
#             for(i in 2:N)
#             {
#               #Antithetics variates
#               mm<-rnorm(M/2, mean=0, sd=1)
#               mm<-rbind(mm, -mm)
#
#               #Construct increments
#               paths[[i]]<- paths[[i-1]] * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm,ncol=1))
#             }
#
#             #Payoffs of a put option
#             for(i in 1:N)
#             {
#               payoffs[[i]]<- exp(-(i-1) * dt * r)*pmax(strike-paths[[i]], 0)
#             }
#
#             actual<-AndersenBroadie()
#             print(actual[[2]])
#             expect_equal(actual[[2]], 4.47, tolerance=0.02)
#             print("")
# })

test_that("1DPut1ex",
{
  set.seed(1)
  print("1DPut1ex")
  actual<-BSOption1D(51, 0.02, 0.2, 0.06, 0, 40, 36, 1e5, option="put")
  print(actual[[1]])
  print(actual[[2]])
  expect_equal(actual[[1]], 4.47, tolerance=0.02)
  expect_equal(actual[[2]], 4.47, tolerance=0.02)
  print("")
})

#
# test_that("1DPut2ex",
# {
#   print("1DPut2ex")
#   actual<-BSOption1D(51, 0.04, 0.4, 0.06, 0, 40, 36, 1e5, option="put")
#   print(actual[[1]])
#   print(actual[[2]])
#   expect_equal(actual[[1]], 8.508, tolerance=0.2)
#   expect_equal(actual[[2]], 8.508, tolerance=0.2)
#   print("")
# })

#
# test_that("1DPut3ex",
# {
#   print("1DPut3ex")
#   actual<-BSOption1D(11, 0.02, 0.2, 0.06, 0, 40, 44, 1e5, option="put")
#   print(actual[[1]])
#   print(actual[[2]])
#   expect_equal(actual[[1]], 1.11, tolerance=0.2)
#   expect_equal(actual[[2]], 1.11, tolerance=0.2)
#   print("")
# })

#test_that("1DPut4ex",
#{
#  print("1DPut4ex")
#  actual<-BSOption1D(51, 0.04, 0.4, 0.06, 0, 40, 44, 1e5, option="put")
#  print(actual[[1]])
#  print(actual[[2]])
#  expect_equal(actual[[1]], 5.647, tolerance=0.2)
#  expect_equal(actual[[2]], 5.647, tolerance=0.2)
#  print("")
#})

test_that("1DPut1exSparse",
{
  set.seed(1)
  print("1DPut1exSparse")
  actual<-BSOption1D(11, 0.04, 0.4, 0.06, 0, 40, 44, 1e5, option="put")
  print(actual[[1]])
  print(actual[[2]])
  expect_equal(actual[[1]], actual[[2]], tolerance=0.04)
  print("")
})

test_that("2DCall",
{
  set.seed(1)
  print("2DCall")
  actual<-BSOption2D(10, 1.0/3.0, 0.2, 0.05, 0.1, 100, 90, 1e5, option="call")
  print(actual)
  expect_equal(actual, 8.08, tolerance=0.06)
})
