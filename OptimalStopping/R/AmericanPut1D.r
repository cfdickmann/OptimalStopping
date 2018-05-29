#' @export
BSOption1D<-function(N, dt, sigma, r, div, strike, S0, M, option="call")
{
  paths<-list()
  payoffs<-list()

  paths[[1]]<- matrix(rep(c(S0), M), ncol=1, nrow=M)

  set.seed(1)

  #Simulating paths
  for(i in 2:N)
  {
    mm<-rnorm(M/2, mean=0, sd=1)
    mm<-rbind(mm, -mm)
    paths[[i]]<- paths[[i-1]] * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm,ncol=1))
  }

  #Moment matching
  for(i in 2:N)
  {
    factor<-S0*exp((i-1)*dt*(r-div))/mean(paths[[i]])
    paths[[i]]<- paths[[i]]/factor
  }

  for(i in 1:N)
  {
    if(option=="put")
      payoffs[[i]]<- exp(-(i-1) * dt * r)*pmax(strike-paths[[i]], 0)

    if(option=="call")
      payoffs[[i]]<- exp(-(i-1) * dt * r)*pmax(paths[[i]]-strike, 0)
  }

  return (AndersenBroadie(paths, payoffs, option,  onlyUseInTheMoneyPaths = TRUE, verbose=FALSE, testRegressionPaths=TRUE, dt, r, div, sigma, strike))
     }
