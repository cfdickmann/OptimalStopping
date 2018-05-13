#' @export
BSOption2D<-function(N,dt,sigma,r,div,strike,S0,M, payoff="call")
{
  paths<-list()
  payoffs<-list()

  paths[[1]]<- matrix(rep(c(S0), M*2), ncol=2, nrow=M)

  set.seed(1)

  for(i in 2:N)
  {
    mm<-matrix(rnorm(M, mean=0, sd=1), ncol=2)
    mm<-rbind(mm,-mm)
    paths[[i]]<- paths[[i-1]] * exp((r - div - 0.5 * sigma * sigma) * dt + sqrt(dt) * sigma * mm)
  }

  for(i in 1:N)
  {
    if(payoff=="call")
      payoffs[[i]]<- exp(-dt*(i-1) * r) * pmax( apply(FUN=max, X=paths[[i]], MARGIN=1) - strike, 0)

    if(payoff=="put")
      payoffs[[i]]<- exp(-dt*(i-1) * r) * pmax( strike - apply(FUN=max, X=paths[[i]], MARGIN=1), 0)
  }

  return (LongstaffSchwartz(paths,payoffs, verbose=FALSE, onlyUseInTheMoneyPaths = FALSE, testRegressionPaths=TRUE))
}
