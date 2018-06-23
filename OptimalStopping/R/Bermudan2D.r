
#' @title BSOption2D
#' @description Lower and upper bound to a two-dimensional stopping problem in case of geometric Brownian Motion and a plain vanilla payoff (i.e. call or put)as in Glasserman (2013).
#' @references Longstaff, Francis A., and Eduardo S. Schwartz. "Valuing American options by simulation: a simple least-squares approach." The review of financial studies 14.1 (2001): 113-147.
#' @references Glasserman, Paul. Monte Carlo methods in financial engineering. Vol. 53. Springer Science & Business Media, 2013.
#' @param N The number of time steps
#' @param dt The size of time steps
#' @param sigma The volatility
#' @param r The interest rate
#' @param div The dividend yield
#' @param strike The strike price
#' @param S0 The spot
#' @param M The number of simulations
#' @param option The type of option ("call" of "put")
#' @return A lower bound to the true value of the optimal stopping problem
#' @examples
#' #Glasserman example
#' lower <- BSOption2D(10, 1.0/3.0, 0.2, 0.05, 0.1, 100, 90, 1e5, option="call")
#' print(lower)
#' @export
BSOption2D<-function(N, dt, sigma, r, div, strike, S0, M, option="call")
{
  #Allocate memory for paths and payoffs
  paths<-list()
  payoffs<-list()
  paths[[1]]<- matrix(rep(c(S0), M*2), ncol=2, nrow=M)

  #Construct paths
  for(i in 2:N)
  {
    #Draw random numbers
    mm<-matrix(rnorm(M, mean=0, sd=1), ncol=2)

    #Antithetics
    mm<-rbind(mm,-mm)

    #Build increments
    paths[[i]]<- paths[[i-1]] * exp((r - div - 0.5 * sigma * sigma) * dt + sqrt(dt) * sigma * mm)
  }

  #Calculate payoffs
  for(i in 1:N)
  {
    if(option=="call")
      payoffs[[i]]<- exp(-dt*(i-1) * r) * pmax( apply(FUN=max, X=paths[[i]], MARGIN=1) - strike, 0)

    if(option=="put")
      payoffs[[i]]<- exp(-dt*(i-1) * r) * pmax( strike - apply(FUN=max, X=paths[[i]], MARGIN=1), 0)
  }

  return (LongstaffSchwartz(paths,payoffs, verbose=FALSE, onlyUseInTheMoneyPaths = FALSE, testRegressionPaths=TRUE)[[1]])
}
