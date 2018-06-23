#' @title BSOption1D
#' @description Lower and upper bound for a one-dimensional stopping problem in case of geometric Brownian Motion and a plain vanilla payoff (i.e. call or put) as in Longstaff & Schwartz (2001).
#' @references Andersen, Leif, and Mark Broadie. "Primal-dual simulation algorithm for pricing multidimensional American options." Management Science 50.9 (2004): 1222-1234.
#' @references Longstaff, Francis A., and Eduardo S. Schwartz. "Valuing American options by simulation: a simple least-squares approach." The review of financial studies 14.1 (2001): 113-147.
#' @references Glasserman, Paul. Monte Carlo methods in financial engineering. Vol. 53. Springer Science & Business Media, 2013.
#' @author Fabian Dickmann
#' @param N The number of time steps
#' @param dt The size of time steps
#' @param sigma The volatility
#' @param r The interest rate
#' @param div The dividend yield
#' @param strike The strike price
#' @param S0 The spot
#' @param M The number of simulations (The number of subsimulations in the ANdersen Broadie algorithm is 1000 times smaller)
#' @param option The type of option ("call" of "put")
#' @return A lower and an upper bound for the true value of the optimal stopping problem
#' @examples
#' #Longstaff & Schwartz example
#' bounds <- BSOption1D(51, 0.02, 0.2, 0.06, 0, 40, 36, 1e5, option="put")
#' lower <- bounds[1]
#' upper <- bounds[2]
#' print(lower)
#' print(upper)
#' @export
BSOption1D<-function(N, dt, sigma, r, div, strike, S0, M, option="call")
{
  #Allocate memory for paths and payoffs
  paths <- list()
  payoffs <- list()
  paths[[1]] <- matrix(rep(c(S0), M), ncol=1, nrow=M)

  #Simulating paths
  for(i in 2:N)
  {
    mm<-rnorm(M/2, mean=0, sd=1)
    mm<-rbind(mm, -mm)
    paths[[i]] <- paths[[i-1]] * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm,ncol=1))
  }

  #Moment matching
  for(i in 2:N)
  {
    #Factor to divide by
    factor <- S0*exp((i-1)*dt*(r-div))/mean(paths[[i]])

    #Adjust paths
    paths[[i]] <- paths[[i]]/factor
  }

  #Calculate payoffs
  for(i in 1:N)
  {
    if(option=="put")
      payoffs[[i]] <- exp(-(i-1) * dt * r) * pmax( strike - paths[[i]], 0)

    if(option=="call")
      payoffs[[i]] <- exp(-(i-1) * dt * r) * pmax( paths[[i]] - strike, 0)
  }

  return (AndersenBroadie(paths, payoffs, option, dt, r, div, sigma, strike))
}
