#Payoff function for plain vanilla style options for use in the Andersen Broadie algorithm
ABpayoff<-function(x, strike, option)
{
  if(option=="put")
    payoffs<- pmax(strike - x, 0)

  if(option=="call")
    payoffs<- pmax(x - strike, 0)

  return(payoffs)
}

#Subsimulations for the Andersen Broadie algorithm
BS <- function(S0, MM, dt, r, div, sigma)
{
  #Draw random number
  mm<-rnorm(MM/2, mean=0, sd=1)

  #Antithetics
  mm<-rbind(mm, -mm)

  #Simulate paths
  p<- S0 * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm, ncol=1))

  return (p)
}

#' @description Andersen & Brodie algorithm (martingales from continuation values) for one-dimensional geometric Brownian Motion in case of a plain vanilla payoff
#' @title Andersen % Broadie algorithm (martingales from continuation values)
#' @references Andersen, Leif, and Mark Broadie. "Primal-dual simulation algorithm for pricing multidimensional American options." Management Science 50.9 (2004): 1222-1234.
#' @references Longstaff, Francis A., and Eduardo S. Schwartz. "Valuing American options by simulation: a simple least-squares approach." The review of financial studies 14.1 (2001): 113-147.
#' @references Glasserman, Paul. Monte Carlo methods in financial engineering. Vol. 53. Springer Science & Business Media, 2013.
#' @param paths The paths
#' @param payoffs The payoffs realized along the paths (first argument)
#' @param option Whether it is a "call" or "put" option
#' @param dt Size of time step
#' @param r Interest rate
#' @param div Dividend yield
#' @param sigma Volatility
#' @param strike The strike
#' @return A lower and an upper bound for the solution to the optimal stopping problem
#' @examples
#' #Andersen & Broadie algorithm
#' paths<-list()
#' payoffs<-list()
#' M<-10000
#' N<-51
#' dt<-0.02
#' r<-0.06
#' div<-0
#' sigma<-0.2
#' strike<-40
#' S0<-36
#' option<-"put"
#'
#' #Simulating paths
#' paths[[1]]<- matrix(rep(c(S0), M), ncol=1, nrow=M)
#'
#' for(i in 2:N)
#' {
#'   #Antithetics variates
#'   mm<-rnorm(M/2, mean=0, sd=1)
#'   mm<-rbind(mm, -mm)
#'
#'   #Construct increments
#'   paths[[i]]<- paths[[i-1]] * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm,ncol=1))
#' }
#'
#' #Payoffs of a put option
#' for(i in 1:N)
#' {
#'   payoffs[[i]]<- exp(-(i-1) * dt * r)*pmax(strike-paths[[i]], 0)
#' }
#'
#' bounds<-AndersenBroadie(paths, payoffs, option, dt, r, div, sigma, strike)
#' low<-bounds[1]
#' high<-bounds[2]
#' print(low)
#' print(high)
#' @export
AndersenBroadie<-function(paths, payoffs, option="put", dt, r, div, sigma, strike)
{
   #Perform Longstaff & Schwartz in order to get continuation value approximations
   LSresults<-LongstaffSchwartz(paths, payoffs, onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=TRUE, returnContinuationValues=TRUE)

   #Read continuation value approximations
   reg<-LSresults[[2]]

   #Number of paths
   M<-length(payoffs[[1]]) / 1000

   #Number of subsims
   MM<-1000

   #Number of time steps
   N<-length(payoffs)

   #Allocate memory for martingales
   martingales<-list()

   for(i in 1:M)
   {
      #Allocate space for martingale
      mar<-numeric(N)
      mar[1]<-0

      #Martingales from continuation values
      for(n in 1:(N-1))
      {
        #Value of trajectory
        pp<- paths[[n+1]][i]

        #Insert into continuation value approximation
        if(n == N-1)
          cont<-0
        else
          cont<-predict(reg[[n+1]], data.frame(x1=payoffs[[n+1]][i], x2=pp, x3=pp^2, x4=pp^3, x5=pp^4, x6=pp^5))

        #Approximated true value
        Vhut<-max(cont, payoffs[[n+1]][i])

        #Generate subsimulations
        subsims<-BS(paths[[n]][i][1], MM, dt, r, div, sigma)

        #Payoffs of subsimulations
        Zpay<-exp( -r * dt * (n-1)) * ABpayoff(subsims, strike, option)

        #estimate continuation values on subsimulations
        if(n==(N-1))
          Zcont<-0
        else
          Zcont<-predict(reg[[n+1]], data.frame(x1=Zpay, x2=subsims, x3=subsims^2, x4=subsims^3, x5=subsims^4, x6=subsims^5))

        #Build martingale as martingale part of approximated true value process
        ZVhut<-pmax(Zpay,Zcont)

        #Martingale increment
        delta<-Vhut-mean(ZVhut)

        #Martingale value
        mar[n+1] <- mar[n] + delta
      }

     #Store martingale
     martingales[[i]]<-mar
   }

   #Allocate memory for results of dual formulation
   ergs <- numeric(M)

   for(i in 1:M)
   {
     #Vector of payoffs
     pay <- numeric(N)
     for (k in 1:N)
       pay[k] <- payoffs[[k]][i]

     #Dual formulation
     ergs[i] <- max(pay - martingales[[i]])
   }

   results<-list()
   results[[1]]<-LSresults[[1]]
   results[[2]]<-mean(ergs)

   return (results)
}
