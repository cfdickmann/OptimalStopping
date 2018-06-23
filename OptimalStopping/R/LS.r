#' @title Longstaff & Schwartz algorithm
#' @description The model-free Longstaff & Schwartz algorithm for lower bounds to optimal stopping problems.
#' @references Longstaff, Francis A., and Eduardo S. Schwartz. "Valuing American options by simulation: a simple least-squares approach." The review of financial studies 14.1 (2001): 113-147.
#' @param paths The paths (i.e. trajectories) to base the stopping decision on. In general, this can be any information useful for the stopping decision.
#' @param payoffs The payoffs realized along the paths (first argument)
#' @param onlyUseInTheMoneyPaths Whether to only use in-the-money paths for regression
#' @param verbose Whether to print an output
#' @param testRegressionPaths Whether to test paths based on the new rule before preceeding to the previous time step
#' @param returnContinuationValues Whether to return the approximations to the continuation values as well.
#' @return A lower bound to the true value of the optimal stopping problem
#' @examples
#' #Longstaff & Schwartz algorithm
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
#' low<-LongstaffSchwartz(paths, payoffs)
#' print(low)
#' @export
LongstaffSchwartz<-function(paths, payoffs, onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=FALSE, returnContinuationValues=FALSE)
{
  #Number of time steps
  N <- length(payoffs)

  #Number of simulations
  M <- length(payoffs[[1]])

  #Allocate memory for realized approximations to continuation values
  contImp<-payoffs[[N]]

  #List of approximations to continuation values
  reg<-list()

  for(i in ((N-1):1))
  {
    if(verbose)
        print(i)

    #Testing paths are all paths with even index numbers
    test <- (1:M) * ((1:M) %% 2==0)
    test<-test[!test %in% c(0)]

    #Training paths are all paths with odd index numbers
    train <- (1:M) * ((1:M) %% 2==1)
    train<-train[!train %in% c(0)]

    #In-the-money paths
    ind<-(payoffs[[i]]>0)*(1:M) * ((1:M) %% 2==1)
    ind<-ind[!ind %in% c(0)]
    if(length(ind)==0)
      ind<-train

    if(ncol(paths[[i]])==1)
    {
      #Fifth order polynomial and payoff as basis function set for the one-dimensional case
      pp<-paths[[i]][, 1]
      xframe<-data.frame(x1=payoffs[[i]], x2=pp, x3=pp^2, x4=pp^3, x5=pp^4, x6=pp^5)
    }
    else
    {
      #Second order polynomial and payoff as basis function set for the multidimensional case
      xx<-matrix(payoffs[[i]], ncol=1)

      for(k in (1:ncol(paths[[i]])))
        xx<-cbind(xx,paths[[i]][,k])

      for(k in (1 : ncol(paths[[i]])))
        for(j in (k : ncol(paths[[i]])))
          xx<-cbind(xx,(paths[[i]][ , k]) * (paths[[i]][ , j]))

      #for(k in (1 : ncol(paths[[i]])))
      #  for(j in (k : ncol(paths[[i]])))
      #    for(r in (j : ncol(paths[[i]])))
      #      xx<-cbind(xx,(paths[[i]][ , k]) * (paths[[i]][ , j])* (paths[[i]][ , r]))

      xframe<-data.frame(xx)
    }

    #Perform regression for continuation values on paths
    if(onlyUseInTheMoneyPaths)
      reg[[i]]<-lm(contImp ~ ., xframe, subset=ind)
    else
      reg[[i]]<-lm(contImp ~ ., xframe, subset=train)

    #Predict contImp for next (proceeding) time step
    if(i>1) continuationValues<-predict(reg[[i]], xframe)
      else continuationValues<-mean(contImp)

    if(testRegressionPaths)
      #Apply estimated stopping rule on training paths
      contImp<-(payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp
    else
    {
      #Apply estimated stopping rule on testing paths
      contImp[test] <-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp)[test]

      #Apply continuation value estimation on training paths
      contImp[train]<-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*continuationValues)[train]
    }
  }

  #Return results
  if(returnContinuationValues)
  {
    #Return approximations of continuation values if desired
    results<-list()
    results[[1]]<-mean(contImp[test])
    results[[2]]<-reg
    return (results)
  }
  else
  {
    #Return lower bound
    return (mean(contImp[test]))
  }
}
