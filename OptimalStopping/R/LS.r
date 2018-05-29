
#' Longstaff Schwartz
#' Longstaff, Francis A., and Eduardo S. Schwartz. "Implementation of the Longstaff-Schwartz interest rate model." The Journal of Fixed Income 3.2 (1993): 7-14.
#' @param x A number
#' @param y A number
#' @return The sum of \code{x} and \code{y}
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
LongstaffSchwartz<-function(paths, payoffs, onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=FALSE)
{
  N <- length(payoffs)
  M <- length(payoffs[[1]])

  contImp<-payoffs[[N]]

  reg<-list()

  for(i in ((N-1):1))
  {
    if(verbose)
        print(i)

    test <- (1:M) * ((1:M) %% 2==0)
    test<-test[!test %in% c(0)]

    train <- (1:M) * ((1:M) %% 2==1)
    train<-train[!train %in% c(0)]

    ind<-(payoffs[[i]]>0)*(1:M) * ((1:M) %% 2==1)
    ind<-ind[!ind %in% c(0)]
    if(length(ind)==0)
      ind<-train

    pp<-paths[[i]]

    xx<-matrix(payoffs[[i]], ncol=1)

    for(k in (1:ncol(paths[[i]])))
      xx<-cbind(xx,paths[[i]][,k])

    for(k in (1 : ncol(paths[[i]])))
      for(j in (k : ncol(paths[[i]])))
        xx<-cbind(xx,(paths[[i]][ , k]) * (paths[[i]][ , j]))

    for(k in (1 : ncol(paths[[i]])))
      for(j in (k : ncol(paths[[i]])))
        for(r in (j : ncol(paths[[i]])))
          xx<-cbind(xx,(paths[[i]][ , k]) * (paths[[i]][ , j])* (paths[[i]][ , r]))

        xframe<-data.frame(x1=xx[,1], x2=xx[,2], x3=xx[,3], x4=xx[,4])

    if(onlyUseInTheMoneyPaths)
      reg[[i]]<-lm(contImp ~ ., xframe, subset=ind)
    else
      reg[[i]]<-lm(contImp ~ ., xframe, subset=train)

    #Predict contImp for next (proceeding) time step
    if(i>1) continuationValues<-predict(reg[[i]], xframe)
      else continuationValues<-mean(contImp)

    if(testRegressionPaths)
      contImp<-(payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp
    else
    {
      contImp[test] <-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp)[test]
      contImp[train]<-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*continuationValues)[train]
    }
  }

  results<-list()
  results[[1]]<-mean(contImp[test])
  results[[2]]<-reg
  return (results)
}
