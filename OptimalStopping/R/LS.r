#' @export
LongstaffSchwartz<-function(paths, payoffs, onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=FALSE)
{
  N <- length(payoffs)
  M <- length(payoffs[[1]])

  contImp<-payoffs[[N]]

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

    pp<-paths[[i]]

    xx<-matrix(payoffs[[i]], ncol=1)

    for(k in (1:ncol(paths[[i]])))
      xx<-cbind(xx,paths[[i]][,k])

    for(k in (1 : ncol(paths[[i]])))
      for(j in (k : ncol(paths[[i]])))
        xx<-cbind(xx,(paths[[i]][ , k]) * (paths[[i]][ , j]))

    xframe<-data.frame(xx)

    if(onlyUseInTheMoneyPaths)
      reg<-lm(contImp ~ ., xframe, subset=ind)
    else
      reg<-lm(contImp ~ ., xframe, subset=train)

    #Predict contImp for next (proceeding) time step
    if(i>1) continuationValues<-predict(reg, xframe)
      else continuationValues<-mean(contImp)

    if(testRegressionPaths)
      contImp<-(payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp
    else
    {
      contImp[test] <-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*contImp)[test]
      contImp[train]<-((payoffs[[i]]>continuationValues)*payoffs[[i]] + (payoffs[[i]]<=continuationValues)*continuationValues)[train]
    }
  }

  return (mean(contImp[test]))
}
