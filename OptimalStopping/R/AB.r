

ABpayoff<-function(x, strike, option)
{
  if(option=="put")
    payoffs<- pmax(strike-x, 0)

  if(option=="call")
    payoffs<- pmax(x-strike, 0)

  return(payoffs)
}

BS <- function(S0, MM, dt, r, div, sigma)
{
  mm<-rnorm(MM/2, mean=0, sd=1)
  mm<-rbind(mm, -mm)
  p<- S0 * exp(dt * (r - div -0.5 * sigma * sigma) + sqrt(dt) * sigma * matrix(mm, ncol=1))
  return (p)
}

 #' @export
 AndersenBroadie<-function(paths, payoffs, option="put", onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=FALSE, dt, r, div, sigma, strike)
 {
   LSresults<-LongstaffSchwartz(paths,payoffs,onlyUseInTheMoneyPaths=FALSE, verbose=FALSE, testRegressionPaths=TRUE)
   reg<-LSresults[[2]]

   M<-length(payoffs[[1]])/1000
   MM<-1000 #Number of subsims
   N<-length(payoffs)

   martingales<-list()

   for(i in 1:M)
   {
      mar<-numeric(N)
      mar[1]<-0

      for(n in 1:(N-1))
      {
        xx<-matrix(payoffs[[n+1]][i], ncol=1)
        xx<-cbind(xx, paths[[n+1]][i])
        xx<-cbind(xx, paths[[n+1]][i] * paths[[n+1]][i])
        xx<-cbind(xx, paths[[n+1]][i] * paths[[n+1]][i]* paths[[n+1]][i])

        if(n==(N-1))
          cont<-0
        else
          cont<-predict(reg[[n+1]], data.frame(x1=xx[,1], x2=xx[,2], x3=xx[,3], x4=xx[,4]))
        Vhut<-max(cont, payoffs[[n+1]][i])

        subsims<-BS(paths[[n]][i][1], MM, dt, r, div, sigma)
        Zpay<-exp(-r*dt*(n-1))*ABpayoff(subsims, strike, option)
        Zxx<-matrix(Zpay, ncol=1)
        Zxx<-cbind(Zxx, subsims)
        Zxx<-cbind(Zxx, subsims * subsims)
        Zxx<-cbind(Zxx, subsims * subsims* subsims)
        if(n==(N-1))
          Zcont<-0
        else
          Zcont<-predict(reg[[n+1]], data.frame(x1=Zxx[,1], x2=Zxx[,2], x3=Zxx[,3], x4=Zxx[,4]))
        ZVhut<-pmax(Zpay,Zcont)

        delta<-Vhut-mean(ZVhut)
        mar[n+1] <- mar[n] + delta
      }

     martingales[[i]]<-mar
   }

    ergs<-numeric(M)
    for(i in 1:M)
    {
      pay<-numeric(N)
      for (k in 1:N) {
       pay[k]<-payoffs[[k]][i]
      }

      ergs[i]<-max(pay-martingales[[i]])
    }

   results<-list()
   results[[1]]<-LSresults[[1]]
   results[[2]]<-mean(ergs)
   return (results)
 }


