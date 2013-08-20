gl<-function(data,epsilon=.001){
  pscounts<-data.frame()
  for (k in as.double(names(table(data$id)))){
    ds<-data[data$id==k,]
    if (ds$study[1]=="cc"|ds$study[1]==1){ 
    fun<-function(aa){
      delta<-1
      repeat{
      a0<-sum(ds$cases)-sum(aa)  
      cov<-1/c(a0,aa)+1/(ds$n-c(a0,aa))  
      e<-ds$logrr[ds$v!=0]+log(a0)+log(ds$n[ds$v!=0]-aa)-
        log(aa)-log(ds$n[ds$v==0]-a0)
      H<-diag(cov[-1]+cov[1])
      for (i in 1:(length(ds$id)-1)){
        for (j in 1:(length(ds$id)-1)){
          if (i!=j) H[i,j]<-cov[1]
        }}
      aa1<-aa+solve(H)%*%e
      delta<-sum((aa1-aa)^2)
      if (delta<epsilon) break
      aa<-aa1
      }
      c(a0,aa)
    }
    cases1<-round(fun(ds$cases[ds$v!=0]),0)
    controls1<-ds$n-cases1
    }
    if ((ds$study[1]=="ir"|ds$study[1]=="ci")|(ds$study[1]==2|ds$study[1]==3)){ 
      fun<-function(aa){
        delta<-1
        repeat{
          a0<-sum(ds$cases)-sum(aa)  
          cov<-1/c(a0,aa)  
          e<-ds$logrr[ds$v!=0]+log(a0)+log(ds$n[ds$v!=0])-
            log(aa)-log(ds$n[ds$v==0])
          H<-diag(cov[-1]+cov[1])
          for (i in 1:(length(ds$id)-1)){
            for (j in 1:(length(ds$id)-1)){
              if (i!=j) H[i,j]<-cov[1]
            }}
          aa1<-aa+solve(H)%*%e
          delta<-sum((aa1-aa)^2)
          if (delta<epsilon) break
          aa<-aa1
        }
        c(a0,aa)
      }
      cases1<-round(fun(ds$cases[ds$v!=0]),0)
      controls1<-ds$n-cases1
    }
    pscounts<-rbind(pscounts,data.frame(id=ds$id,cases1,controls1))
  }
  pscounts
}