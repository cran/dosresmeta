hamling <-
function(data){
  pscounts<-data.frame()
  for (k in as.double(names(table(data$id)))){
    ds<-data[data$id==k,]
    #defining p and z
    if (ds$study[1]=="cc"|ds$study[1]==1){ 
      ds$controls<-ds$n-ds$cases
      p<-ds$controls[ds$v==0]/sum(ds$controls)
      z<-sum(ds$controls)/sum(ds$cases)
      st.in<-c(ds$cases[ds$v==0],
               ds$controls[ds$v==0])
    }
    if ((ds$study[1]=="ir"|ds$study[1]=="ci")|(ds$study[1]==2|ds$study[1]==3)){ 
      p<-ds$n[ds$v==0]/sum(ds$n)
      z<-sum(ds$n)/sum(ds$cases)
      st.in<-c(ds$cases[ds$v==0],
               ds$n[ds$v==0])
    }
    cases1<-c()
    controls1<-c()
    fun<-function(cc,ris=F){
      A<-cc[1];B<-cc[2]
      for (i in 1:length(ds$id)){
        if (ds$study[1]=="cc"|ds$study[1]==1){ 
          cases1[i]<-abs((1+(A/B)*exp(ds$logrr[i]))/(ds$v[i]-1/A-1/B))
          controls1[i]<-abs((1+B/(A*exp(ds$logrr[i])))/(ds$v[i]-1/A-1/B))
        }
        if ((ds$study[1]=="ir"|ds$study[1]=="ci")|(ds$study[1]==2|ds$study[1]==3)){ 
          cases1[i]<-abs((1-(A/B)*exp(ds$logrr[i]))/(ds$v[i]-1/A+1/B))
          controls1[i]<-abs((B/(A*exp(ds$logrr[i])-1))/(ds$v[i]-1/A+1/B))
        } 
      }
      ifelse(ris==F,{
        p1<-controls1[1]/sum(controls1)
          z1<-sum(controls1)/sum(cases1)
          ((p-p1)/p)^2+((z-z1)/z)^2},
        return(data.frame(id=ds$id,cases1,controls1)))
    }
    st.fin<-optim(st.in,fun)$par
    pscounts<-rbind(pscounts,round(fun(st.fin,ris=T),0))
 }
  pscounts
}
