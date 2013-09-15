hamling <-
function(cases,n,logrr,v,study,id,data){

if(!missing(data)){
      data<-data.frame(data)
mf<-match.call()
mf.cases<-mf[[match("cases", names(mf))]]
mf.n<-mf[[match("n", names(mf))]]
mf.logrr<-mf[[match("logrr", names(mf))]]
mf.v<-mf[[match("v", names(mf))]]
mf.study<-mf[[match("study", names(mf))]]
mf.id<-mf[[match("id", names(mf))]]
cases<-eval(mf.cases, data, enclos = sys.frame(sys.parent()))
n<-eval(mf.n, data, enclos = sys.frame(sys.parent()))
logrr<-eval(mf.logrr, data, enclos = sys.frame(sys.parent()))
v<-eval(mf.v, data, enclos = sys.frame(sys.parent()))
study<-eval(mf.study, data, enclos = sys.frame(sys.parent()))
id<-eval(mf.id, data, enclos = sys.frame(sys.parent()))
if(is.null(id)) id<-1
}
if(missing(id)) id<-1

fun<-function(cc,ris=F){
A<-cc[1];B<-cc[2]
        if ((study[id==j][1]=="cc"|study[id==j][1]=="ir")|
(study[id==j][1]==1|study[id==j][1]==2)){
          cases1<-abs((1+(A/B)*exp(logrr[id==j]))/(v[id==j]-1/A-1/B))
          n1<-abs((1+B/(A*exp(logrr[id==j])))/(v[id==j]-1/A-1/B))
        }
        if (study[id==j][1]=="ci"|study[id==j][1]==3){ 
          cases1<-abs((1-(A/B)*exp(logrr[id==j]))/(v[id==j]-1/A+1/B))
          n1<-abs((B/(A*exp(logrr[id==j])-1))/(v[id==j]-1/A+1/B))
        }
       ifelse(ris==F,{
        p1<-n1[1]/sum(n1)
        z1<-sum(n1)/sum(cases1)
          ((p-p1)/p)^2+((z-z1)/z)^2},
        return(data.frame(cases1,n1)))
}

pscounts<-data.frame()
for (j in 1:max(id)){
if (study[id==j][1]=="cc"|study[id==j][1]==1){
      controls<-n[id==j]-cases[id==j]
      p<-controls[v[id==j]==0]/sum(controls)
      z<-sum(controls)/sum(cases[id==j])
      st.in<-c(cases[id==j & v==0],controls[v[id==j]==0])
   }
      if ((study[id==j][1]=="ir"|study[id==j][1]=="ci")|
(study[id==j][1]==2|study[id==j][1]==3)){
       p<-n[id==j & v==0]/sum(n[id==j])
      z<-sum(n[id==j])/sum(cases[id==j])
      st.in<-c(cases[id==j & v==0],n[id==j & v==0])
}
st.fin<-optim(st.in,fun)$par
pscounts<-rbind(pscounts,round(fun(st.fin,ris=T),0))
    }
as.data.frame(pscounts)
}
