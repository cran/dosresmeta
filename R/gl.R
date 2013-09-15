gl <-
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

fun<-function(AA){
	epsilon<-1e-5
   	repeat{
   		A0=sum(cases[id==j])-sum(AA)
		if ((study[id==j][1]=="cc"|study[id==j][1]=="ci")|
			(study[id==j][1]==1|study[id==j][1]==3)){
   			cov<-1/c(A0,AA)+1/(n[id==j]-c(A0,AA))  
   			e<-logrr[id==j & v!=0]+log(A0)+log(n[id==j & v!=0]-AA)-
			log(AA)-log(n[id==j & v==0]-A0)
		}
		if (study[1]=="ir"|study[1]==2){
			cov<-1/c(A0,AA)
   			e<-logrr[id==j & v!=0]+log(A0)+log(n[id==j & v!=0])-
   	  			log(AA)-log(n[id==j & v==0])
		}
		H<-diag(cov[-1]+cov[1])
   		for (i in 1:(length(n[id==j])-1)){
   			for (l in 1:(length(n[id==j])-1)){
       				if (i!=l) H[i,l]<-cov[1]
   			}
		}
		AA1<-AA+solve(H)%*%e
		delta<-2*sum((AA1-AA)^2)
   		if (delta<epsilon) break
		AA<-AA1
	}
	cases1<-round(c(A0,AA),0)
	if (study[id==j][1]=="cc"|study[id==j][1]==1) n1<-n[id==j]-cases1
       		else n1<-n[id==j]
	data.frame(cases1,n1)
}
	
pscounts<-data.frame()
for (j in 1:max(id)){
	pscounts<-rbind(pscounts,fun(cases[id==j & v!=0]))
}
as.data.frame(pscounts)
}