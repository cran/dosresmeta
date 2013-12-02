dosresmeta <-
function(formula, study, se, cov, data,
method="reml", psmethod="h", alpha=.05,
   fcov,ucov){

      data <- data.frame(data)
mf <- match.call()
mf.study <- mf[[match("study", names(mf))]]
mf.cov <- mf[[match("cov", names(mf))]]
mf.se <- mf[[match("se", names(mf))]]

formula <- update(formula, . ~ . +0)
mfm <- model.frame(formula, data)
logrr <- as.matrix(model.response(mfm, "numeric"))
x <- model.matrix(attr(mfm, "terms"), data=mfm)

if (length(mf.study)==1){
study <- eval(mf.study, data, enclos = sys.frame(sys.parent()))
id <- 1
} else{
id <- eval(as.name(as.character(mf.study[2])), data, enclos = sys.frame(sys.parent()))
study <- eval(as.name(as.character(mf.study[3])), data, enclos = sys.frame(sys.parent()))
}
id <- as.numeric(id)
if (is.factor(study) | is.character(study)) study <- tolower(study)
cases <- eval(as.name(as.character(mf.cov[2])), data, enclos = sys.frame(sys.parent()))
n <- eval(as.name(as.character(mf.cov[3])), data, enclos = sys.frame(sys.parent()))
if (length(mf.se)==1){
se <- eval(mf.se, data, enclos = sys.frame(sys.parent()))
v <- se^2
} else{
loglb <- eval(as.name(as.character(mf.se[2])),data,enclos=sys.frame(sys.parent()))
logub <- eval(as.name(as.character(mf.se[3])),data,enclos=sys.frame(sys.parent()))
v <-( (logub - loglb)/(2*qnorm(alpha/2)) )^2
}
v[is.na(v)] <- 0

param <- data.frame()
fit.stat <- data.frame()
tmfm <- data.frame()
Ccov <- list()
covparam <- list()

if (psmethod=="h"){
cases1 <- hamling(cases, n, logrr, v, study, id)$cases1
n1 <- hamling(cases ,n ,logrr ,v ,study ,id)$n1
} 
if (psmethod=="gl"){
cases1 <- grl(cases, n, logrr, v, study, id)$cases1
n1 <- grl(cases, n, logrr, v, study, id)$n1
}

for (j in unique(id)){
if (psmethod=="h" | psmethod=="gl"){
if (study[id==j][1] == "cc" | study[id==j][1] ==1 ){
s0 <- 1/cases1[id==j & v==0] + 1/n1[id==j & v==0]
si <- s0 + 1/cases1[id==j & v!=0]+ 1/n1[id==j & v!=0]
}
if (study[id==j][1] == "ir" | study[id==j][1] == 2){
s0 <- 1/cases1[id==j & v==0]
si <- s0 + 1/cases1[id==j & v!=0]
}
if (study[id==j][1] == "ci" | study[id==j][1] == 3){
s0 <- 1/cases1[id==j & v==0] - 1/n1[id==j & v==0]
si <- s0 + 1/cases1[id==j & v!=0] - 1/n1[id==j & v!=0]
}
rcorr <- s0/(si%*%t(si))^.5
diag(rcorr) <- 1
ccov <- (v[id==j & v!=0]%*%t(v[id==j & v!=0]))^.5*rcorr
}
if (psmethod=="fl"){
ccov <- matrix(fcov, ncol=length(v[id==j & v!=0]),
nrow=length(v[id==j & v!=0]))
diag(ccov) <- v[id==j & v!=0]
}
if (psmethod=="user"){
ccov <- ucov[[j]]
}

vb <- solve(t(x[id == j & v != 0,]) %*% solve(ccov) %*% (x[id == j & v != 0,]))
beta<-vb%*%(t(x[id == j & v != 0,]) %*% solve(ccov) %*% (logrr[id==j & v!=0]))
beta.score <- beta/(diag(vb)^.5)
beta.pvalue <- 2*(1 - pnorm(abs(beta.score)))
Q <- t(logrr[id==j & v!=0] - x[id==j & v!=0,] %*% beta
) %*% solve(ccov) %*% (logrr[id==j & v!=0] - x[id==j & v!=0,] %*% beta)
Q.pvalue <- 1 - pchisq(Q, length(logrr[id==j & v!=0]) - dim(beta)[1])
ll<- -.5*(length(logrr[id==j & v!=0]))*log(2*pi) - .5*log(det(ccov)) - .5*Q

param <- rbind(param,data.frame(id=j,coef=colnames(x),value=beta,
se <- diag(vb)^.5,z=beta.score,pvalue=beta.pvalue))
fit.stat <- rbind(fit.stat,data.frame(id=j,Q=Q,pvalue=Q.pvalue,logll=ll))
covparam <- c(covparam,list(vb))
Ccov <- c(Ccov,list(ccov))

#Transforming variables
tmfm <- rbind(tmfm, data.frame(id=j,
solve(t(chol(ccov))) %*% as.matrix(mfm[id==j & v!=0,])))
}

if (length(unique(id))==1){
procedure <- list()        
} else {
coeff <- data.frame()
    for (j in unique(id)){
      coeff <- rbind(coeff, param[param$id==j,]$value)
}
colnames(coeff) <- param$coef[param$id==1]
procedure <- mvmeta(as.matrix(coeff), covparam, method=method)     

#R2 tool
mod <- lm(formula,data=tmfm)
colnames(tmfm)[-1] <- paste("t",colnames(tmfm)[-1],sep="")
procedure$tdata <- tmfm
procedure$R2 <- summary(mod)$r.squared
procedure$R2adj <- summary(mod)$adj.r.squared
}

colnames(param)<-list("id","","Estimate","Std. Error","z value","Pr(>|z|)" )  
colnames(fit.stat)<-list("id","Q","Pr(>chi2)","log ll")

procedure$Param<-param
procedure$fit.stat<-fit.stat
procedure$Ccov<-Ccov
procedure$covparam<-covparam
procedure$psnumber<-data.frame(cases1,n1)

procedure
}
