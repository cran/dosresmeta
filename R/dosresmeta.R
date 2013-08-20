dosresmeta <-
function(formula,study,se,cov,data,fcov,
               method="reml",psmethod="h",alpha=.05){
  
  #setting working data.frame
  if (length(study)==1){
    id<-1
    if (is.factor(study)) study<-tolower(study)
  }   
  if (length(study)==2){
    id<-eval(as.name(study[1]),data)
    study<-eval(as.name(study[2]),data)
    if (is.factor(study)) study<-tolower(study)
  }    
  n<-eval(as.name(cov[1]),data)
  cases<-eval(as.name(cov[2]),data)
  ds<-data.frame(id=as.double(id),study,cases,n)
  mf<-model.frame(formula, data=data)
  ds$logrr=model.response(mf, "numeric")
  ds=cbind(ds,model.matrix(attr(mf, "terms"), data=mf))
  if (length(se)==1){
    ds$v=eval(as.name(se),data)^2
  }
  if (length(se)==2){
    loglb=eval(as.name(se[1]),data)
    logub=eval(as.name(se[2]),data)
    v=((logub-loglb)/(2*qnorm(alpha/2)))^2
    ds=cbind(ds,v)
  }
  #for reference category se() must not be missing (set to 0)
  ds$v[is.na(ds$v)]=0

  #results  
  param<-data.frame()
  fit.stat<-data.frame()
  Ccov<-list()
  covparam<-list()
  tmf<-data.frame()

  formula=as.formula(paste('logrr~0+',as.name(as.character(formula[3]))))
  if (psmethod=="h") ds<-cbind(ds,hamling(data=ds)[,-1])
  if (psmethod=="gl") ds<-cbind(ds,gl(data=ds)[,-1])
  for (num in as.double(names(table(ds$id)))){
    dsk<-ds[ds$id==num,]
    #approximating correlation matrix between logRR
    if (psmethod!="f"){
    matr<-dsk[dsk$v!=0,]
    matr1<-dsk[dsk$v==0,]
    rcorr<-matrix(,,ncol=dim(matr)[1],nrow=dim(matr)[1])
    if (ds$study[1]=="cc"|ds$study[1]==1){ 
      for (i in 1:nrow(matr)){
        for (j in 1:nrow(matr)){
          rcorr[i,j]<-ifelse(i!=j,(1/matr1$cases1+1/matr1$controls1)/(
            ((1/matr1$cases1+1/matr1$controls1+1/matr$cases1[i]+1/matr$controls1[i])^.5)*
              ((1/matr1$cases1+1/matr1$controls1+1/matr$cases1[j]+1/matr$controls1[j])^.5)),1)
        }}
    }
    if ((ds$study[1]=="ir"|ds$study[1]=="ci")|(ds$study[1]==2|ds$study[1]==3)){ 
      for (i in 1:nrow(matr)){
        for (j in 1:nrow(matr)){
          rcorr[i,j]<-ifelse(i!=j,(1/matr1$cases1-1/matr1$controls1)/(
            ((1/matr1$cases1-1/matr1$controls1+1/matr$cases1[i]-1/matr$controls1[i])^.5)*
              ((1/matr1$cases1-1/matr1$controls1+1/matr$cases1[j]-1/matr$controls1[j])^.5)),1)
        }}
    }
    #approximating covariance matrix between logRR
    ccov<-rcorr
    for (i in 1:nrow(matr)){
      for (j in 1:nrow(matr)){
        ccov[i,j]<-rcorr[i,j]*((matr$v[i]*matr$v[j])^.5)
      }}
    }
    if(psmethod=="f") ccov<-fcov[[num]]

    mf<-model.frame(formula, data=matr)
    x<-t(model.matrix(attr(mf, "terms"), data=mf))

    #Parameter estimates
    vb<-solve(x%*%solve(ccov)%*%t(x))
    beta<-vb%*%(x%*%solve(ccov)%*%t(t(matr$logrr)))
    betascore<-c()
    betapv<-c()
    for (i in 1:length(beta)){
     betascore[i]<-beta[i]/(vb[i,i]^.5)
     betapv[i]<-2*(1-pnorm(abs(betascore[i])))
    }
    
    #Transforming variables
    tmf<-rbind(tmf,cbind(id=num,as.data.frame(solve(t(chol(ccov)))%*%as.matrix(mf)))
    );colnames(tmf)<-c("id",colnames(mf))
    
    
    #GOF parameters
    Q<-(matr$logrr-t(beta)%*%x)%*%solve(ccov)%*%t(matr$logrr-t(beta)%*%x)
    pq<-1-pchisq(Q,dim(matr)[1]-dim(beta)[1])
    ll<--.5*(dim(matr)[1])*log(2*pi)-5*log(det(ccov))-.5*Q
    
    value<-c()
    coeff<-c()
    se<-c()
    for(l in 1:length(beta)){
      coeff[l]<-paste("beta",l,sep="")
      value[l]<-beta[l,1]
      se[l]<-vb[l,l]^.5
    }
    param<-rbind(param,data.frame(id=dsk$id[1],coeff,value,se,score=betascore,pvalue=betapv))
    fit.stat<-rbind(fit.stat,c(dsk$id[1],Q=Q,pq=pq,logll=ll));names(fit.stat)<-c("id","Q","pvalue","Logll") 
    #covariance matrix of d-r estimates
    covparam<-c(covparam,list(vb))
    #estimated covariance matrix among logRR
    Ccov<-c(Ccov,list(ccov))
  
  }
  
  if (psmethod!="f") psnumber<-data.frame(id=ds$id,cases1=ds$cases1,controls1=ds$controls1)
  
  procedure<-list()
  if(max(psnumber$id)!=1){
    coeff<-data.frame()
    for (num in as.double(names(table(param$id)))){
      coeff<-rbind(coeff,param[param$id==num,]$value)
    }
    colnames(coeff)<-names(table(param$coeff))
    procedure<-mvmeta(as.matrix(coeff),covparam,method=method) 
  }
  
  #R2 tool
  mod<-lm(formula,data=tmf)
  colnames(tmf)[-1]<-paste("t",colnames(tmf)[-1],sep="")
  procedure$R2<-summary(mod)$r.squared
  procedure$R2adj<-summary(mod)$adj.r.squared
  
  procedure$psnumber<-psnumber
  procedure$param<-param 
  procedure$covparam<-covparam
  procedure$Ccov<-Ccov
  procedure$fit.stat<-fit.stat
  procedure$tdata<-tmf

  if(max(psnumber$id)==1){
    procedure=list(Param=procedure$param,Sigma=as.matrix(procedure$Ccov[[1]]),
         Fit.Stat=procedure$fit.stat)
  } 
  procedure
}
