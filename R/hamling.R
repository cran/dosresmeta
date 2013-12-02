hamling <-
function (cases , n, logrr, v, study, id, data) {
## If data argument is provided, it evaluates the other arguments in the actual data
## n is defined as total subject for cc and ci, person-time for ir
if( !missing(data) ) {
   data <- data.frame(data)
mf <- match.call()
mf.cases <- mf[[match("cases", names(mf))]]
mf.n <- mf[[match("n", names(mf))]]
mf.logrr <- mf[[match("logrr", names(mf))]]
mf.v <- mf[[match("v", names(mf))]]
mf.study <- mf[[match("study", names(mf))]]
mf.id <- mf[[match("id", names(mf))]]
cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
logrr <- eval(mf.logrr, data, enclos = sys.frame(sys.parent()))
v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
study <- eval(mf.study, data, enclos = sys.frame(sys.parent()))
id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
if( is.null(id) ) id <- 1
}
## For a single study id could not be given, so it is set equal to 1
if( missing(id) ) id <- 1

## Function to be optimized
## For more details see Hamling 2007
fun <- function(cc, ris=F){
A <- cc[1]
B <- cc[2]

      if (study[id==j][1] == "cc" | study[id==j][1] == 1) {
        cases1 <- (1 + (A/B) * exp(logrr[id==j & v!=0])) / (v[id==j & v!=0] - 1/A -1/B)
          n1 <- (1 + B / (A * exp(logrr[id==j & v!=0]))) / (v[id==j & v!=0] - 1/A - 1/B)
        }
        if ((study[id==j][1] == "ir" | study[id==j][1] == 2) |
(study[id==j][1] == "ci" | study[id==j][1] == 3)){ 
          cases1 <- 1 / (v[id==j & v!=0] - 1/A)
          n1 <- (B / ( A * exp(logrr[id==j & v!=0]) ) ) / (v[id==j & v!=0] - 1/A )
        }
      if (ris==F) {
        p1 <- B / (B + sum(n1))
        z1<- (B + sum(n1)) / (A + sum(cases1))
(((p - p1) / p)^2 + ((z - z1) / z)^2)
} else {
        return( data.frame(cases1 = c(A, cases1), n1 = c(B, n1)) )
}
}

## Obtained pscounts for several studies
pscounts <- data.frame()
for (j in unique(id)){
if (study[id==j][1] == "cc" | study[id==j][1] == 1){
      controls <- n[id==j] - cases[id==j]
      p <- controls[v[id==j]==0]/sum(controls)
      z <- sum(controls)/sum(cases[id==j])
      st.in <- c(cases[id==j & v==0], controls[v[id==j]==0])
   }
      if ( (study[id==j][1] == "ir" | study[id==j][1] == "ci") |
(study[id==j][1] == 2 | study[id==j][1] == 3) ){
       p <- n[id==j & v==0]/sum(n[id==j])
      z <- sum(n[id==j])/sum(cases[id==j])
      st.in <- c(cases[id==j & v==0], n[id==j & v==0])
}
## Obtaining "optimal" values for A and B
st.fin <- optim(st.in, fun)$par
## Obtaining cases1 & n1 from A & B
pscounts <- rbind(pscounts, fun(st.fin, ris=T))
    }
pscounts
}
