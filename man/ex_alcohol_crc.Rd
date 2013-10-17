\name{ex_alcohol_crc}
\alias{ex_alcohol_crc}
\docType{data}
\title{
Several incidence-rate studies
}
\description{
Summarized data about the relation between total alcohol intake
and colorectal cancer risk in 8 prospective cohort
studies participating in the Pooling Project of Prospective
Studies of Diet and Cancer.
}
\usage{data(ex_alcohol_crc)}
\format{
  A data frame with 48 observations on the following 7 variables.
  \describe{
    \item{\code{study}}{a factor with levels \code{atm} \code{hpm} \code{io} \code{n1} \code{n2} \code{nt} \code{ntm} \code{nym}}
    \item{\code{type}}{a factor with levels \code{ir}}
    \item{\code{dose}}{a numeric vector}
    \item{\code{cases}}{a numeric vector}
    \item{\code{peryears}}{a numeric vector}
    \item{\code{logrr}}{a numeric vector}
    \item{\code{se}}{a numeric vector}
  }
}
\details{
There are no more details required}
\source{
Orsini N, Ruifeng L, Wolk A, Khudyakov P, Spiegelman D.
Meta-analysis for linear and non-linear dose-response relationships: examples, an evaluation of approximations, and software.
American Journal of Epidemiology. 2012; 175(1):66-73. 
}
\references{
There are no references
}
\examples{
data(ex_alcohol_crc)
## maybe str(ex_alcohol_crc) ; plot(ex_alcohol_crc) ...
}
\keyword{datasets}
