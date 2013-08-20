\name{cc_ex}
\alias{cc_ex}
\docType{data}
\title{
Single case-control study
}
\description{
Case-control data on alcohol and breast cancer risk (Rohan and McMichael 1988)
}
\usage{data(cc_ex)}
\format{
  A data frame with 4 observations on the following 12 variables.
  \describe{
    \item{\code{gday}}{a character vector}
    \item{\code{dose}}{a numeric vector}
    \item{\code{case}}{a numeric vector}
    \item{\code{control}}{a numeric vector}
    \item{\code{n}}{a numeric vector}
    \item{\code{crudeor}}{a numeric vector}
    \item{\code{adjrr}}{a numeric vector}
    \item{\code{lb}}{a numeric vector}
    \item{\code{ub}}{a numeric vector}
    \item{\code{logrr}}{a numeric vector}
    \item{\code{loglb}}{a numeric vector}
    \item{\code{logub}}{a numeric vector}
  }
}
\details{
There are no more details required
}
\source{
Orsini, N., Bellocco, R., Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
}
\references{
There are no references
}
\examples{
data(cc_ex)
## maybe str(cc_ex) ; plot(cc_ex) ...
}
\keyword{datasets}
