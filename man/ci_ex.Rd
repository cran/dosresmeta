\name{ci_ex}
\alias{ci_ex}
\docType{data}
\title{
Cumulative incidence data
}
\description{
Cumulative incidence data on high-fat dairy food and colorectal cancer risk (Larsson,Bergkvist, and Wolk 2005)
}
\usage{data(ci_ex)}
\format{
  A data frame with 5 observations on the following 9 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{cases}}{a numeric vector}
    \item{\code{n}}{a numeric vector}
    \item{\code{adjrr}}{a numeric vector}
    \item{\code{lb}}{a numeric vector}
    \item{\code{ub}}{a numeric vector}
    \item{\code{logrr}}{a numeric vector}
    \item{\code{loglb}}{a numeric vector}
    \item{\code{logub}}{a numeric vector}
  }
}
\details{
There are no more details required}
\source{
Orsini, N., Bellocco, R., Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
}
\references{
There are no references
}
\examples{
data(ci_ex)
## maybe str(ci_ex) ; plot(ci_ex) ...
}
\keyword{datasets}
