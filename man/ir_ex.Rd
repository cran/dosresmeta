\name{ir_ex}
\alias{ir_ex}
\docType{data}
\title{
Incidence-rate data
}
\description{
Incidence-rate data on fiber intake and coronary heart disease risk (Wolk et al. 1999)
}
\usage{data(ir_ex)}
\format{
  A data frame with 5 observations on the following 10 variables.
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
    \item{\code{doser}}{a numeric vector}
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
data(ir_ex)
## maybe str(ir_ex) ; plot(ir_ex) ...
}
\keyword{datasets}
