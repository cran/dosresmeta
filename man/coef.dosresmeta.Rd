\name{coef.dosresmeta}
\alias{coef}
\alias{coef.dosresmeta}
\alias{vcov}
\alias{vcov.dosremeta}
\alias{vcov.dosresmeta}
\title{Extract Coefficients and (Co)Variance Matrix from dosresmeta Objects}
\usage{
\method{coef}{dosresmeta}(object, format = c("vector", "matrix"), ...)

\method{vcov}{dosresmeta}(object, ...)
}
\arguments{
  \item{object}{an object of class "\code{dosresmeta}"}

  \item{format}{format of the returned object}

  \item{\dots}{further arguments passed to or from other
  methods.}
}
\value{
For \code{coef}, a vector (default) or matrix with the
estimated (fixed-effects) coefficients.

For \code{vcov}, the (co)variance matrix of the estimated
(fixed-effects) coefficients.
}
\description{
Extract Coefficients and (Co)Variance Matrix from
dosresmeta Objects
}
\examples{
## Load data and run the model
data("alcohol_cvd")
model <- dosresmeta(formula = logrr ~ dose + I(dose^2), type = type, id = id,
                   se = se, cases = cases, n = n, data = alcohol_cvd)

## Fixed-effect coefficients
coef(model)

## Fixed-effect (co)variance matrix
vcov(model)
}
\author{
Alessio Crippa, \email{alessio.crippa@ki.se}
}
\seealso{
\code{\link{dosresmeta}}, \code{\link{coef}},
\code{\link{vcov}}
}

