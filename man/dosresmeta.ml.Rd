% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dosresmeta.R
\name{dosresmeta.ml}
\alias{dosresmeta.ml}
\alias{dosresmeta.reml}
\title{ML and REML Estimators for dosresmeta Models}
\usage{
dosresmeta.ml(Xlist, Zlist, ylist, Slist, nalist, q, nall, control, ...)

dosresmeta.reml(Xlist, Zlist, ylist, Slist, nalist, q, nall, control, ...)
}
\arguments{
\item{Xlist}{a m-dimensional list of study-specific design matrices for the fixed-effects part of the model.}

\item{Zlist}{a m-dimensional list of study-specific design matrices for the random-effects part of the model.}

\item{ylist}{a m-dimensional list of study-specific of vectors of estimated outcomes.}

\item{Slist}{a m-dimensional list of within-study (co)variance matrices of estimated outcomes.}

\item{nalist}{a m-dimensional list of k-dimensional study-specific logical vectors, identifying missing outcomes.}

\item{q, nall}{numeric scalars: number of predictors, number of observations (excluding missing).}

\item{control}{list of parameters for controlling the fitting process, usually internally set to default values by \code{dosresmeta.control}.}

\item{...}{further arguments passed to or from other methods. Currently not used.}
}
\value{
These functions return an intermediate list object, whose components are then processed by \code{\link{dosresmeta.fit}}. Other components are added later 
through \code{\link{dosresmeta}} to finalize an object of class "\code{dosresmeta}".
}
\description{
These functions implement maximum likeliihood (ML) and restricted maximum likelihood (REML) estimators for random-effects dose-response 
meta-analysis. They are meant to be used internally and not directly run by the users.
}
\details{
The estimation involves \eqn{p} fixed-effects coefficients and the \eqn{p(p+1)/2} random-effects parameters defining the between-study (co)variance matrix.
The hybrid estimation procedure is based first on few runs of iterative generalized least square algorithm and then quasi-Newton iterations, 
using specific likelihood functions, until convergence. The estimation algorithm adopts a profiled (or concentrated) approach, that is expressed 
only in terms of the random-effects parameters. Cholesky and and QR decompositions are used for computational stability and efficiency, and for assuring the 
positive-definiteness of the estimated between-study (co)variance matrix. See the help page for the likelihood functions for further details.
}
\examples{

data("alcohol_cvd")

## Random-effect dose-response model assuming linearity, ML estimator
lin.ml <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
                     se = se, cases = cases, n = n, data = alcohol_cvd,
                     , method = "ml")
summary(lin.ml)

## Random-effect dose-response model assuming linearity, REML estimator
lin.reml <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
                       se = se, cases = cases, n = n, data = alcohol_cvd,
                       , method = "reml")
summary(lin.reml)

}
\references{
Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. 
Statistics in Medicine, 31(29), 3821-3839.
}
\seealso{
\code{\link{dosresmeta}}, \code{\link{dosresmeta-package}}, \code{\link{dosresmeta.ml}}
}
\author{
Alessio Crippa, \email{alessio.crippa@ki.se}
}
