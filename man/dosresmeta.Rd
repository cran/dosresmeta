\name{dosresmeta}
\alias{dosresmeta}
\title{
Performing dose-response meta-analysis
}
\description{
The function estimates a dose-response relation for each study included in the analysis, taking into account the correlation among observation. It estimates the covariance matrix according to one of 
the methods described by Greeland S. and Longnecker M., Hamling J., and Easton D. and use it to efficiently estimate the dose-response relation, as described by Orsini N (2011). The estimates are then 
combined according to principle of multivariate random effects meta-analysis to provide an overall pooled functional relation.
}
\usage{
dosresmeta(formula, study, se, cov, data, method = "reml", 
	psmethod = "h", alpha = 0.05, fcov, ucov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
	an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the dose-response relation to be estimated for each study included in the analysis.
	The details of model specification are given under "Details".
}
  \item{study}{
	a 2 or 1-dimensional vector containing the id variable and/or and type (of study) variable. 
}
  \item{se}{
	a 2 or 1-dimensional vector containing containing the variable reporting standar errors, or alternatively the upper and lower bounds of log relative risks.
}
  \item{cov}{
	2 dimensional vector containg the variables needed to approximate the coviariance matrix of log relative risks (cases,n): \code{cases} is the number of cases regardless the type of study;  \code{n} is the size of the study for case-control and cumulative
	incidence data, the Person-Times for incidence rate data.
}
  \item{data}{
	the dataset (or object coercible by "\code{as.data.frame}" to a data frame) containing the variables in the model formula and the variables specified by the previous statements.
}
  \item{method}{
	method used to estimate the pooled dose response relation: "fixed" for fixed-effects models, "ml" or "reml" for random-effects models fitted through (restricted) maximum likelihood, and "mm" for random-
	effects models fitted through method of moments. The default method is "reml". See \code{\link{mvmeta}}.
}
  \item{psmethod}{
	method used to approximate the covariance matrix of log relative risks: "gl" for the method described by Greenland S. and Longnecker M., "h" for the method describe by Hamling J., "fl" for the one
	presented by Easton D, or "user" if specified by the user. In the last case a correct list of matrices bu be specified in the statment "ucov". The default method is "h".
}
  \item{alpha}{
	the scalar indicating the alpha nominal values used in the published data, by defaul equal to .05.
}
  \item{fcov}{
	an optional vector that specifies the covariance (scalar) of log relative risks within the same study. It must be of the same length of the number of studies included in the analyses..
}
  \item{ucov}{
	an optional list of matrices which specifies the covariance matrix of log relative risks within the same study. It must be the same length of the number of studies included in the analyses.
}
}
\details{
The function estimates efficientlly the dose-response relation specified in the formula by estimating the covariance matrix according to the method indicated in the psmethod. This method varies
according to the type of the studies included in the analysis. The variable reporting the type of study must be specified in the \code{study}. The values for the variable must be "\code{cc}", "\code{ir}", "\code{ci}" for
case-control, incidence-rate and cumulative-incidence data respectively.
}
\value{
The function returns an object of class "\code{mvmeta}" in case of multiple studies, a list of objects if estimating a dose-response relation for a single study. 
}
\references{
Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.

Orsini, N., Bellocco, R.,  Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.

Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. Statistics in Medicine, 31(29), 3821-3839.

}
\author{
Alessio Crippa <alessio.crippa@ki.se>
}
\note{
The function required the package \code{\link{mvmeta}} to be installed and loaded.
}


\seealso{
	\code{\link{mvmeta}}
}
\examples{
### FIRST EXAMPLE: Single case-control study
### linear trend using Hamling's method

# Inspect data
data(cc_ex)

# Fitting the model
(mod1<-dosresmeta(formula=logrr~dose,study="cc",
	cov=c(case,n),se=c(loglb,logub),data=cc_ex))

# Parameter Estimate
mod1$Param
# Approximated Covariance Matrix
mod1$Ccov


### SECOND EXAMPLE: Multiple studies
### linear trend using Greenland and Longnecker method
### method of moments

# Inspect data
data(ex_alcohol_crc)

# Fitting the model
(mod2<-dosresmeta(formula=logrr~dose,study=c(study,type),
	cov=c(cases,peryears),se=se,data=ex_alcohol_crc,
	psmethod="gl",method="mm"))

# Summarize the results
summary(mod2)
# residual and fitted values
residuals(mod2)
fitted(mod2)
# Pseudo counts and estimated covariance matrices
mod2$psnumber
mod2$Ccov
}
\keyword{package}
