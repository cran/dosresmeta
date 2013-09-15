\name{hamling}
\alias{hamling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hamling's method
}
\description{
This function approximates the covariance matrix of log relative risks according to the method described by Hamling (2008).
}
\usage{
hamling(cases, n, logrr, v, study, id, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cases}{
	The variable (or vector) reporting the number of cases, regardless of the type of study.
}
  \item{n}{
	The variable(or vector) reporting the number of the size (for each exposure) of the study for case-control and cumulative
	incidence data, the Person-Times (for each exposure) for incidence rate data.
}
  \item{logrr}{
	The variable (or vector) reporting the log relative risks for each exposure
}
  \item{v}{
	The variable (or vector) reporting the variance of published log relative risks. If not already in the dataset, it needs to be calculated by the user.
}
  \item{study}{
	the variable (or vector) or value reporting the type of study. The values for the variable must be "\code{cc}", "\code{ir}", "\code{ci}" for
case-control, incidence-rate and cumulative-incidence data respectively.
}
  \item{id}{
	an optional variable (or vector) reporting the id of the studies. 
}
  \item{data}{
	an optional dataframe containing the variables needed to approximate the covariance matrix.
}
}
\details{
No further details required
}
\value{
A matrix reporting the pseudo-counts originated from the pubblished data.}
\references{
Hamling, J., Lee, P., Weitkunat, R.,  Ambuhl, M. (2008). Facilitating meta-analyses by deriving relative effect and 
precision estimates for alternative comparisons from a set of estimates presented by exposure level or disease category. 
Statistics in medicine, 27(7), 954-970.}
\author{
Alessio Crippa
}
\examples{
data(cc_ex)
cc_ex$v<-((cc_ex$logub-cc_ex$loglb)/(2*qnorm(.975)))^2
hamling(case,n,logrr,v,study="cc",data=cc_ex)
}