\name{gl}
\alias{gl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Greenland and Longnecker method
}
\description{
Approximating covariance matrix according to the method described by Greenland and Longnecker (1992)
}
\usage{
gl(data, epsilon = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the dataset with the information needed to perform the algorithm.
}
  \item{epsilon}{
define how accurate wuold be the solution for the algorithm}
}
\details{
No further details required.}
\value{
A matrix reporting the pseudo-counts originated from the pubblished data.}

\references{
Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, 
with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.}
\author{
Alessio Crippa}