% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\docType{data}
\name{cc_ex}
\alias{cc_ex}
\title{Case-control data on alcohol and breast cancer risk}
\format{
A data frame with 4 observations on the following 10 variables:
\tabular{ll}{
\code{gday} \tab label for exposure levels.\cr
\code{dose} \tab assigned dose levels.\cr
\code{case} \tab number of cases for each exposure level.\cr
\code{control} \tab number of controls for each exposure level.\cr
\code{n} \tab total number of subjects for each exposure level.\cr
\code{crudeor} \tab unadjusted odds ratios for each exposure level.\cr
\code{adjrr} \tab adjusted odds ratios for each exposure level.\cr
\code{lb} \tab lower bound for the confidence limits of the adjusted odds ratios.\cr
\code{ub} \tab upper bound for the confidence limits of the adjusted odds ratios.\cr
\code{logrr} \tab natural logarithm of the adjusted odds ratios.\cr
}
}
\description{
The dataset reports the summarized dose-response results from a case-control study
on alcohol and breast cancer, first presented by Rohan and McMichael.
}
\references{
Rohan, T. E., McMichael, A. J. (1988). Alcohol consumption and risk op breast cancer. International journal of cancer, 41(5), 695-699.

Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis.
American journal of epidemiology, 135(11), 1301-1309.
}
\author{
Alessio Crippa, <\email{alessio.crippa@ki.se}>
}
\keyword{data}
