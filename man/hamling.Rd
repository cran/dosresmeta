\name{hamling}
\alias{hamling}
\title{Approximating effective-counts as proposed by Hamling}
\usage{
hamling(logrr, v, cases, n, type, id, data, se, lb, ub, order = TRUE,
  alpha = 0.05)
}
\arguments{
  \item{logrr}{a vector to specify the reported log
  relative risks.}

  \item{v}{a vector to specify the variances of the
  reported log relative risks. Alternatively the user can
  provide the standard error in the \code{se} argument, or
  the confidence interval for the reported relative risks
  in the \code{lb} and \code{ub} arguments.}

  \item{cases}{a vector to specify the number of cases for
  each exposure level.}

  \item{n}{a vector to specify the total number of subject
  for each exposure level. For incidence-rate data \code{n}
  indicates the amount of person-time for each exposure
  level.}

  \item{type}{a vector (or a string) to specify the
  study-specific design. The values for case-control,
  incidence-rate, and cumulative incidence data are
  \code{cc}, \code{ir}, and \code{ci} (or 1, 2, and 3),
  respectively.}

  \item{id}{an optional vector to specify the id variable
  for the studies included in the analysis.}

  \item{data}{an optional data frame (or object coercible
  by \code{\link{as.data.frame}} to a data frame)
  containing the variables in the previous arguments.}

  \item{se}{an optional vector to specify the standard
  error of the reported log relative risks; needed if
  \code{v} is not provided.}

  \item{lb}{an optional vector to specify the lower bound
  of the confidence interval for the reported relative
  risks; needed if \code{v} and \code{se} are not
  provided.}

  \item{ub}{an optional vector to specify the upper bound
  of the confidence interval for the reported relative
  risks; needed if \code{v} and \code{se} are not
  provided.}

  \item{order}{a logical value to specify if the vectors
  need to be sorted. See details.}

  \item{alpha}{a scalar to specify the alpha nominal value
  used in the published data, by defaul equal to .05.}
}
\value{
The results are returned structured in a data frame.
}
\description{
The function \code{hamling} reconstructs the set of
pseudo-numbers (or "effective" numbers) of cases and
non-cases consistent with the input data (log relative
risks) for either a single or multiple summarized data. The
method was proposed in 2008 by Hamling.
}
\details{
The function reconstructs the effective counts
corresponding to the multivariable adjusted log relative
risks as well as their standard errors. A unique solution
is guaranteed by keeping the ratio non-cases to cases and
the fraction of unexposed subjects equal to the unadjusted
data (Hamling). The function requires the data to be sorted
by id and in such a way that the referent values correspond
to the first record for each study. This is automatically
done by the function when \code{order = TRUE} (default).
}
\examples{
data("alcohol_cvd")
hamling(logrr = logrr, se = se, cases = cases, n = n, type = type,
 id = id, data = alcohol_cvd)
}
\author{
Alessio Crippa, \email{alessio.crippa@ki.se}
}
\references{
Hamling, J., Lee, P., Weitkunat, R., Ambuhl, M. (2008).
Facilitating meta-analyses by deriving relative effect and
precision estimates for alternative comparisons from a set
of estimates presented by exposure level or disease
category. Statistics in medicine, 27(7), 954-970.

Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D.
(2012). Meta-analysis for linear and nonlinear
dose-response relations: examples, an evaluation of
approximations, and software. American journal of
epidemiology, 175(1), 66-73.
}

