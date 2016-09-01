\name{stat.lambda}
\alias{stat.lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute test statistic for effect ratio
}
\description{Computes statistic in used for inference on effect ratio in
  equation (2) of Section 3.3 of Baiocchi (2010).  Called by
  \code{\link{eff.ratio}} when finding estimate and permutation based confidence
  interval for effect ratio.
}
\usage{
stat.lambda(lambda, dta, match)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{Hypothesized value of the effect ratio}
  \item{dta}{Data frame where first column is outcome and second column
    is treatment}
  \item{match}{Data frame where first column contains indices for those individuals
    encouraged into treatment by instrumental variable and second column
    contains indices for those individuals discouraged from treatment by
    instrumental variable
}
}
\value{
Returns the value of the statistic described above
}
\references{
  Baiocchi M, Small D, Lorch S, Rosenbaum P (2010). Building a stronger instrument in
an observational study of perinatal care for premature infants. Journal of the American
Statistical Association, 105(492), 1285-1296.
}
\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{eff.ratio}}}