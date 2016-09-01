\name{eff.ratio}
\alias{eff.ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference for effect ratio
}
\description{Conducts inference on effect ratio as described in Section
  3.3 of Baiocchi (2010), resulting in an estimate and a permutation
  based confidence interval for the effect ratio.
}
\usage{
eff.ratio(dta, match, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame where first column is outcome and second column
    is treatment}
  \item{match}{Data frame where first column contains indices for those individuals
    encouraged into treatment by instrumental variable and second column
    contains indices for those individuals discouraged from treatment by
    instrumental variable}
   \item{alpha}{Level of confidence interval}
}

\value{
  \item{est.emp}{Empirical estimate of effect ratio}
  \item{est.HL}{Hodges-Lehmann type estimate of effect ratio}
  \item{lower}{Lower limit to 1-alpha/2 confidence interval for effect ratio}
  \item{upper}{Upper limit to 1-alpha/2 confidence interval for effect ratio}
}
\references{
  Baiocchi M, Small D, Lorch S, Rosenbaum P (2010). Building a stronger instrument in
an observational study of perinatal care for premature infants. Journal of the American
Statistical Association, 105(492), 1285-1296.
}
\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{eff.ratio}}}