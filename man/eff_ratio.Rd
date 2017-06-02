\name{eff_ratio}
\alias{eff_ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference for effect ratio
}
\description{Conducts inference on effect ratio as described in Section
  3.3 of Baiocchi (2010), resulting in an estimate and a permutation
  based confidence interval for the effect ratio.
}
\usage{
eff_ratio(dta, match, outc, trt, alpha)
}

\arguments{
  \item{dta}{The name of the data frame object}
  \item{match}{Data frame where first column contains indices for those
    individuals encouraged into treatment by instrumental variable and
    second column contains indices for those individuals discouraged
    from treatment by instrumental variable; returned by both
    \code{\link{opt_nearfar}} and \code{\link{matches}}}
  \item{outc}{The name of the outcome variable in quotes, e.g., ``wages''}
  \item{trt}{The name of the treatment variable, e.g., ``educ''}
   \item{alpha}{Level of confidence interval}
}

\value{
  \item{est.emp}{Empirical estimate of effect ratio}
  \item{est.HL}{Hodges-Lehmann type estimate of effect ratio}
  \item{lower}{Lower limit to 1-alpha/2 confidence interval for effect ratio}
  \item{upper}{Upper limit to 1-alpha/2 confidence interval for effect ratio}
}
\references{
  Baiocchi M, Small D, Lorch S, Rosenbaum P (2010). Building a stronger
  instrument in an observational study of perinatal care for premature
  infants. Journal of the American Statistical Association, 105(492), 1285-1296.
}
\author{Joseph Rigdon \email{jrigdon@stanford.edu}}


\examples{
k2 = matches(dta=mtcars, covs=c("cyl", "disp"), sinks=0.2, iv="carb",
    cutpoint=2, imp.var=c("cyl"), tol.var=0.03)

eff_ratio(dta=mtcars, match=k2, outc="wt", trt="gear", alpha=0.05)
}