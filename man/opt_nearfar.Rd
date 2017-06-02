\name{opt_nearfar}
\alias{opt_nearfar}
\title{Finds optimal near-far match}
\description{
Discovers optimal near-far matches using the partial F statistic (for
continuous treatments) or partial deviance (for binary and treatments)
}
\usage{
opt_nearfar(dta, trt, covs, iv, trt.type = "cont", imp.var = NA,
tol.var = NA, adjust.IV = TRUE, sink.range = c(0, 0.5), cutp.range = NA,
max.time.seconds = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dta}{The name of the data frame on which matching was performed}
  \item{trt}{The name of the treatment variable, e.g., ``educ''}
  \item{iv}{The name of the instrumental variable, e.g., iv="QOB"}
  \item{covs}{A vector of the names of the covariates to make ``near'',
    e.g., covs=c("age", "sex", "race")}
  \item{trt.type}{Treatment variable type: ``cont'' for continuous, or
    ``bin'' for binary}
  \item{imp.var}{A list of (up to 5) named variables to prioritize in
    the ``near'' matching}
  \item{tol.var}{A list of (up to 5) tolerances attached to the
    prioritized variables where 0 is highest penalty for mismatch}
  \item{adjust.IV}{if TRUE, include measured confounders in treatment~IV
    model that is optimized; if FALSE, exclude}
  \item{sink.range}{A two element vector of (min, max) for range of
    sinks over which to optimize in the near-far match; default (0, 0.5)
    such that maximally 50\% of observations can be removed}
  \item{cutp.range}{a two element vector of (min, max) for range of
    cutpoints (how far apart the IV will become) over which to optimize
    in the near-far match; default is (one SD of IV, range of IV)}
  \item{max.time.seconds}{How long to let the optimization algorithm
    run; default is 300 seconds = 5 minutes}
}

\value{
  \item{n.calls}{Number of calls made to the objective function}
  \item{sink.range}{A two element vector of (min, max) for range of
    sinks over which to optimize in the near-far match; default
    (0, 0.5) such that maximally 50\% of observations can be removed}
  \item{cutp.range}{a two element vector of (min, max) for range of
    cutpoints (how far apart the IV will become) over which to optimize
    in the near-far match; default is (one SD of IV, range of IV)}
  \item{pct.sink}{Optimal percent sinks}
  \item{cutp}{Optimal cutpoint}
  \item{maxF}{Highest value of partial F-statistic (continuous treatment)
    or residual deviance (binary treatment) found by simulated annealing
    optimizer}
  \item{match}{A two column matrix where the first column is the index
    of an ``encouraged'' individual and the second column is the index
    of the corresponding ``discouraged'' individual from the pair matching}
  \item{summ}{A table of mean variable values for both the
    ``encouraged'' and ``discouraged'' groups across all variables
    plus absolute standardized differences for each variable}

}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\references{
Lu B, Greevy R, Xu X, Beck C (2011). Optimal nonbipartite matching and
its statistical applications. The American Statistician, 65(1), 21-30.

Xiang Y, Gubian S, Suomela B, Hoeng J (2013). Generalized Simulated
Annealing for Efficient Global Optimization: the GenSA Package for R.
The R Journal, 5(1). URL http://journal.r-project.org/.
}

\examples{
k = opt_nearfar(dta=mtcars, trt="drat", covs=c("cyl", "disp"),
    trt.type="cont", iv="carb", imp.var=NA, tol.var=NA, adjust.IV=TRUE,
    max.time.seconds=2)
summary(k)
}