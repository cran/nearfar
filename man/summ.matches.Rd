\name{summ.matches}
\alias{summ.matches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes table of absolute standardized differences
}
\description{
Computes absolute standardized differences for both
continuous and binary variables.  Called by \code{\link{opt.nearfar}} to
summarize results of near-far match.
}
\usage{
summ.matches(dta, match)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame with n rows (observations) and p columns (variables}
  \item{match}{A two-column matrix of row indices of paired matches}
}
\value{
A table of mean variable values for both the ``encouraged'' and
``discouraged'' groups across all variables plus absolute standardized
differences for each variable
}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{opt.nearfar}}}