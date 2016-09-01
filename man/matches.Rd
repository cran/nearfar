\name{matches}
\alias{matches}

\title{
Function to find pair matches using a distance matrix.  Called by
\code{\link{opt.nearfar}} to discover optimal near-far matches.
}
\description{
Given values of percent sinks and cutpoint, this function will find the
corresponding near-far match
}
\usage{
matches(X, imp.var = NA, tol.var = NA, sinks = 0, IV = NA, cutpoint = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix of measured confounders (with column names) on which to make ``near'' in the matching }
  \item{imp.var}{A list of (up to 5) named variables to prioritize in the ``near'' matching}
  \item{tol.var}{A list of (up to 5) tolerances attached to the
    prioritized variables where 0 is highest penalty for mismatch}
  \item{sinks}{Percentage of the data to match to sinks (and thus remove) if desired; default is 0}
  \item{IV}{Vector of instrumental variable values on which to make ``far'' in the matching}
  \item{cutpoint}{Value below which individuals are too similar on IV; increase to make individuals more ``far'' in match}
}
\details{
Default settings yield a "near" match on only observed confounders in X;
add IV, sinks, and cutpoint to get near-far match.
}
\value{
A two-column matrix of row indices of paired matches
}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{opt.nearfar}}}

\references{
Lu B, Greevy R, Xu X, Beck C (2011). Optimal nonbipartite matching and its statistical
applications. The American Statistician, 65(1), 21-30.
}




