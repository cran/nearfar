\name{calipers}
\alias{calipers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matching priority function
}
\description{
Updates given distance matrix to prioritize specified measured
confounders in a pair match.  Used in consort with
\code{\link{matches}} function to prioritize specific measured
confounders in a near-far match in the \code{\link{opt.nearfar}} function.
}
\usage{
calipers(distmat, variable, tolerance = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distmat}{An object of class distance matrix}
  \item{variable}{Named variable from list of measured confounders}
  \item{tolerance}{Penalty to apply to mismatched observations; values
    near 0 penalize mismatches more}
}

\value{
Returns an updated distance matrix
}

\seealso{\code{\link{matches}}, \code{\link{opt.nearfar}}}