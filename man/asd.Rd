\name{asd}
\alias{asd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute absolute standardized difference of a variable
between two groups
}
\description{
Called by the \code{\link{summ.matches}} function to compute absolute
standardized differences between pair matches for each variable in a data frame.
}
\usage{
asd(x, match)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A continuous or binary variable}
  \item{match}{A two-column matrix of row indices of paired matches}
}

\value{
Function returns the absolute standardized difference
}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{summ.matches}}, \code{\link{opt.nearfar}}}

