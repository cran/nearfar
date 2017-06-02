\name{smahal}
\alias{smahal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute rank-based Mahalanobis distance matrix between each pair}
\description{
This function computes the rank-based Mahalanobis distance matrix
between each pair of observations in the data set.  Called by
\code{\link{matches}} (and ultimately \code{\link{opt_nearfar}})
function to set up a distance matrix used to create pair matches.
}
\usage{
smahal(X)
}

\arguments{
  \item{X}{A matrix of observed confounders with n rows (observations) and p columns (variables)}
}
\value{
Returns the rank-based Mahalanobis distance matrix between every pair of observations
}

\examples{
smahal(mtcars[1:4, 2:3])
}
