\name{summ_matches}
\alias{summ_matches}
\title{
Computes table of absolute standardized differences
}
\description{
Computes absolute standardized differences for both
continuous and binary variables.  Called by \code{\link{opt_nearfar}} to
summarize results of near-far match.
}
\usage{
summ_matches(dta, iv, covs, match)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dta}{The name of the data frame on which matching was performed}
  \item{iv}{The name of the instrumental variable, e.g., iv="QOB"}
  \item{covs}{A vector of the names of the covariates to make ``near'',
    e.g., covs=c("age", "sex", "race")}
  \item{match}{A two-column matrix of row indices of paired matches}
}
\value{
A table of mean variable values for both the ``encouraged'' and
``discouraged'' groups across all variables plus absolute standardized
differences for each variable
}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\seealso{\code{\link{opt_nearfar}}}

\examples{
k2 = matches(dta=mtcars, covs=c("cyl", "disp"), sinks=0.2, iv="carb",
     cutpoint=2, imp.var=c("cyl"), tol.var=0.03)
summ_matches(dta=mtcars, iv="carb", covs=c("cyl", "disp"), match=k2)

}