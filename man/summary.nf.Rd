\name{summary.nf}
\alias{summary.nf}
\title{Summary method for object of class ``nf''}
\description{Displays key information, e.g., number of matches tried,
  and post-match balance, for \code{\link{opt_nearfar}} function}

\usage{
\method{summary}{nf}(object, \dots)
}

\arguments{
  \item{object}{Object of class ``nf'' returned by \code{\link{opt_nearfar}}}
  \item{\dots}{additional arguments affecting the summary produced}
}
\value{
Returns a summary of results from \code{\link{opt_nearfar}} function
}
\author{Joseph Rigdon \email{jrigdon@wakehealth.edu}}

\seealso{\code{\link{opt_nearfar}}}

\examples{
k = opt_nearfar(dta=mtcars, trt="drat", covs=c("cyl", "disp"),
    trt.type="cont", iv="carb", imp.var=NA, tol.var=NA, adjust.IV=TRUE,
    max.time.seconds=1)
summary(k)
}


