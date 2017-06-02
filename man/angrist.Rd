\name{angrist}
\alias{angrist}
\docType{data}
\title{Angrist data set for education and wages}
\description{A random sample of 1000 observations from the data set
  used by Angrist and Krueger in their investigation of the impact '
  of education on future wages.}
\format{
  A data frame with 1000 observations on the following 7 variables.
  \describe{
    \item{\code{wage}}{a numeric vector}
    \item{\code{educ}}{a numeric vector}
    \item{\code{qob}}{a numeric vector}
    \item{\code{IV}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{married}}{a numeric vector}
    \item{\code{race}}{a numeric vector}
  }
}
\details{
This data set is a random sample of 1000 observations from the URL
listed below.
}
\source{
\url{http://economics.mit.edu/faculty/angrist/data1/data/angkru1991}
}
\references{
  Angrist JD, Krueger AB (1991). Does Compulsory School Attendance
  Affect Schooling and Earnings? The Quarterly Journal of Economics,
  106(4), 979-1014.
}
\examples{
library(nearfar)
str(angrist)
## maybe str(angrist) ; plot(angrist) ...
}
\keyword{datasets}
