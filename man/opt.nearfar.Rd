\name{opt.nearfar}
\alias{opt.nearfar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds optimal near-far match
}
\description{
Discovers optimal near-far matches using the partial F statistic (for
continuous treatments) or residual deviance (for binary treatments);
PLEASE NOTE required variable ordering
}
\usage{
opt.nearfar(dta, trt.bin = FALSE, imp.var = NA, tol.var = NA,
 adjust.IV = TRUE, sink.range = c(0, 0.5), cutp.range = NA, max.time.seconds = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame wherein first column is outcome, second column is treatment, third column is IV, and fourth through last columns are measured confounders}
  \item{trt.bin}{TRUE if treatment of interest is binary; FALSE otherwise}
  \item{imp.var}{A list of (up to 5) named variables to prioritize in the ``near'' matching}
  \item{tol.var}{A list of (up to 5) tolerances attached to the prioritized variables where 0 is highest priority}
  \item{adjust.IV}{if TRUE, include measured confounders in treatment~IV model that is optimized; if FALSE, exclude}
  \item{sink.range}{A two element vector of (min, max) for range of
    sinks over which to optimize in the near-far match; default (0, 0.5) such that maximally 50\% of observations can be removed}
  \item{cutp.range}{a two element vector of (min, max) for range of cutpoints (how far apart the IV will become) over which to optimize in the near-far match; default is (one SD of IV, range of IV)}
  \item{max.time.seconds}{How long to let the optimization algorithm run; default is 300 seconds = 5 minutes}
}
\details{PLEASE NOTE required variable ordering for input data frame
  such that first column is outcome, second column is treatment, third
  column is IV, and fourth through last columns are measured confounders
- otherwise your results will not make sense!  Additionally, if any
absolute standardadized differences are greater than 0.2 in the measured
confounders, you may want to further restrict the sink and cutpoint
ranges over which to optimize.}
\value{
  \item{n.calls}{Number of calls made to the objective function}
  \item{sink.range}{A two element vector of (min, max) for range of
    sinks over which to optimize in the near-far match; default (0, 0.5) such that maximally 50\% of observations can be removed}
  \item{cutp.range}{a two element vector of (min, max) for range of cutpoints (how far apart the IV will become) over which to optimize in the near-far match; default is (one SD of IV, range of IV)}
  \item{pct.sink}{Optimal percent sinks}
  \item{cutp}{Optimal cutpoint}
  \item{maxF}{Highest value of partial F-statistic (continuous treatment)
    or residual deviance (binary treatment) found by simulated annealing optimizer}
  \item{match}{A two column matrix where the first column is the index of an ``encouraged'' individual and the second column is the index of the corresponding ``discouraged'' individual from the pair matching}
  \item{summ}{A table of mean variable values for both the ``encouraged'' and
``discouraged'' groups across all variables plus absolute standardized
differences for each variable}

}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\references{
Lu B, Greevy R, Xu X, Beck C (2011). Optimal nonbipartite matching and its statistical
applications. The American Statistician, 65(1), 21-30.

Xiang Y, Gubian S, Suomela B, Hoeng J (2013). Generalized Simulated Annealing for
Efficient Global Optimization: the GenSA Package for R. The R Journal, 5(1). URL
http://journal.r-project.org/.
}

\examples{
#Generate data
set.seed(81)
dta = mvrnorm(100,c(10,10,10),matrix(c(1,-0.5,0.5,-0.5,1,0.5,0.5,0.5,1),3,3))
Zstar = dta[,1] #Part of Z that is correlated with unmeas conf
X.unmeas = dta[,2] #Unmeas conf
X.meas = dta[,3] #Meas conf
IV = rnorm(100,10,1) #Instrumental variable
Z = 1+5*Zstar+3*X.meas+1*IV+rnorm(100,0,10) #Observed treatment
Y = 1+1*Z+1*X.meas+5*X.unmeas+rnorm(100,0,20) #Outcome
df.sim = data.frame(Y=Y,Z=Z,IV=IV,X=X.meas) #set up for near-far match 
#(X is measured confounder)

#Execute near-far match (just for illustration)
#The default setting is max.time.seconds=300
nf = opt.nearfar(dta=df.sim,trt.bin=FALSE,imp.var=NA,tol.var=NA,
    adjust.IV=TRUE,max.time.seconds=3)

#Look at absolute standardized difference of variables after match
nf$summ

#Illustrate inference using effect ratio
eff.ratio(dta=df.sim,match=nf$match,alpha=0.05)

#Illustrate inference using 2SLS post-near-far match
df.sim3 = df.sim[as.numeric(nf$match),]
m.post = ivreg(Y~X+Z|IV+X,data=df.sim3)
summary(m.post)$coeff
}
