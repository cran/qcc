\name{process.capability.sixpack}
\alias{process.capability.sixpack}
\title{Process capability sixpack plots}
\description{Draw a process capability sixpack plots.}
}
\usage{
process.capability.sixpack(object, spec.limits, target, nsigmas, std.dev)
}
\arguments{
\item{object}{a `qcc' object of type \code{"xbar"}}
\item{spec.limits}{a vector specifying the lower and upper specification limits.}
\item{target}{a value specifying the target of the process. If missing the value from the `qcc' object is used if not \code{NULL}, otherwise the target is set at the middle value bewteen specification limits.}
\item{nsigmas}{a numeric value specifying th number of sigmas to use. If not provided is taken from the `qcc' object.}
\item{std.dev}{a value specifying the within-group standard deviation. If not provided is taken from the `qcc' object.}
}
\details{The process capability sixpack plots are: \cr
a xbar chart \cr
a R or S chart (if sample sizes > 10) \cr
a run chart \cr
a histogram \cr
a normal Q-Q plot \cr
a capability plot 
}
%\value{}
\references{
Minitab (2000) \emph{Users's Guide 2: Data Analysis and Quality Tools}, rel. 13, Chapter 14. 
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
\seealso{\code{\link{process.capability}}, \code{\link{qcc}}}
\examples{
x <- matrix(rnorm(100), ncol=5)
q <- qcc(x, type="xbar", plot=FALSE)
process.capability.sixpack(q, spec.limits = c(-2,2), target=0)

data(pistonrings)
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)
q <- qcc(diameter[1:25,], type="xbar", nsigmas=3, plot=FALSE)
process.capability.sixpack(q, spec.limits=c(73.95,74.05))
}
\keyword{htest}
\keyword{hplot}

