\name{stats.R}
\alias{stats.R}
\alias{sd.R}
\alias{limits.R}

\title{Functions to plot Shewhart R chart}
\description{These functions are used to compute statistics required by the R chart.}
\usage{
stats.R(data, sizes)
sd.R(data, sizes, std.dev, equal.sd = TRUE)
limits.R(center, std.dev, sizes, conf)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes. Optional}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.R} function, required for \code{limits.R}}
  \item{equal.sd}{a logical indicating whether the standard deviation is constant for all samples/groups.}
  \item{conf}{a numeric value used to compute control limits, specifying the number of standard deviations (if \code{conf} > 1) or the confidence level (if 0 < \code{conf} < 1).}
}
%\details{}
\value{
The function \code{stats.R} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.R} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.R} returns a matrix with lower and upper control limits. 
}
\references{
Montgomery, D.C. (1991)  \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process
Control}, New York, Chapman and Hall.
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{}
\seealso{\code{\link{qcc}} }
%\examples{}
\keyword{htest}
\keyword{hplot}

