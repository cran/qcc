\name{stats.xbar}
\alias{stats.xbar}
\alias{sd.xbar}
\alias{limits.xbar}

\title{Statistics used in computing and drawing a Shewhart xbar chart}
\description{These functions are used to compute statistics required by the xbar chart.}
\usage{
stats.xbar(data, sizes)
sd.xbar(data, sizes, std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"))
limits.xbar(center, std.dev, sizes, conf)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic}
  \item{sizes}{samples sizes. Optional}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.xbar} function, required for \code{limits.xbar}. See details.}
  \item{conf}{a numeric value used to compute control limits, specifying the number of standard deviations (if \code{conf} > 1) or the confidence level (if 0 < \code{conf} < 1).}
}
\details{
Methods  available for estimating the process standard deviation:
\tabular{ll}{
Method \tab Description \cr
\code{"UWAVE-R"} \tab UnWeighted AVErage of subgroup estimates \cr
                 \tab based on subgroup Ranges \cr
\code{"UWAVE-SD"} \tab UnWeighted AVErage of subgroup estimates  \cr
                 \tab based on subgroup Standard Deviations \cr
\code{"MVLUE-R"} \tab Minimum Variance Linear Unbiased Estimator  \cr
                 \tab computed as a weighted average of subgroups \cr
                 \tab estimates based on subgroup Ranges \cr
\code{"MVLUE-SD"} \tab Minimum Variance Linear Unbiased Estimator  \cr
                 \tab computed as a weighted average of subgroup  \cr
                 \tab estimates based on subgroup Standard Deviations \cr
\code{"RMSDF"} \tab Root-Mean-Square estimator computed as a weighted average of \cr
                 \tab subgroup estimates based on subgroup Standard Deviations
}

\tabular{llll}{
Method \tab \code{"xbar"} \tab \code{"R"} \tab \code{"S"} \cr
\code{"UWAVE-R"} \tab default \tab default \tab not available \cr
\code{"UWAVE-SD"} \tab  \tab not available \tab default \cr
\code{"MVLUE-R"} \tab \tab \tab not available  \cr
\code{"MVLUE-SD"} \tab \tab not available  \tab \cr
\code{"RMSDF"} \tab \tab not available \tab
}

Detailed definitions of formulae implemented are available in the SAS/QC 9.2 User's Guide.
}
\value{
The function \code{stats.xbar} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.xbar} returns \code{std.dev} the standard deviation of the statistic charted. This is based on results from Burr (1969).

The function \code{limits.xbar} returns a matrix with lower and upper control limits. 
}
\references{
Burr, I.W. (1969) Control charts for measurements with varying sample sizes. \emph{Journal of Quality Technology}, 1(3), 163-167. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{}
\seealso{\code{\link{qcc}} }
%\examples{}
\keyword{htest}
\keyword{hplot}

