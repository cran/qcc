\name{stats.xbar.one}
\alias{stats.xbar.one}
\alias{sd.xbar.one}
\alias{limits.xbar.one}

\title{Functions to plot Shewhart xbar chart for one-at-time data}
\description{These functions are used to compute statistics required by the xbar chart for one-at-time data.}
\usage{
stats.xbar.one(data, sizes)
sd.xbar.one(data, sizes, std.dev, k=2)
limits.xbar.one(center, std.dev, sizes, conf)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes. Not needed, size=1 is used.}
  \item{k}{number of successive pairs of observations for computing the standard deviation based on moving ranges of k points.}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.xbar.one} function, required for \code{limits.xbar.one}}
  \item{conf}{a numeric value used to compute control limits, specifying the number of standard deviations (if \code{conf} > 1) or the confidence level (if 0 < \code{conf} < 1).}
}
%\details{}
\value{
The function \code{stats.xbar.one} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.xbar.one} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.xbar.one} returns a matrix with lower and upper control limits. 
}
\references{
Montgomery, D.C. (2000) \emph{Introduction to Statistical Quality Control}, 4th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
# Water content of antifreeze data (Wetherill and Brown, 1991, p. 120)
x <- c(2.23, 2.53, 2.62, 2.63, 2.58, 2.44, 2.49, 2.34, 2.95, 2.54, 2.60, 2.45,
       2.17, 2.58, 2.57, 2.44, 2.38, 2.23, 2.23, 2.54, 2.66, 2.84, 2.81, 2.39,
       2.56, 2.70, 3.00, 2.81, 2.77, 2.89, 2.54, 2.98, 2.35, 2.53)
sigma <- NA
k <- 2:24
for (j in k)
    sigma[j] <- sd.xbar.one(x, k=j)
# plot estimates of sigma for different values of k
plot(k, sigma[k], type="b")
# "as the size inceases further, we would expect sigma-hat to settle down
#  at a value close to the overall sigma-hat" (p. 121)
abline(h=sd(x), col=2, lty=2)

# the Shewhart control chart for one-at-time data
qcc(x, type="xbar.one", data.name="Water content (in ppm) of batches of antifreeze")
}
\keyword{htest}
\keyword{hplot}

