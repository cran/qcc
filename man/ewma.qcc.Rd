\name{ewma.qcc}
\alias{ewma}
%\alias{ewma.default}
\alias{ewma.qcc}
\title{EWMA chart}
\description{Draw an EWMA chart from an object of class `qcc'.}
\usage{
\method{ewma}{qcc}(object, lambda=0.2, nsigmas=object$nsigmas, 
     add.stats = TRUE, xlab, ylab, ylim = NULL, axes.las = 0, 
     restore.par = TRUE, \dots)
}
\arguments{
\item{object}{an object of class `qcc'.}

\item{lambda}{the smoothing parameter \eqn{0 \le \lambda \le 1}{0 <= lambda <= 1}}

\item{nsigmas}{a numeric value specifying th number of sigmas to use for computing control limits.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{xlab}{a string giving the label for the x-axis.}

\item{ylab}{a string giving the label for the y-axis.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}

\item{\dots}{}
}
\details{EWMA chart smooths a series of data based on a moving average with weights which decay exponentially. Useful to detect small and permanent variation on the mean of the process.
}
\value{Returns an object of class `ewma' which inherits from the `qcc' object. No methods are specifically defined for the `ewma' class.}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}, New York, Chapman and Hall.
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}, \code{\link{ewmaSmooth}}, \code{\link{cusum}}}
\examples{
data(pistonrings)
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)
q <- qcc(diameter[1:25,], type="xbar", nsigmas=3, plot=FALSE)
ewma(q, lambda=0.2)
q <- qcc(diameter[1:25,], newdata=diameter[26:40,], type="xbar", plot=FALSE)
ewma(q, lambda=0.2, nsigmas=2.7)
}
\keyword{htest}
\keyword{hplot}
