\name{cusum.qcc}
\alias{cusum}
\alias{cusum.default}
\alias{cusum.qcc}
\title{Cusum chart}
\description{Draw a Cusum chart from an object of class `qcc'.}
\usage{
\method{cusum}{qcc}(object, decision.int = 5, se.shift = 1,
      label.bounds = c("LDB", "UDB"), add.stats = TRUE, chart.all = TRUE, 
      ylim = NULL, axes.las = 0, restore.par = TRUE, \dots)
}
\arguments{
\item{object}{an object of class `qcc'.}

\item{decision.int}{A numeric value specifying the number of standard errors of the summary statistics at which the cumulative sum is out of control.}

\item{se.shift}{The amount of shift to detect in the process, measured in standard errors of the summary statistics.}

\item{label.bounds}{a character vector specifying the labels for the the decision interval boundaries.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}
\item{\dots}{}
}
\details{Cusum charts display how the group summary statistics deviate above or below the process center or target value, relative to the standard errors of the summary statistics. Useful to detect small and permanent variation on the mean of the process.
}
\value{Returns an object of class `cusum' which inherits from the `qcc' object. No methods are specifically defined for the `cusum' class.}
\references{
Montgomery, D.C. (2000) \emph{Introduction to Statistical Quality Control}, 4th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}, \code{\link{ewma}}}
\examples{
data(pistonrings)
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)
q <- qcc(diameter[1:25,], type="xbar", nsigmas=3, plot=FALSE)
cusum(q)
q <- qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], nsigmas=3, plot=FALSE)
cusum(q)
cusum(q, chart.all=FALSE)
cusum(qcc(diameter, type="xbar", nsigmas=3, target=74, std.dev=0.02, plot=FALSE))
}
\keyword{htest}
\keyword{hplot}
