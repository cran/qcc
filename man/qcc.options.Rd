\name{qcc.options}
\alias{qcc.options}
\alias{.qcc.options}

\title{Set or return options for the `qcc' package.}

\description{This function can be used to control the behavior of the
`qcc' library such as the background color, out-of-control points appearance,
 and many others.}
\usage{qcc.options(\dots)}
\arguments{
\item{\dots}{the option to be set or retrieved. See details.}
}
\details{
The available options are:
\describe{

\item{\code{exp.R.unscaled}}{a vector specifying, for each sample size, the expected value of the relative range (i.e. \eqn{R/\sigma}) for a normal distribution. This appears as \eqn{d_2} on most tables containing factors for the construction of control charts.}

\item{\code{se.R.unscaled}}{a vector specifying, for each sample size, the standard error of the relative range (i.e. \eqn{R/\sigma}) for a normal distribution. This appears as \eqn{d_3} on most tables containing factors for the construction of control charts.}

\item{\code{beyond.limits$pch}}{plotting character used to highlight points beyond control limits.}  

\item{\code{beyond.limits$col}}{color used to highlight points beyond control limits.}  

\item{\code{violating.runs$pch}}{plotting character used to highlight points violating runs.}

\item{\code{violating.runs$col}}{color used to highlight points violating runs.}  

\item{\code{run.length}}{the maximum value of a run before to signal a point as out of control.}  

\item{\code{bg.margin}}{background color used to draw the margin of the charts.} 

\item{\code{bg.figure}}{background color used to draw the figure of the charts.}

\item{\code{cex}}{character expansion used to draw labels, title, statistics at the bottom of plot produced by this package.}
} 
}
\value{
If the functions is called with no argument return a list of available options. 

If an option argument is provided the corresponding value is returned. 

If a value is associated with an option argument, such option is set and the list of updated option values is invisibly returned.
In this case the list \code{.qcc.options} is modified and any modification will remain in effect for the rest of the session.
}

\author{Luca Scrucca \email{luca@stat.unipg.it}}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
qcc.options()
qcc.options("cex")
qcc.options("cex"=2)
qcc.options(bg.margin="yellow")
}
\keyword{htest}
\keyword{hplot}

