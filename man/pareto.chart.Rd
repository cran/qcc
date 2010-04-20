\name{pareto.chart}
\alias{pareto.chart}
\title{Pareto chart}
\description{Plot a Pareto chart.}
\usage{
pareto.chart(x, ylab = "Frequency", xlab, ylim, main, 
             col = heat.colors(length(x)), \dots)
}
\arguments{
\item{x}{a vector of values. \code{names(x)} are used for labelling the bars.}
\item{ylab}{a string specifying the label for the y-axis.}
\item{xlab}{a string specifying the label for the x-axis.}
\item{ylim}{a numeric vector specifying the limits for the y-axis.}
\item{main}{a string specifying the main title to appear on the plot.}   \item{col}{a value for the color, a vector of colors, or a palette for the bars. See the help for \code{colors} and \code{palette}.}
\item{\dots}{other graphical arguments to be passed to the \code{barplot} function.}
}
\value{Returns a matrix containing the descriptive statistics used to draw the Pareto chart.}
\details{A Pareto chart is a barplot where the categories are ordered in non increasing order, and a line is also added to show the cumulative sum.}
%\value{}
\references{
Montgomery, D.C. (2000) \emph{Introduction to Statistical Quality Control}, 4th ed. New York: John Wiley & Sons. 
}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
\seealso{\code{\link{barplot}}}
\examples{
defect <- c(80, 27, 66, 94, 33)
names(defect) <- c("price code", "schedule date", "supplier code", "contact num.", "part num.")
pareto.chart(defect, ylab = "Error frequency")
pareto.chart(defect, ylab = "Error frequency", xlab = "Error causes", las=1)
pareto.chart(defect, ylab = "Error frequency", col=rainbow(length(defect)))
}
\keyword{hplot}
