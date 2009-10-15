#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

pareto.chart <- function(x, ylab = "Frequency", xlab, ylim, main, col = heat.colors(length(x)), ...)
{
  call <- match.call(expand.dots = TRUE)
  varname <- deparse(substitute(x))
  x <- as.table(x)
  if (length(dim(x))>1) 
     stop("only one-dimensional object (table, vector, etc.) may be provided")

  x <- sort(x, na.last=TRUE)
  missing <- is.na(names(x))
  x <- c(rev(x[!missing]), x[missing])
  missing <- is.na(names(x))
  cumsum.x <- cumsum(x[!missing]) 
  q <- seq(0, max(cumsum.x), length=5)

  if (missing(xlab)) 
     xlab <- ""
  if (missing(ylim)) 
     ylim <- c(0, max(cumsum.x)*1.05)
  if (missing(main)) 
     main <- paste("Pareto Chart for", varname)
  if (missing(col))
     col <- heat.colors(length(x))

  # set las and mar if not provided by user
  w <- max(sapply(names(x), nchar))
  if (is.null(call$las)) las <- 3 else las <- call$las
  if (is.null(call$mar))
     { if (las==1) mar <- c(0,1,0,2)  
       else        mar <- c(log(max(w),2),1,0,2) }
  else mar <- call$mar
  oldpar <- par(mar = par("mar")+mar, 
                las = las, 
                cex = qcc.options("cex"),
                no.readonly = TRUE)
  on.exit(par(oldpar))

  pc <- barplot(x, width = 1, space = 0.2, main = main, 
                ylim = ylim, ylab = ylab, xlab = xlab, col = col, ...)
  # adding line for percentage level overwrite bars...
  abline(h=q[2:5], col="lightgrey")
  # ... so we redraw bars (not nice but works!)
  rect(pc-0.5, rep(0,length(x)), pc+0.5, x, col = col)
  lines(pc[!missing], cumsum.x, type="b", cex=0.7)
  box()
  axis(4, at=q, las=3, labels=paste(seq(0,1,length=5)*100,"%",sep=""))
  mtext("Cumulative Percentage", 4, line=2.5, las=3) 

  tab <- cbind(x[!missing], cumsum.x, 
               x[!missing]/max(cumsum.x)*100, 
               cumsum.x/max(cumsum.x)*100) 
  colnames(tab) <- c("Frequency", "Cum.Freq.", 
                     "Percentage", "Cum.Percent.")
  names(dimnames(tab)) <- c("", paste("\nPareto chart analysis for", varname))
  return(tab)
}

