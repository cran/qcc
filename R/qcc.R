#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL CHARTS IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Luca Scrucca                                                   #
#              Department of Statistics                                       #
#              University of Perugia, ITALY                                   #
#              luca@stat.unipg.it                                             #
#                                                                             #
#-----------------------------------------------------------------------------#

#
#  Main function to create a `qcc' object
#

"qcc" <- function(data, type, sizes, center, std.dev, limits, target, data.name, labels, newdata, newsizes, newlabels, nsigmas = 3, confidence.level, rules = shewhart.rules, plot = TRUE, ...)
{
  call <- match.call()
  if (missing(data))
     stop("'data' argument is not specified")
  if (missing(data.name)) 
     data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(type))
     stop("'type' of chart is not specified")

  if (missing(sizes)) 
     { if (any(type==c("p", "np", "u")))
          stop(paste("sample 'sizes' must be given for a", type, "Chart"))
       else
          sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
     { if (length(sizes)==1)
          sizes <- rep(sizes, nrow(data))
       else if (length(sizes) != nrow(data))
                stop("sizes length doesn't match with data") }
         
  if (missing(labels))
     { if (is.null(rownames(data))) labels <- 1:nrow(data)
       else                         labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if (!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if (missing(center)) center <- stats$center
  
  if (missing(std.dev)) 
     { sd <- paste("sd.", type, sep = "")
       if (!exists(sd, mode="function"))
          stop(paste("function", sd, "is not defined!"))
       std.dev <- do.call(sd, list(data, sizes))  }
  else 
     { if (!mode(std.dev) == "numeric")
          stop("if provided the argument 'std.dev' must be a function or a numeric element")  }
    
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  
  object <- list(call = call, data.name = data.name, data = data, 
                 type = type, statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)

  # check for new data provided and update object
  if (!missing(newdata))
     { newdata.name <- deparse(substitute(newdata))
       newdata <- data.matrix(newdata)
       if (missing(newsizes))
          { if (any(type==c("p", "np", "u")))
               stop(paste("sample sizes must be given for a", type, "Chart"))
            else
               newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
       else
          { if (length(newsizes)==1)
               newsizes <- rep(newsizes, nrow(newdata))
            else if (length(newsizes) != nrow(newdata))
                     stop("newsizes length doesn't match with newdata") }
       stats <- paste("stats.", type, sep = "")
       if (!exists(stats, mode="function"))
          stop(paste("function", stats, "is not defined"))
       newstats <- do.call(stats, list(newdata, newsizes))$statistics
       if (missing(newlabels))
          { if (is.null(rownames(newdata)))
               { start <- length(statistics)
                 newlabels <- seq(start+1, start+length(newstats)) }
            else
               { newlabels <- rownames(newdata) }
          }
       names(newstats) <- newlabels
       object$new.statistics <- newstats
       object$newdata  <- newdata
       object$newsizes <- newsizes
       object$newdata.name <- newdata.name
       
       statistics <- c(statistics, newstats)
       sizes <- c(sizes, newsizes)
     }

  conf <- nsigmas
  if (!missing(confidence.level))
     conf <- confidence.level
  if (conf >= 1)
     { object$nsigmas <- conf }
  else
     if (conf > 0 & conf < 1)
        { object$confidence.level <- conf } 

  # get control limits      
  if (missing(limits))
     { limits <- paste("limits.", type, sep = "")
       if (!exists(limits, mode="function"))
          stop(paste("function", limits, "is not defined"))
       limits <- do.call(limits, list(center = center, std.dev = std.dev,
                                      sizes = sizes, conf = conf)) 
     }
  else 
     { if (!missing(std.dev))
          warning("'std.dev' is not used when limits is given")
       if (!is.numeric(limits))
          stop("'limits' must be a vector of length 2 or a 2-columns matrix")
       limits <- matrix(limits, ncol = 2)
       dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
     }

  lcl <- limits[,1]
  ucl <- limits[,2]
  object$limits <- limits

  if (!missing(target)) 
     if (!is.null(target))
        object$target <- target 

  if (is.function(rules)) violations <- rules(object)
  else                    violations <- NULL
  object$violations <- violations
         
  class(object) <- "qcc"
             
  if(plot) plot(object, ...) 

  return(object)
}

"summary.qcc" <- function(object, ...) print(object, ...)

"print.qcc" <- function(x, ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics), ...)
  sizes <- object$sizes
  if(!any(diff(sizes)))
     sizes <- sizes[1]
  if(length(sizes) == 1)
     cat("\nGroup sample size: ", format(sizes))
  else {
         cat("\nSummary of group sample sizes: ")
         tab <- table(sizes)
         print(matrix(c(as.numeric(names(tab)), tab), 
                      ncol = length(tab), byrow = TRUE, 
                      dimnames = list(c("  sizes", "  counts"),
                      character(length(tab)))), ...)
        }
  cat("\nNumber of groups: ", length(statistics))
  center <- object$center
  cat("\nCenter of group statistics: ", format(center))
  sd <- object$std.dev
  cat("\nStandard deviation: ", format(sd), "\n")

  newdata.name <- object$newdata.name
  new.statistics <- object$new.statistics
  if (!is.null(new.statistics)) 
     { cat(paste("\nSummary of group statistics in ", 
                 newdata.name, ":\n", sep = ""))
       print(summary(new.statistics), ...)
       newsizes <- object$newsizes
       if (!any(diff(newsizes)))
          newsizes <- newsizes[1]
       if (length(newsizes) == 1)
             cat("\nGroup sample size: ", format(newsizes))
       else 
           { cat("\nSummary of group sample sizes:\n")
             new.tab <- table(newsizes)
             print(matrix(c(as.numeric(names(new.tab)), new.tab),
                          ncol = length(new.tab), byrow = TRUE, 
                          dimnames = list(c("  sizes", "  counts"),
                                          character(length(new.tab)))), ...)
           }
       cat("\nNumber of groups: ", length(new.statistics), "\n")
     }
  
  if (!is.null(object$target))
     cat(paste("\nTarget value:", format(object$target), "\n"))

  limits <- object$limits
  if (!is.null(limits)) 
     { cat("\nControl limits:\n")
       print(limits, ...) }

  invisible()        
}


"plot.qcc" <- function(x, add.stats = TRUE, chart.all = TRUE, 
                       label.limits = c("LCL ", "UCL"),
                       title, xlab, ylab, ylim, axes.las = 0,
                       restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class `qcc' is required")

  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$new.statistics
  newdata.name <- object$newdata.name
  violations <- object$violations
  if (chart.all) 
     statistics <- c(stats, newstats)
  else
     if (is.null(newstats))
        statistics <- stats
     else
        statistics <- newstats   
  
  if (missing(title))
     { if (is.null(newstats))
            main.title <- paste(type, "Chart\nfor", data.name)
       else if (chart.all)
                 main.title <- paste(type, "Chart\nfor", data.name, 
                                     "and", newdata.name)
            else main.title <- paste(type, "Chart\nfor", newdata.name) }
  else main.title <- paste(title)
  
  indices <- 1:length(statistics)
  indcs <- c(indices - 1/2, max(indices) + 1/2)

  oldpar <- par(bg  = .qcc.options$bg.margin, 
                cex = .qcc.options$cex,
                mar = if(add.stats)
                           par("mar")+c(max(4, length(violations)), 0,0,0)
                      else par("mar"),
                no.readonly = TRUE)
  if (restore.par) on.exit(par(oldpar))

  # plot Shewhart chart
  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits, center),
       xlim = range(indcs), 
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE, main = main.title)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = .qcc.options$bg.figure)
  lines(indices, statistics, type = "b", pch=20) 
  axis(2, las = axes.las)
  nm.grp <- names(statistics)
  if (!length(nm.grp))
     nm.grp <- as.character(indices)
  axis(1, at = indices, labels = nm.grp, las = axes.las)
  box()

  if (length(center) == 1)
       abline(h = center)
  else lines(indcs, c(center, center[length(center)]), type="s")

  if (length(lcl) == 1) 
     { abline(h = lcl, lty = 2)
       abline(h = ucl, lty = 2) }
  else 
     { if (chart.all)
          { lines(indcs, c(lcl, lcl[length(lcl)]), type="s", lty = 2)
            lines(indcs, c(ucl, ucl[length(ucl)]), type="s", lty = 2) }
       else     
          { abline(h = lcl[length(lcl)], lty = 2)
            abline(h = ucl[length(ucl)], lty = 2) }
     }      
  text(rep(par("usr")[1], 2), c(lcl[1], ucl[1]), label.limits, 
       pos=4, col="darkgray")

  if (is.null(.qcc.options$beyond.limits))
     stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
  if (length(violations$beyond.limits))
     { v <- violations$beyond.limits
       points(indices[v], statistics[v], 
              col=.qcc.options$beyond.limits$col, 
              pch=.qcc.options$beyond.limits$pch) }
  if (is.null(.qcc.options$violating.runs))
     stop(".qcc.options$violating.runs undefined. See help(qcc.options).")
  if (length(violations$violating.runs))
     { v <- violations$violating.runs
       points(indices[v], statistics[v], 
              col=.qcc.options$violating.runs$col, 
              pch=.qcc.options$violating.runs$pch) }
  
  if (chart.all & (!is.null(newstats)))
     { len.obj.stats <- length(object$statistics)
       len.new.stats <- length(statistics) - len.obj.stats
       abline(v = len.obj.stats + 0.5, lty = 3)
       mtext(paste("Calibration data in", data.name), cex=0.8, 
             at = len.obj.stats/2, adj = 0.5)
       mtext(paste("New data in", object$newdata.name), cex=0.8, 
             at = len.obj.stats + len.new.stats/2, adj = 0.5)
     }

  if (add.stats) 
     { # computes the x margins of the figure region
       plt <- par()$plt; usr <- par()$usr
       px <- diff(usr[1:2])/diff(plt[1:2])
       xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
       at.col <- xfig[1] + diff(xfig[1:2])*c(0.10, 0.40, 0.65)
       # write info at bottom
       mtext(paste("Number of groups = ", length(statistics), sep = ""), 
             side = 1, line = 5, adj = 0, font=2, at=at.col[1])
       center <- object$center
       if (length(center) == 1)
          mtext(paste("Center = ", signif(center, options()$digits), sep = ""),
                side = 1, line = 6, adj = 0, font=2, at=at.col[1])
       else 
          mtext("Center is variable", 
                side = 1, line = 6, adj = 0, font=2, at=at.col[1])
       mtext(paste("StdDev = ", signif(std.dev, options()$digits), sep = ""),
             side = 1, line = 7, adj = 0, font=2, at=at.col[1])
       target <- object$target
       if (!is.null(target))
          { if (length(target) == 1)
                mtext(paste("Target = ", signif(target, options()$digits), 
                            sep = ""),
                      side = 1, line = 5, adj = 0, font=2, at=at.col[2])
            else 
                mtext("Target is variable", 
                      side = 1, line = 5, adj = 0, font=2, at=at.col[2])
          }
       if (length(unique(lcl)) == 1)
          mtext(paste("LCL = ", signif(lcl, options()$digits), sep = ""), 
                side = 1, line = 6, adj = 0, font=2, at=at.col[2])
       else 
          mtext("LCL is variable", side = 1, line = 6, adj = 0, font=2,
                at=at.col[2])
       if (length(unique(ucl)) == 1)
          mtext(paste("UCL = ", signif(ucl, options()$digits), sep = ""),
                side = 1, line = 7, adj = 0, font=2, at=at.col[2])
       else 
          mtext("UCL is variable", side = 1, line = 7, adj = 0, font=2,
                at=at.col[2])
       if (!is.null(violations))
          { mtext(paste("Number beyond limits =",
                        length(unique(violations$beyond.limits))), 
                  side = 1, line = 6, adj = 0, font=2,
                  at = at.col[3])
            mtext(paste("Number violating runs =",
                        length(unique(violations$violating.runs))), 
                  side = 1, line = 7, adj = 0, font=2,
                  at = at.col[3])
           }
     }

  invisible() 
}

#
#  Functions used to compute Shewhart charts statistics
#

# xbar

"stats.xbar" <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  statistics <- apply(data, 1, mean, na.rm=TRUE)
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

"sd.xbar" <- function(data, sizes, std.dev, equal.sd = TRUE)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if (any(sizes == 1))
     stop("group sizes must be larger than one")
  if (missing(std.dev))
     var.within <- apply(data, 1, var, na.rm=TRUE)
  else 
     var.within <- std.dev^2
  var.df <- sum(sizes - 1)
  c4 <- function(n)
        sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2))
  if (equal.sd) 
     { std.dev <- sqrt(sum((sizes - 1) * var.within)/var.df) / c4(var.df + 1) }
  else 
     { c <- c4(sizes)/(1 - c4(sizes)^2)
       std.dev <- sum(c * sqrt(var.within))/sum(c * c4(sizes)) }
  return(std.dev)
}

"limits.xbar" <- function(center, std.dev, sizes, conf)
{
  if (!any(diff(sizes)))
     sizes <- sizes[1]
  se.stats <- std.dev/sqrt(sizes)
  if (conf >= 1) 
     { lcl <- center - conf * se.stats
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid 'conf' argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# S chart

"stats.S" <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  statistics <- sqrt(apply(data, 1, var, na.rm=TRUE))
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

"sd.S" <- function(data, sizes, std.dev, equal.sd = TRUE)
{
  sd.xbar(data, sizes, std.dev, equal.sd)
}

"limits.S" <- function(center, std.dev, sizes, conf)
{
  if (!any(diff(sizes)))
     sizes <- sizes[1]
  c4 <- function(n)
        sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2))
  se.stats <- std.dev * sqrt(1 - c4(sizes)^2)
  if (conf >= 1) 
     { lcl <- pmax(0, center - conf * se.stats)
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- std.dev * sqrt(qchisq(1 - (1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
            lcl <- std.dev * sqrt(qchisq((1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
          }
          else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  limits
}

# R Chart 

"stats.R" <- function(data, sizes)
{
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  statistics <- apply(data, 1, function(x) diff(range(x, na.rm=TRUE)))
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

"sd.R" <- function(data, sizes, std.dev, equal.sd = TRUE)
{
  sd.xbar(data, sizes, std.dev, equal.sd)
}

"limits.R" <- function(center, std.dev, sizes, conf)
{
  if (!any(diff(sizes)))
     sizes <- sizes[1]
  se.R.unscaled <- .qcc.options$se.R.unscaled
  Rtab <- length(se.R.unscaled)
  if (conf >= 1) 
     { if (any(sizes > Rtab))
          stop(paste("group size must be less than", 
                      Rtab + 1, "when giving nsigmas"))
       se.R <- se.R.unscaled[sizes] * std.dev
       lcl <- pmax(0, center - conf * se.R)
       ucl <- center + conf * se.R
     }
  else 
     { if (conf > 0 && conf < 1) 
          { ucl <- qtukey(1 - (1 - conf)/2, sizes, 1e100) * std.dev
            lcl <- qtukey((1 - conf)/2, sizes, 1e100) * std.dev
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# xbar Chart for one-at-time data

"stats.xbar.one" <- function(data, sizes)
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

"sd.xbar.one" <- function(data, sizes, std.dev, k=2)
{
  data <- as.vector(data)
  n <- length(data)
  d2 <- .qcc.options$exp.R.unscaled
  if (is.null(d2))
     stop(".qcc.options$exp.R.unscaled is null")
  if (missing(std.dev))
     { d <- 0
       for (j in k:n)
           d <- d+abs(diff(range(data[c(j:(j-k+1))])))
       var.within <- d/(n-k+1) 
       std.dev <- var.within/d2[k] }
  else 
     { std.dev <- std.dev/d2[k] }
  return(std.dev)
}

       
"limits.xbar.one" <- function(center, std.dev, sizes, conf)
{
  se.stats <- std.dev
  if (conf >= 1) 
     { lcl <- center - conf * se.stats
       ucl <- center + conf * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}


# p Chart

"stats.p" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  list(statistics = data/sizes, center = pbar)
}

"sd.p" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(pbar * (1 - pbar))
  return(std.dev)
}

"limits.p" <- function(center, std.dev, sizes, conf)
{ 
  # use normal approximation for computing control limits
  limits <- limits.xbar(center, std.dev, sizes, conf)
  limits[limits < 0] <- 0
  limits[limits > 1] <- 1
  return(limits)
}

# np Chart

"stats.np" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  center <- sizes * pbar
  if (length(unique(center)) == 1)
     center <- center[1]
  list(statistics = data, center = center)
}

"sd.np" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(sizes * pbar * (1 - pbar))
  if (length(unique(std.dev)) == 1)
     std.dev <- std.dev[1]
  return(std.dev)
}

"limits.np" <- function(center, std.dev, sizes, conf)
{ 
  # use normal approximation for computing control limits
  if (length(unique(sizes)) == 1)
     sizes <- sizes[1]
  limits <- limits.xbar(center, std.dev*sqrt(sizes), sizes, conf)
  limits[limits < 0] <- 0
  if (length(sizes)==1)
      limits[limits > sizes] <- sizes
  else
      limits[limits > sizes] <- sizes[limits > sizes]
  return(limits)
}

# c Chart

"stats.c" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
     stop("all sizes must be be equal for a c chart")
  statistics <- data
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

"sd.c" <- function(data, sizes)
{
  data <- as.vector(data)
  std.dev <- sqrt(mean(data))
  return(std.dev)
}

"limits.c" <- function(center, std.dev, sizes, conf)
{
  if (conf >= 1) 
     { lcl <- center - conf * sqrt(center)
       lcl[lcl < 0] <- 0
       ucl <- center + conf * sqrt(center)
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- qpois(1 - (1 - conf)/2, center)
            lcl <- qpois((1 - conf)/2, center)
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# u Chart

"stats.u" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  statistics <- data/sizes
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

"sd.u" <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  std.dev <- sqrt(sum(sizes * data)/sum(sizes))
  return(std.dev)
}

"limits.u" <- function(center, std.dev, sizes, conf)
{
  sizes <- as.vector(sizes)
  if (!any(diff(sizes)))
     sizes <- sizes[1]
  if (conf > 0 & conf < 1)  
     conf <- qnorm(1 - (1 - conf)/2)
  lcl <- center - conf * sqrt(center/sizes)
  lcl[lcl < 0] <- 0
  ucl <- center + conf * sqrt(center/sizes)
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

#
# Functions used to signal points out of control 
#

"shewhart.rules" <- function(object, run.length = qcc.options("run.length"))
{
# Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object)
  vr <- violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}

"beyond.limits" <- function(object)
{
# Return cases beyond limits
  statistics <- c(object$statistics, object$new.statistics) 
  lcl <- object$limits[,1]
  ucl <- object$limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

"violating.runs" <- function(object, run.length = qcc.options("run.length"))
{
# Return indices of points violating runs
  center <- object$center
  statistics <- c(object$statistics, object$new.statistics)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  index.lenghts <- (1:length(runs$lengths))[runs$lengths >= run.length]
  index.stats   <- 1:length(statistics)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] 
                    - (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] 
                    - (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
     { for (i in 1:length(vbeg.above))
           violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
     { for (i in 1:length(vbeg.below))
           violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}

#-------------------------------------------------------------------#
#                                                                   #
#                      CUSUM CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

"cusum" <- function(object, ...) UseMethod("cusum")

"cusum.default" <- function(object, ...) object

"cusum.qcc" <- function(object, decision.int = 5, se.shift = 1, label.bounds = c("LDB", "UDB"), add.stats = TRUE, chart.all = TRUE, ylim = NULL, axes.las = 0, restore.par = TRUE, ...)
{
  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")

  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- if (is.null(object$target)) object$center
            else                        object$target
  data.name <- object$data.name

  set.to.zero <- 0
  if (chart.all) 
     { statistics <- c(object$statistics, object$new.statistics)
       sizes <- c(object$sizes, object$newsizes)
       set.to.zero <- length(object$statistics)  }
  else 
     { if (is.null(object$new.statistics)) 
          { statistics <- object$statistics
            sizes <- object$sizes }
       else 
          { statistics <- object$new.statistics
            sizes <- object$newsizes  }
     }

  n <- length(statistics)
  z <- (statistics - center)/(std.dev/sqrt(sizes))
  ldb <-  - decision.int
  udb <- decision.int
  indices <- 1:length(statistics)
  #
  z.f <- z - se.shift/2
  cusum.pos <- rep(NA, n)
  cusum.pos[1] <- max(0, z.f[1])
  for (i in 2:n)
      cusum.pos[i] <- max(0, cusum.pos[i-1]+z.f[i])
  #
  z.f <- z + se.shift/2
  cusum.neg <- rep(NA, n)
  cusum.neg[1] <- max(0, -z.f[1])
  for (i in 2:n)
      cusum.neg[i] <- max(0, cusum.neg[i-1]-z.f[i])
  cusum.neg <- -cusum.neg

  if (is.null(object$new.statistics))
     { main.title <- paste("Cusum Chart\nfor", data.name) }
  else
     { if (chart.all)
            main.title <- paste("Cusum Chart\nfor", data.name,
                                "and", object$newdata.name)
       else main.title <- paste("Cusum Chart\nfor", object$newdata.name)
     }

  oldpar <- par(bg  = .qcc.options$bg.margin, 
                cex = .qcc.options$cex,
                mar = if(add.stats) pmax(par("mar"), c(9,0,0,0))
                      else          par("mar"),
                no.readonly = TRUE)
  if (restore.par) on.exit(par(oldpar))

  plot(indices, cusum.pos, type="n",
       ylim = range(cusum.pos, cusum.neg, ldb, udb, ylim),
       ylab = "", xlab = "Group", axes = FALSE, main=main.title)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = .qcc.options$bg.figure)
  abline(h = 0, lwd = 2)
  abline(h = c(ldb, udb), lty = 2)
  nm.grp <- names(statistics)
  if (!length(nm.grp))
     nm.grp <- as.character(indices)
  axis(1, at = indices, labels = nm.grp, las = axes.las)
  axis(2, at=round(seq(par("usr")[3],par("usr")[4],by=1)), las = axes.las)
  box()
  
  if (chart.all) 
     { lind <- length(indices)
       lines(indices[1:set.to.zero], 
             cusum.pos[1:set.to.zero], 
             type = "b", pch=20)
       lines(indices[(set.to.zero + 1):lind], 
             cusum.pos[(set.to.zero + 1):lind],
             type = "b", pch=20)
       lines(indices[1:set.to.zero], 
             cusum.neg[1:set.to.zero], 
             type = "b", pch=20)
       lines(indices[(set.to.zero + 1):lind], 
             cusum.neg[(set.to.zero + 1):lind],
             type = "b", pch=20)
     }
  else 
     { lines(indices, cusum.pos, type = "b", pch=20)
       lines(indices, cusum.neg, type = "b", pch=20) }

  mtext("Cumulative Sum", line=3, side=2)
  mtext("Above Target", srt=90, line=2, side=2, at=0+par("usr")[4]/2)
  mtext("Below Target", srt=90, line=2, side=2, at=0+par("usr")[3]/2)
  text(rep(par("usr")[1],2), c(ldb, udb), label.bounds, pos=4, col="darkgray")

  if (is.null(.qcc.options$beyond.limits))
     stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
  beyond.bounds <- list(upper=(cusum.pos > udb), lower=(cusum.neg < ldb))
  n.beyond.bounds <- sum(beyond.bounds$upper) + sum(beyond.bounds$lower)
  if (n.beyond.bounds > 0)
     { v <- beyond.bounds$upper
       points(indices[v], cusum.pos[v], 
              col=.qcc.options$beyond.limits$col, 
              pch=.qcc.options$beyond.limits$pch)
       v <- beyond.bounds$lower
       points(indices[v], cusum.neg[v], 
              col=.qcc.options$beyond.limits$col, 
              pch=.qcc.options$beyond.limits$pch) }

  if (chart.all & !is.null(object$new.statistics))
     { len.obj.stats <- length(object$statistics)
       len.new.stats <- length(statistics) - len.obj.stats
       abline(v = len.obj.stats + 0.5, lty = 3)
       mtext(paste("Calibration Data in", data.name), cex=0.8,
             at = len.obj.stats/2, line = 0, adj = 0.5)
       mtext(paste("New Data in", object$newdata.name), cex=0.8, 
             at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
     }

    if (add.stats) 
       { # computes the x margins of the figure region
         plt <- par()$plt; usr <- par()$usr
         px <- diff(usr[1:2])/diff(plt[1:2])
         xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
         at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.50)
         # write info at bottom
         mtext(paste("Number of groups = ", length(statistics), sep = ""),
               side = 1, line = 5, adj = 0, font = 2, at=at.col[1])
         if (length(center) == 1)
            mtext(paste("Target = ", 
                        signif(center,options()$digits), sep = ""),
                  side = 1, line = 6, adj = 0, font = 2, at=at.col[1])
         else 
            mtext("Target is variable", side = 1, line = 6, adj = 0)
         mtext(paste("StdDev = ", signif(std.dev, options()$digits), sep = ""),
               side = 1, line = 7, adj = 0, font = 2, at=at.col[1])
         mtext(paste("Decision boundaries (std. err.) =", decision.int),
               side = 1, line = 5, adj = 0, font = 2, at=at.col[2])
         mtext(paste("Shift detection (std. err.) =", 
                     signif(se.shift, digits = options()$digits)), 
               side = 1, line = 6, adj = 0, font = 2, at=at.col[2])
         mtext(paste("No. of points beyond boundaries =", n.beyond.bounds),
               side = 1, line = 7, adj = 0, font = 2, at=at.col[2])
       }
            
  object$cusum <- list(x = indices, 
                       pos = cusum.pos, 
                       neg = cusum.neg,                 
                       decision.int = decision.int, 
                       se.shift = se.shift)
                                                        
  class(object) <- c("cusum", "qcc")
  invisible(object)
}


#-------------------------------------------------------------------#
#                                                                   #
#                       EWMA CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

"ewmaSmooth" <- function(x, y, lambda=0.20, start, ...)
{
#
# Exponential-Weighted Moving Average 
# 
# Return smooth values based on 
# 
# z_t = lambda*y_t + (1-lambda)*z_t-1      
# 
# where 0<= lambda <=1 is the parameter which controls the weights applied 
# to the data, and start is the starting value.
# Returns a list with elements:
# x = ordered x-values
# y = smoothed fitted values of y
# 
  if (length(y)!=length(x))
     stop("x and y must have the same length!")
  if (abs(lambda)>1)
     stop("lambda parameter must be between 0 and 1")
  ord <- order(x) 
  x <- x[ord]
  y <- y[ord]
  n <- length(y)
  if (missing(start)) start <- y[1]
  
  S1 <- diag(rep(1,n))
  for (i in 1:(n-1))
      {for (j in i:n)
           { S1[j,i] <- (1-lambda)^(j-i) }}
            
  S2 <- (1-lambda)^seq(1,n)
  z <- lambda*(S1%*%y) + S2*start
  list(x=x, y=z, lambda=lambda, start=start)
}

"ewma" <- function (object, ...) UseMethod("ewma")

# "ewma.default" <- function (object, ...) ewmaSmooth(object, ...)

"ewma.qcc" <- function(object, lambda = 0.2, nsigmas = object$nsigmas, 
                       add.stats = TRUE, xlab, ylab, ylim = NULL, 
                       axes.las = 0, restore.par = TRUE, ...)
{
# Exponential-Weighted Moving Average Plot

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class `qcc' is required")
  if (is.null(nsigmas))
     stop("the object of class `qcc' does not contain nsigmas, so it must be provided as argument.")

  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- if (is.null(object$target)) object$center
            else                        object$target
  data.name <- object$data.name

  statistics <- c(object$statistics, object$new.statistics)
  sizes <- c(object$sizes, object$newsizes)

  n <- length(statistics)
  indices <- 1:length(statistics)
         
  ewma <- ewmaSmooth(indices, statistics, lambda=lambda, start=center)
  sigma2 <- std.dev^2/sizes * 
            ((lambda/(2-lambda))*(1-(1-lambda)^(2*(1:n))))

  ucl <- center + nsigmas*sqrt(sigma2)
  lcl <- center - nsigmas*sqrt(sigma2)

  if (is.null(object$new.statistics))
     main.title <- paste("EWMA Chart\nfor", data.name)
  else
     main.title <- paste("EWMA Chart\nfor", data.name, 
                         "and", object$newdata.name)
     
  oldpar <- par(bg  = .qcc.options$bg.margin, 
                cex = .qcc.options$cex,
                mar = if(add.stats) pmax(par("mar"), c(9,0,0,0))
                      else          par("mar"),
                no.readonly = TRUE)
  if (restore.par) on.exit(par(oldpar))

  plot(indices, statistics, type="n",
       ylim = range(statistics, lcl, ucl, center, ylim), 
       xlim = range(indices), 
       ylab = ifelse(missing(ylab), "Group Summary Statistics", ylab),
       xlab = ifelse(missing(xlab), "Group", xlab),
       axes = FALSE, main = main.title)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = .qcc.options$bg.figure)
  axis(1, at = indices, labels = indices, las = axes.las)
  axis(2, las = axes.las)
  box()
  points(indices, statistics, pch=20) 
  lines(ewma$x, ewma$y)
  lines(indices, ucl, lty=2)
  lines(indices, lcl, lty=2)
  abline(h=center, lty=1)
  
  if (!is.null(object$new.statistics))
     { len.obj.stats <- length(object$statistics)
       len.new.stats <- length(statistics) - len.obj.stats
       abline(v = len.obj.stats + 0.5, lty = 3)
       mtext(paste("Calibration Data in", data.name), cex=0.8,
             at = len.obj.stats/2, line = 0, adj = 0.5)
       mtext(paste("New Data in", object$newdata.name), cex=0.8, 
             at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
     }

  if (add.stats) 
     { # computes the x margins of the figure region
       plt <- par()$plt; usr <- par()$usr
       px <- diff(usr[1:2])/diff(plt[1:2])
       xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
       at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.60)
       # write info at bottom
       mtext(paste("Number of groups = ", length(statistics), sep = ""),
             side = 1, line = 5, adj = 0, font = 2, at=at.col[1])
       if (length(center) == 1)
          mtext(paste("Target = ", 
                      signif(center,options()$digits), sep = ""),
                side = 1, line = 6, adj = 0, font = 2, at=at.col[1])
       else 
          mtext("Target is variable", side = 1, line = 6, adj = 0)
       mtext(paste("StdDev = ", signif(std.dev, options()$digits), sep = ""),
             side = 1, line = 7, adj = 0, font = 2, at=at.col[1])
       mtext(paste("Smoothing parameter = ", signif(lambda,5)),
             side = 1, line = 5, adj = 0, font = 2, at = at.col[2])
       mtext(paste("Control limits at ", nsigmas, "*sigma", sep=""),
             side = 1, line = 6, adj = 0, font = 2, at = at.col[2])
       # mtext(substitute(lambda==l, list(l=signif(lambda,5))),
       #       side = 1, line = 6, adj = 0, font = 2,
       #       at = par("usr")[1] + 0.5 * diff(par("usr")[1:2]))
     }

  object$ewma <- list(x = ewma$x, 
                      y = { y <- as.vector(ewma$y)
                            names(y) <- c(names(object$statistics),
                                          names(object$new.statistics))
                            y },                       
                      sigma = sqrt(sigma2), nsigmas = nsigmas, 
                      limits = { limits <- cbind(lcl,ucl)
                                 colnames(limits) <- c("LCL", "UCL")
                                 limits })
  class(object) <- c("ewma", "qcc")
  invisible(object)
} 



#-------------------------------------------------------------------#
#                                                                   #
#          Operating Characteristic Function                        #
#                                                                   #
#-------------------------------------------------------------------#

"oc.curves" <- function(object, ...)
{
# Draws the operating characteristic curves for the qcc object 

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")

  size <- unique(object$sizes)
  if (length(size)>1)
     stop("Operating characteristic curves available only for equal sample sizes!")

  if (object$type=="p" | object$type=="np")
     beta <- oc.curves.p(object, ...)
  else
     if (object$type=="c" | object$type=="u")
        beta <- oc.curves.c(object, ...)
     else
        if (object$type=="xbar")
           beta <- oc.curves.xbar(object, ...)
     else
           stop("Operating characteristic curves not available for this type of chart.")

  invisible(beta)
}

"oc.curves.xbar" <- function(object, n, c = seq(0, 5, length=101), nsigmas = object$nsigmas, identify=FALSE)
{
# Draw the operating-characteristic curves for the xbar-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a shift of c*sigma in the mean on the first sample following the shift.

  type <- object$type
  if (!(object$type=="xbar"))
     stop("not a `qcc' object of type \"xbar\".")

  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
  if (missing(n))
     n <- unique(c(size, c(1,5,10,15,20)))

  beta <- matrix(NA, length(n), length(c))
  for (i in 1:length(n))
      beta[i,] <- pnorm(nsigmas-c*sqrt(n[i])) - pnorm(-nsigmas-c*sqrt(n[i]))
  rownames(beta) <- paste("n=",n,sep="")
  colnames(beta) <- c

  oldpar <- par(bg=.qcc.options$bg.margin, no.readonly = TRUE)
  on.exit(par(oldpar))

  plot(c, beta[1,], type="n",
       ylim = c(0,1), xlim = c(0,max(c)),
       xlab = "Process shift (std.dev)",
       ylab = "Prob. type II error ",
       main = paste("OC curves for", object$type, "chart"))
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = .qcc.options$bg.figure)
  for (i in 1:length(n))
      lines(c, beta[i,], type = "l", lty=i)
  beta <- t(beta)
  names(dimnames(beta)) <- c("shift (std.dev)", "sample size")
  
  if (identify)
     { cs <- rep(c,length(n))
       betas <- as.vector(beta)
       labels <- paste("c=", formatC(cs, 2, flag="-"), 
                       ": beta=", formatC(betas, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-betas), 2, flag="-"), sep="")
       i <- identify(cs, betas, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  else
     { legend(max(c), 1, legend = paste("n =", n), 
              bg = .qcc.options$bg.figure,
              lty = 1:length(n), xjust = 1, yjust = 1)
     }  
  invisible(beta)
}

"oc.curves.p" <- function(object, nsigmas = object$nsigmas, identify=FALSE)
{
  type <- object$type
  if (!(object$type=="p" | object$type=="np"))
     stop("not a `qcc' object of type \"p\" or \"np\".")
     
  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample sizes!")
     
  if (is.null(object$limits))
     stop("the `qcc' object does not have control limits!")
  limits <- object$limits
  p <- seq(0, 1, length=101)
  
  if (object$type=="p") 
     { UCL <- min(floor(size*limits[,2]), size)
       LCL <- max(floor(size*limits[,1]), 0) }
  else
     { UCL <- min(floor(limits[,2]), size)
       LCL <- max(floor(limits[,1]), 0) }
  beta <- pbinom(UCL, size, p) - pbinom(LCL-1, size, p)
  names(beta) <- p
  
  oldpar <- par(bg=.qcc.options$bg.margin, no.readonly = TRUE)
  on.exit(par(oldpar))

  plot(p, beta, type = "n", 
       ylim = c(0,1), xlim = c(0,1),
       main = paste("OC curves for", object$type, "chart"), 
       xlab = expression(p), 
       ylab = "Prob. type II error ")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = .qcc.options$bg.figure)
  lines(p, beta)
  lines(rep(p[which.max(beta)], 2), c(0, max(beta)), lty = 2)
  
  warning("Some computed values for the type II error have been rounded due to the discreteness of the binomial distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("p=", formatC(p, 2, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(p, beta, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  invisible(beta)  
}

"oc.curves.c" <- function(object, nsigmas = object$nsigmas, identify=FALSE)
{
  type <- object$type
  if (!(object$type=="c" | object$type=="u"))
     stop("not a `qcc' object of type \"c\" or \"u\".")
     
  size <- unique(object$sizes)
  if (length(size) > 1)
     stop("Operating characteristic curves available only for equal sample size!")
        
  if (is.null(object$limits))
     stop("the `qcc' object does not have control limits!")
  limits <- object$limits
  CL  <- object$center
  std.dev <- object$std.dev
  if (object$type=="c") 
     { max.lambda <- ceiling(CL+10*std.dev)
       UCL <- floor(limits[1,2])
       LCL <- floor(limits[1,1])
     }      
  else
     { max.lambda <- ceiling(CL*size+10*std.dev)[1]
       UCL <- floor(size*limits[1,2])
       LCL <- floor(size*limits[1,1])
     }  
  lambda <- seq(0, max.lambda)
  beta <- ppois(UCL, lambda) - ppois(LCL-1, lambda)
  names(beta) <- lambda

  oldpar <- par(bg=.qcc.options$bg.margin, no.readonly = TRUE)
  on.exit(par(oldpar))
  
  plot(lambda, beta, type = "n", 
       ylim = c(0,1), xlim = range(lambda),
       main = paste("OC curves for", object$type, "chart"), 
       xlab = "Mean", 
       ylab = "Prob. type II error ")
  lines(lambda, beta)     
  lines(rep(lambda[which.max(beta)], 2), c(0, max(beta)), lty = 2)

  warning("Some computed values for the type II error have been rounded due to the discreteness of the Poisson distribution. Thus, some ARL values might be meaningless.")

  if (identify)
     { labels <- paste("lambda=", formatC(lambda, 0, flag="-"), 
                       ": beta=", formatC(beta, 4, flag="-"), 
                       ", ARL=", formatC(1/(1-beta), 2, flag="-"), sep="")
       i <- identify(lambda, beta, labels, pos=4, offset=0.2)
       apply(as.matrix(labels[i$ind]), 1, cat, "\n")
     }
  invisible(beta)  
}


#-------------------------------------------------------------------#
#                                                                   #
#    Process Capability Analysis                                    #
#                                                                   #
#-------------------------------------------------------------------#

"process.capability" <- function(object, spec.limits, target, std.dev, nsigmas, confidence.level = 0.95, breaks="scott", add.stats=TRUE, print=TRUE, restore.par=TRUE)
{
# Computes process capability indices for a qcc object of type "xbar" 
# and plot the histogram

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!(object$type=="xbar" | object$type=="xbar.one"))
     stop("Process Capability Analysis only available for charts type \"xbar\" and \"xbar.one\" charts")

  x <- as.vector(object$data)
  x <- x[!is.na(x)]
  sizes <- object$sizes
  center <- object$center
  if (missing(std.dev))
     std.dev <- object$std.dev
  n <- length(x)
  title <- paste("Process Capability Analysis\nfor", object$data.name)
  if (missing(spec.limits))
     stop("specification limits must be provided")
  if (!length(spec.limits)==2)
     stop("wrong specification limits format")
  LSL <- min(spec.limits, na.rm=TRUE)
  USL <- max(spec.limits, na.rm=TRUE)
  
  if (missing(target))
     { if (is.null(object$target))
          target <- mean(spec.limits, na.rm=TRUE)
       else
          target <- object$target }
  if (target < LSL | target > USL)
     warning("target value is not within specification limits...") 
         
  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas
  
  if (confidence.level < 0 | confidence.level > 1)
     stop("the argument confidence.level must be a value between 0 and 1") 

  # computes process capability indices
  
  Cp <- (USL - LSL) / (2*nsigmas*std.dev)
  Cp.u <- (USL-center)/(nsigmas*std.dev)
  Cp.l <- (center-LSL)/(nsigmas*std.dev)
  Cp.k <- min(Cp.u, Cp.l)
  # Cpm <- (USL - LSL) / (2*nsigmas*sqrt(sum((x-target)^2)/(n-1)))
  Cpm <- Cp / sqrt(1+((center-target)/std.dev)^2)

  # compute confidence limits 
  
  alpha <- 1-confidence.level
  Cp.limits   <- Cp*sqrt(qchisq(c(alpha/2, 1-alpha/2), n-1)/(n-1))
  Cp.u.limits <- Cp.u*(1+c(-1,1)*qnorm(confidence.level)*
                       sqrt(1/(9*n*Cp.u^2)+1/(2*(n-1))))
  Cp.l.limits <- Cp.l*(1+c(-1,1)*qnorm(confidence.level)*
                       sqrt(1/(9*n*Cp.l^2)+1/(2*(n-1))))
  Cp.k.limits <- Cp.k*(1+c(-1,1)*qnorm(1-alpha/2)*
                       sqrt(1/(9*n*Cp.k^2)+1/(2*(n-1))))
  df <- n*(1+((center-target)/std.dev)^2)/
          (1+2*((center-target)/std.dev)^2)
  Cpm.limits <- Cpm*sqrt(qchisq(c(alpha/2, 1-alpha/2), df)/df)
  names(Cp.limits) <- names(Cp.k.limits) <- names(Cpm.limits) <- 
    c(paste(round(100*alpha/2, 1), "%", sep=""),
      paste(round(100*(1-alpha/2), 1), "%", sep=""))
    
  exp.LSL <- pnorm((LSL-center)/std.dev) * 100
  if (exp.LSL < 0.01) exp.LSL <- 0
  exp.USL <- (1-pnorm((USL-center)/std.dev)) * 100
  if (exp.USL < 0.01) exp.USL <- 0
  obs.LSL <- sum(x<LSL)/n * 100
  obs.USL <- sum(x>USL)/n * 100
  xlim <- range(x, USL, LSL, target)
  xlim <- xlim+diff(xlim)*c(-0.1,0.1)
  xx <- seq(min(xlim), max(xlim), length=100)
  dx <- dnorm(xx, center, std.dev)
  h <- hist(x, breaks = breaks, plot=FALSE) # compute histogram
  ylim <- range(h$density, dx)
  ylim <- ylim+diff(ylim)*c(0,0.05)

  tab <- cbind(c(Cp, Cp.l, Cp.u, Cp.k, Cpm),
               rbind(Cp.limits, Cp.l.limits, Cp.u.limits, 
                     Cp.k.limits, Cpm.limits))
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")
  colnames(tab) <- c("Value", names(Cp.limits))

  oldpar <- par(bg  = .qcc.options$bg.margin, 
                cex = .qcc.options$cex,
                mar = if(add.stats) 
                            c(9+is.null(center)*-1, 2, 4, 2) + 0.1
                      else  par("mar"),
                no.readonly = TRUE)
  if (restore.par) on.exit(par(oldpar))

  plot(0, 0, type="n", xlim = xlim, ylim = ylim,
       axes = FALSE, ylab="", xlab = "", main = title)
  usr <- par()$usr
  rect(usr[1], usr[3], usr[2], usr[4], col=.qcc.options$bg.figure)
  axis(1); box()
  plot(h, add=TRUE, freq=FALSE) # draw histogram
  # add graphical info
  abline(v=c(LSL,USL), col=2, lty=3, lwd=2)
  text(LSL, usr[4], "LSL", pos=1, col="darkgray", cex=0.8)
  text(USL, usr[4], "USL", pos=1, col="darkgray", cex=0.8)
  if (!is.null(target))
     { abline(v=target, col=2, lty=2, lwd=2)
       text(target, usr[4], "Target", pos=1, col="darkgray", cex=0.8) }
  lines(xx, dx, lty=2)
  
  if (add.stats) 
     { 
       # computes the x margins of the figure region
       plt <- par()$plt
       px <- diff(usr[1:2])/diff(plt[1:2])
       xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
       at.col <- xfig[1] + diff(xfig[1:2])*c(0.07, 0.35, 0.56, 0.75)
       # write info at bottom
       #-- 
       mtext(paste("Number of obs = ", n, sep = ""), 
             side = 1, line = 3, adj = 0, font=2, at = at.col[1])
       mtext(paste("Center = ", signif(center, options()$digits), sep = ""), 
             side = 1, line = 4, adj = 0, font=2, at = at.col[1])
       mtext(paste("StdDev = ", signif(std.dev, options()$digits), sep = ""), 
             side = 1, line = 5, adj = 0, font=2, at = at.col[1])
       #--
       if (!is.null(target))
          msg <- paste("Target = ", signif(target, options()$digits), sep = "")
       else
          msg <- paste("Target = ", sep = "")     
       mtext(msg, side = 1, line = 3, adj = 0, font=2, at=at.col[2]) 
       mtext(paste("LSL = ", signif(LSL, options()$digits), sep = ""), 
             side = 1, line = 4, adj = 0, font=2, at = at.col[2])
       mtext(paste("USL = ", signif(USL, options()$digits), sep = ""), 
             side = 1, line = 5, adj = 0, font=2, at = at.col[2])
       #--
       mtext(paste("Cp     = ", signif(Cp, 3), sep = ""), 
             side = 1, line = 3, adj = 0, font=2, at = at.col[3])
       mtext(paste("Cp_l  = ", signif(Cp.l, 3), sep = ""), 
             side = 1, line = 4, adj = 0, font=2, at = at.col[3])
       mtext(paste("Cp_u = ", signif(Cp.u, 3), sep = ""), 
             side = 1, line = 5, adj = 0, font=2, at = at.col[3])
       mtext(paste("Cp_k = ", signif(Cp.k, 3), sep = ""), 
             side = 1, line = 6, adj = 0, font=2, at = at.col[3])
       if (!is.null(target))
          mtext(paste("Cpm  = ", signif(Cpm, 3), sep = ""), 
                side = 1, line = 7, adj = 0, font=2, at = at.col[3])
       #--
       mtext(paste("Exp<LSL ", signif(exp.LSL, 2), "%", sep = ""), 
             side = 1, line = 3, adj = 0, font=2, at = at.col[4])
       mtext(paste("Exp>USL ", signif(exp.USL, 2), "%", sep = ""), 
             side = 1, line = 4, adj = 0, font=2, at = at.col[4])
       mtext(paste("Obs<LSL ", signif(obs.LSL, 2), "%", sep = ""), 
             side = 1, line = 5, adj = 0, font=2, at = at.col[4])
       mtext(paste("Obs>USL ", signif(obs.USL, 2), "%", sep = ""), 
             side = 1, line = 6, adj = 0, font=2, at = at.col[4])
     }

  if (print)
     { cat("\nProcess Capability Analysis\n")
       cat("\nCall:\n", deparse(match.call()), "\n\n", sep = "")

       cat(paste(formatC("Number of obs = ", width=16), 
                 formatC(n, width=12, flag="-"), 
                 formatC("Target = ", width=10), 
                 formatC(target, digits=options()$digits, flag="-"),
                 "\n", sep=""))
       cat(paste(formatC("Center = ", width=16), 
                 formatC(center, digits=options()$digits, width=12, flag="-"),
                 formatC("LSL = ", width=10), 
                 formatC(LSL, digits=options()$digits, flag="-"),
                 "\n", sep=""))
       cat(paste(formatC("StdDev = ", width=16), 
                 formatC(std.dev, digits=options()$digits, width=12, flag="-"),
                 formatC("USL = ", width=10), 
                 formatC(USL, digits=options()$digits, flag="-"),
                 "\n", sep=""))

       cat("\nCapability indices:\n\n")
       print(tab, digits=4, na.print="", print.gap=2)
       cat("\n")
       cat(paste("Exp<LSL ", format(exp.LSL, digits=2), "%   ", 
                 "Obs<LSL ", format(obs.LSL, digits=2), "% \n", sep=""))
       cat(paste("Exp>USL ", format(exp.USL, digits=2), "%   ", 
                 "Obs>USL ", format(obs.USL, digits=2), "% \n", sep=""))
     }
   
  invisible(list(nobs = n, center = center, std.dev = std.dev, 
                 target=target, 
                 spec.limits = { sl <- c(LSL, USL)
                                 names(sl) <- c("LSL", "USL")
                                 sl },
                 indices = tab, 
                 exp = { exp <- c(exp.LSL, exp.USL)/100
                         names(exp) <- c("< LSL", "> USL")
                         exp }, 
                 obs = { obs <- c(obs.LSL, obs.USL)/100
                         names(obs) <- c("< LSL", "> USL")
                         obs }
                 ))
}
               
    
"process.capability.sixpack" <- function(object, spec.limits, target, nsigmas, std.dev)
{
# Process capability sixpack plots (based on an idea in Minitab)

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!object$type=="xbar")
     stop("Sixpack analysis only available for 'qcc' object of type \"xbar\"")
  if (missing(std.dev))
     std.dev <- object$std.dev
  if (missing(spec.limits))
     stop("specification limits must be provided as vector of length two")
  if (!length(spec.limits)==2)
     stop("specification limits must be provided as vector of length two")
  if (missing(target))
     { if (is.null(object$target))
          target <- mean(spec.limits, na.rm=TRUE)
       else
          target <- object$target }
  else
     { if (target < min(spec.limits) | target > max(spec.limits))
          warning("target value is not within specification limits...") }

  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas
  
  bg.margin.old <- .qcc.options$bg.margin
  bg.figure.old <- .qcc.options$bg.figure
  .qcc.options$bg.margin <<- .qcc.options$bg.figure <<- "white"

  oldpar <- par(no.readonly = TRUE)
  on.exit({ .qcc.options$bg.margin <<- bg.margin.old
            .qcc.options$bg.figure <<- bg.figure.old 
            par(oldpar) })

  layout(matrix(c(1,2,3,4,5,6),3,2), widths = c(2,1), heights=c(1,1,1))
  # layout.show()
  par(mar=c(5,4,2,1), cex=.qcc.options$cex)

  # 1)
  Object <- qcc(object$data, type=object$type, center=object$center,
                std.dev=std.dev, target=target, nsigmas=nsigmas, 
                data.name = object$data.name, plot=FALSE)
  plot(Object, add.stats = FALSE, chart.all = TRUE, restore.par = FALSE)
  
  # 2)
  if (min(Object$sizes)>10) type.varchart <- "S" 
  else                      type.varchart <- "R" 
  q <- qcc(Object$data, type=type.varchart, std.dev=std.dev, nsigmas=nsigmas,
           data.name = Object$data.name, plot=FALSE)
  plot(q, add.stats = FALSE, restore.par = FALSE)
  
  # 3)
  runs <- max(1,nrow(Object$data)-25):nrow(Object$data)
  matplot(Object$data[runs,], col=1, pch=1, cex=0.7, 
          ylab=Object$data.name, xlab="Group", main="Run chart")
  abline(h=Object$center, lty=2)

  # 4)
  pcap <- process.capability(Object, target=target, std.dev=std.dev,
                             nsigmas=nsigmas, spec.limits=spec.limits,
                             add.stats=FALSE, print=FALSE, restore.par=FALSE)
  # 5)
  qqnorm(Object$data)
  qqline(Object$data)
  # 6)
  tol.limits <- c(Object$center-nsigmas*std.dev, 
                  Object$center+nsigmas*std.dev)
  par(mar=c(4,1,3,1))
  xlim <- range(tol.limits, spec.limits, target)
  xlim[1] <- xlim[1]-diff(xlim)
  plot(0:7,0:7, type="n", xlim=xlim, xlab="", ylab="", axes=FALSE)
  title("Capability plot")
  axis(1)
  text(mean(spec.limits), 2, "Specification limits", cex=1, font=1)
  arrows(pcap$spec.limits[1], 1, pcap$spec.limits[2], 
         1, length=0.1, angle=90, code=3)
  if (!is.null(target))
     points(target, 1, pch=3, cex=1)
  text(mean(tol.limits), 4, "Process tolerance", cex=1, font=1)
  arrows(tol.limits[1], 3,  tol.limits[2], 3, length=0.1, angle=90, code=3)
  points(Object$center, 3, pch=3, cex=1)
  usr <- par()$usr
  text(usr[1], 6, paste("Center = ", signif(Object$center, options()$digits), sep = ""), pos=4, font=1, cex=1)
  text(usr[1], 5, paste("StdDev = ", signif(std.dev, options()$digits), sep = ""), pos=4, font=1, cex=1)
  if (!is.null(target))
     text(usr[1], 4, paste("Target = ", signif(target, options()$digits), sep = ""), pos=4, font=1, cex=1)
  text(usr[1], 3, paste("Cp   = ", signif(pcap$indices[1,1], 3), sep = ""),
       pos=4, font=1, cex=1)
  text(usr[1], 2, paste("Cp_k = ", signif(pcap$indices[4,1], 3), sep = ""),
       pos=4, font=1, cex=1)
  if (!is.null(target))
     text(usr[1], 1, paste("Cpm  = ", signif(pcap$indices[5,1], 3), sep = ""),
          pos=4, font=1, cex=1)

  return(invisible())
}




#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

"pareto.chart" <- function(x, ylab = "Frequency", xlab, ylim, main, col = heat.colors(length(x)), ...)
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
  oldpar <- par(mar = par("mar")+mar, las = las, no.readonly = TRUE)
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


#-------------------------------------------------------------------#
#                                                                   #
#                  Cause-and-Effect diagram                         #
#                                                                   #
#-------------------------------------------------------------------#

"cause.and.effect" <- function(cause, effect, title = "Cause-and-Effect diagram", cex=c(1,0.9,1), font=c(1,3,2))
{

  # running mean of successive pairs of obs
  mean2 <- function(x)
  { m <- rep(NA, length(x)-1)
    for (i in 1:(length(x)-1))
         m[i] <- mean(x[c(i,i+1)])
    return(m)
  }

  nc <- length(cause)
  ncup <- nc - round(nc/2)
  nclo <- nc - ncup
  ncc <- max(sapply(cause, length))

  oldpar <- par(mar=c(1,1,3,1), no.readonly = TRUE)
  on.exit(par(oldpar))

  plot(0:100, 0:100, type="n", xlab="", ylab="", axes=FALSE, main=title)
  usr <- par("usr")
  we <- strwidth(effect, units="user")*1.1
  wc <- max(unlist(sapply(cause, strwidth, units="user")))
  hc <- max(strheight(effect, units="user"),
            unlist(sapply(cause, strheight, units="user")))

  # draw effect
  arrows(0, 50, usr[2]-we-1, 50, code=2, length=0.1, angle=20)
  text(usr[2]-we,50, effect, adj=c(0,0.5), cex=cex[3], font=font[3])

  # draw branches and cause labels
  a <- (usr[2]-we)/(max(ncup,nclo)+1)
  ac <-  a*(0:(max(ncup,nclo)))
  for (i in 1:(length(ac)-1))
      { segments(mean2(ac)[i], 95, ac[i+1], 50) 
        text(mean2(ac)[i], 96, names(cause)[i], 
             pos=3, offset=0.5, cex=cex[1], font=font[1]) 
        if (i <= nclo)
           { segments(mean2(ac)[i], 5, ac[i+1], 50)
             text(mean2(ac)[i], 4, names(cause)[[ncup+i]], 
                  pos=1, offset=0.5, cex=cex[1], font=font[1]) }
      }

  # draw labels for upper branches
  for (j in 1:ncup)
      { b <- (50-95)/(ac[j+1]-mean2(ac)[j])
        a <- 95-b*mean2(ac)[j]
        y <- rev(50+cumsum((95-50)/(ncc+1))*(1:(ncc)))
        x <- (y-a)/b
        for (i in 1:length(y))
            { label <- cause[[j]][i]
              if (!is.na(label))
                 text(x[i], y[i], label, pos=4, 
                      offset=0.2, cex=cex[2], font=font[2]) 
             }
      }       
  # draw labels for lower branches
  for (j in 1:ncup)
      { b <- (50-5)/(ac[j+1]-mean2(ac)[j])
        a <- 5-b*mean2(ac)[j]
        y <- cumsum((95-50)/(ncc+1))*(1:(ncc))
        x <- (y-a)/b
        if (j <= nclo)
           for (i in 1:length(y))
               { label <- cause[[ncup+j]][i]
                 if (!is.na(label))
                    text(x[i], y[i], label, pos=4, 
                         offset=0.2, cex=cex[2], font=font[2])
               }
      }

  invisible()
}

"qcc.groups" <- function(data, sample)
{
  if(length(data)!=length(sample))
    stop("data and sample must be vectors of equal length")
  x <- lapply(split(data, sample), as.vector)
  lx <- sapply(x, length)
  for(i in which(lx != max(lx)))
      x[[i]] <- c(x[[i]], rep(NA, max(lx)-lx[i]))
  x <- t(sapply(x, as.vector))
  return(x)
}  


"qcc.overdispersion.test" <- function(x, size, 
                            type=ifelse(missing(size), "poisson", "binomial"))
{
  type <- match.arg(type, c("poisson", "binomial"))
  if (type=="binomial" & missing(size))
     stop("binomial data require argument \"size\"")
  if (!missing(size))
     if (length(x) != length(size))   
        stop("arguments \"x\" and \"size\" must be vector of same length")

  n <- length(x)
  obs.var <- var(x)
  if (type=="binomial")
     { p <- sum(x)/sum(size)
       theor.var <- mean(size)*p*(1-p) }
  else if (type=="poisson")
          { theor.var <- mean(x) }
       else
          stop("invalid \"type\" argument. See help.")
     
  D <- (obs.var * (n-1)) / theor.var
  p.value <- 1-pchisq(D, n-1)
  
  out <- matrix(c(obs.var/theor.var, D, signif(p.value,5)), 1, 3)
  rownames(out) <- paste(type, "data")
  colnames(out) <- c("Obs.Var/Theor.Var", "Statistic", "p-value") 
  names(dimnames(out)) <- c(paste("Overdispersion test"), "")
  return(out)
}  

#
# Options retrieval and setting
#
".qcc.options" <- list(exp.R.unscaled = c(NA, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078, 3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588, 3.640, 3.689, 3.735, 3.778, 3.819, 3.858, 3.895, 3.931),
                      se.R.unscaled = c(NA, 0.8525033, 0.8883697, 0.8798108, 0.8640855, 0.8480442, 0.8332108, 0.8198378, 0.8078413, 0.7970584, 0.7873230, 0.7784873, 0.7704257, 0.7630330, 0.7562217, 0.7499188, 0.7440627, 0.7386021, 0.7334929, 0.7286980, 0.7241851, 0.7199267, 0.7158987, 0.7120802, 0.7084528, 0.7050004, 0.7017086, 0.6985648, 0.6955576, 0.6926770, 0.6899137, 0.6872596, 0.6847074, 0.6822502, 0.6798821, 0.6775973, 0.6753910, 0.6732584, 0.6711952, 0.6691976, 0.6672619, 0.6653848, 0.6635632, 0.6617943, 0.6600754, 0.6584041, 0.6567780, 0.6551950, 0.6536532, 0.6521506),
                      beyond.limits = list(pch=19, col="red"),
                      violating.runs = list(pch=19, col="orange"),
                      run.length = 5,
                      bg.margin = "lightgrey",
                      bg.figure = "white",
                      cex = 0.8)

"qcc.options" <- function (...)
{
  if(nargs() == 0) return(.qcc.options)
  current <- .qcc.options
  temp <- list(...)
  if(length(temp) == 1 && is.null(names(temp))) 
    { arg <- temp[[1]]
      switch(mode(arg),
             list = temp <- arg,
             character = return(.qcc.options[[arg]]),
             stop(paste("invalid argument:", arg))) }
  if(length(temp) == 0) return(current)
  name <- names(temp)
  if(is.null(name)) stop("options must be given by name")
  changed <- current[name]
  current[name] <- temp
  if(sys.parent() == 0) 
       env <- pos.to.env( match("package:qcc", search()) )
  else env <- parent.frame()
  assign(".qcc.options", current, envir = env)
  invisible(current)
}
