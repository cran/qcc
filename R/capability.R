#-------------------------------------------------------------------#
#                                                                   #
#    Process Capability Analysis                                    #
#                                                                   #
#-------------------------------------------------------------------#

process.capability <- function(object, spec.limits, target, std.dev, nsigmas, confidence.level = 0.95, breaks="scott", add.stats=TRUE, print=TRUE, restore.par=TRUE)
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
     { target <- mean(spec.limits, na.rm=TRUE) }
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

  oldpar <- par(bg  = qcc.options("bg.margin"), 
                cex = qcc.options("cex"),
                mar = if(add.stats) 
                            c(9+is.null(center)*-1, 2, 4, 2) + 0.1
                      else  par("mar"),
                no.readonly = TRUE)
  if (restore.par) on.exit(par(oldpar))

  plot(0, 0, type="n", xlim = xlim, ylim = ylim,
       axes = FALSE, ylab="", xlab = "", main = title)
  usr <- par()$usr
  rect(usr[1], usr[3], usr[2], usr[4], col = qcc.options("bg.figure"))
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

  if(add.stats) 
    { # computes the x margins of the figure region
      plt <- par()$plt
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.07, 0.35, 0.56, 0.75)
      # write info at bottom
      #--
      mtext(paste("Number of obs = ", n, sep = ""), 
            side = 1, line = 3, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Center = ", signif(center, options()$digits), sep = ""), 
            side = 1, line = 4, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("StdDev = ", signif(std.dev, options()$digits), sep = ""), 
            side = 1, line = 5, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      #--
      if(!is.null(target))
         msg <- paste("Target = ", signif(target, options()$digits), sep = "")
      else
         msg <- paste("Target = ", sep = "")     
      mtext(msg, side = 1, line = 3, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("LSL = ", signif(LSL, options()$digits), sep = ""), 
            side = 1, line = 4, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("USL = ", signif(USL, options()$digits), sep = ""), 
            side = 1, line = 5, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      #--
      mtext(paste("Cp     = ", signif(Cp, 3), sep = ""), 
            side = 1, line = 3, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Cp_l  = ", signif(Cp.l, 3), sep = ""), 
            side = 1, line = 4, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Cp_u = ", signif(Cp.u, 3), sep = ""), 
            side = 1, line = 5, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Cp_k = ", signif(Cp.k, 3), sep = ""), 
            side = 1, line = 6, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      if(!is.null(target))
         mtext(paste("Cpm  = ", signif(Cpm, 3), sep = ""), 
               side = 1, line = 7, adj = 0, at = at.col[3],
               font = qcc.options("font.stats"),
               cex = qcc.options("cex.stats"))
      #--
      mtext(paste("Exp<LSL ", signif(exp.LSL, 2), "%", sep = ""), 
            side = 1, line = 3, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Exp>USL ", signif(exp.USL, 2), "%", sep = ""), 
            side = 1, line = 4, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Obs<LSL ", signif(obs.LSL, 2), "%", sep = ""), 
            side = 1, line = 5, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
      mtext(paste("Obs>USL ", signif(obs.USL, 2), "%", sep = ""), 
            side = 1, line = 6, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = qcc.options("cex.stats"))
    }

  if(print)
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

process.capability.sixpack <- function(object, spec.limits, target, nsigmas, std.dev)
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
     { target <- mean(spec.limits, na.rm=TRUE) }
  else
     { if (target < min(spec.limits) | target > max(spec.limits))
          warning("target value is not within specification limits...") }

  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas

  bg.margin.old <- qcc.options("bg.margin")
  bg.figure.old <- qcc.options("bg.figure")
  qcc.options("bg.margin" = "white")
  qcc.options("bg.figure" = "white")

  oldpar <- par(no.readonly = TRUE)
  on.exit({ qcc.options("bg.margin" = bg.margin.old)
            qcc.options("bg.figure" = bg.figure.old)
            par(oldpar) })

  layout(matrix(c(1,2,3,4,5,6),3,2), widths = c(2,1), heights=c(1,1,1))
  # layout.show()
  par(mar=c(5,4,2,1), cex = qcc.options("cex"))

  # 1)
  Object <- qcc(object$data, type=object$type, center=object$center,
                std.dev=std.dev, nsigmas=nsigmas, 
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

