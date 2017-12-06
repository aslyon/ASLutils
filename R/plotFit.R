#Plot a binding model

#' Plot a binding model
#'
#' Given a fit model (typically returned by \code{fitBinding}), plot a curve showing the fit evaluated over a range of x-axis values.
#'
#' @param model The model object for which to plot a curve. Typically returned by \code{fitBinding}, though models returned by \code{lm}, \code{nls} are likely to work as well as long as the fitting function was run with \code{model=TRUE}.
#' @param xRange Length-two numeric vector. The x-axis values between which the fit model will be evaluated and a curve plotted.
#' @param xlog Logical. Should the curve be drawn with a logarithmic x-axis scale?
#' @param logBase Numeric. The base of the logarithmic transform used if \code{xlog=TRUE}.
#' @param add Logical. Should the curve be added to an existing graphics device? If \code{FALSE}, the curve is drawn in a new graphics device.
#' @param col Any valid R color specification. The color of the curve.
#' @param lwd Numeric. The line width. \code{lwd=1} is 1/96 of an inch. \code{lwd=0.75} is one point.
#' @param lty Integer. The type of line to draw. 1 is solid, 2 is dashed, 3 is dotted. See \code{lty} under \code{?par} for more details.
#' @param xlab. Character. The x-axis label.
#' @param ylab. Character. The y-axis label.
#' @param xlim. Length-two numeric vector. The x-axis limits of the plot window for \code{add=FALSE}.
#' @param ylim. Length-two numeric vector. The y-axis limits of the plot window for \code{add=FALSE}.
#' @param main. Character. The title of the plot.
#'
#' @return No return value, changes state of graphics device.
#'
#' @export
#'
#' @examples
#' #Fit a model to some data
#' fit <- fitBinding(TuSC_Spc110_binding)
#' #Plot the binding data
#' plotBinding(TuSC_Spc110_binding)
#' #Plot the model curve
#' plotFit(fit)
#
plotFit <- function(model, xRange=NULL, xlog=FALSE, logBase=10, add=TRUE, col=1, lwd=1, lty=1, xlab="Concentration (nM)", ylab="FRET", xlim=NULL, ylim=NULL, main=NULL, ...) {
	#Make sure model was fit using model=TRUE
	if(is.null(model$model)) {
		stop(paste("Model frame not found in \"", deparse(substitute(model)), "\". You should re-run the fitting function with model=TRUE.", sep=""))
	}
	#Get the limits of the x-axis values at which to evaluate model
	if(is.null(xRange) & add==TRUE) {
		#If adding to existing graphics device and user doesn't specify xRange, get from graphics device
		x0 <- par("usr")[1]
		x1 <- par("usr")[2]
	}
	#If plotting in new graphics device and user doesn't specify xRange, get from model object
	if(is.null(xRange) & add==FALSE) {
		if(xlog == TRUE & min(model$model[[2]]) == 0) {
			x0 <- sort(unique(model$model[[2]]))[2]
		} else {
			x0 <- min(model$model[[2]])
		}
		x1 <- max(model$model[[2]])
	}
	#If user specifies xRange, get from that vector
	if(!is.null(xRange)) {
		x0 <- xRange[1]
		x1 <- xRange[2]
	}
	#Get the x-axis values at which to evaluate the model
	#For logarithmic x-axis scaling
	if(xlog == TRUE) {
		#Make sure x0 and x1 are in linear scale
		#If user specifies xRange and adding to existing graphics device, convert to logarithmic
		if(!is.null(xRange) & add==TRUE) {
			x0 <- log(x0, logBase)
			x1 <- log(x1, logBase)
		}
		#If creating new graphics device, x-axis comes from either user-specified xRange in linear units or from model frame in linear units. So if add=F, always convert to log
		if(add==FALSE) {
			x0 <- log(x0, logBase)
			x1 <- log(x1, logBase)
		}
		#To ensure even spacing in log-scale, take exponential of linear sequence from x0 to x1
		xNew <- list(logBase ^ seq(from=x0, to=x1, l=500))
		names(xNew) <- names(model$model)[2]
		#Convert user-specified xlim to logarithmic scale
		if(!is.null(xlim)) {
			xlim <- log(xlim, logBase)
		}
	} else { # For linear x-axis scaling
		#Negative x-axis values are non-sensical, so truncate x0 to 0
		x0 <- ifelse(x0 < 0, 0, x0)
		xNew <- list(seq(from=x0, to=x1, l=500))
		names(xNew) <- names(model$model)[2]
	}
	yNew <- predict(model, xNew)
	if(add == TRUE) {
		if(xlog == TRUE) {
			lines(log(xNew[[1]], 10), yNew, col=col, lwd=lwd, lty=lty, ...)
		} else {
			lines(xNew[[1]], yNew, col=col, lwd=lwd, lty=lty, ...)
		}
		
	} else { #If creating new graphics device
		#Store the current graphic parameters for mar, las, and bty, which will be changed then reset
		mar <- par("mar")
		las <- par("las")
		bty <- par("bty")
		#Set mar, las, and bty
		par(mar=c(4, 4, 1, 0.5), las=1, bty='l')
		if(xlog == TRUE) {
			#Get the positions for x-axis ticks - will be logBase ^ integer, with integer depending on range of independent variable in model
			if(!is.null(xlim)) {
				seqRange <- seq(from=xlim[1], to=xlim[2], l=500)
			} else {
				seqRange <- seq(from=x0, to=x1, l=500)
			}
			at <- unique(floor(seqRange[seqRange != 0]))
			#Get labels for x-axis ticks
			labels <- as.character(round(logBase ^ at, 1))
			#Get minor x-axis tick positions (In log-10, these will be 10, 30, 30, etc., 100, 200, 300, etc...)
			minor.at <- c()
			for(i in 1 : (length(at)-1)) {
				minor.at <- c(minor.at, seq(from = logBase^at[i], to = logBase^at[i+1], by = logBase^at[i])[2:9])
			}
			foo <- round(x1/logBase^max(at), 1) * logBase^max(at)
			last.few <- seq(from = logBase^max(at), to=round(logBase^max(seqRange)/logBase^max(at), 1) * logBase^max(at), by=logBase^max(at))
			minor.at <- c(minor.at, last.few)
			#Setup the the plot and add axes
			plot(log(xNew[[1]], logBase), yNew, col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, main=main, ..., axes=F, type="l")
			axis(1, at=at, labels=labels, tcl=-0.5)
			axis(1, at=log(minor.at, logBase), labels=FALSE, tcl=-0.3)
			axis(2)
			box(bty='l')
		} else { #If linear x-axis scale
			plot(xNew[[1]], yNew, type='l', col=col, lwd=lwd, lty=lty, xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)
		}
		par(mar=mar, las=las, bty=bty)
	}
}