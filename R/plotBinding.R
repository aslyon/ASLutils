#Functions for plotting binding data

#' Plot binding data
#'
#' Plot binding data as scatter plot with error bars.
#'
#' @param data Dataframe with at least two columns, typically returned by \code{fretFit}. The two required columns must contain independent variable and one observations. Uncertainties in the observed values may be included as an optional third column.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param ylim Length two numeric vector. The lower and upper boundaries of the vertical axis of the plot.
#' @param xlim Length two numeric vector. The left and right boundaries of the horizontal axis of the plot.
#' @param pch Any valid plot character (numeric or character).
#' @param cex Plot character magnification relative to default size (e.g. \code{cex=2} will make the plot character two times larger than default).
#' @param xlog Logical. Should the x-axis values be logarithmically transformed?
#' @param logBase Numeric. If \code{xlog=TRUE}, the base of the logarithmic transform of the x-axis values.
#' @param col Any valid color specification. The foreground color of the plot characters.
#' @param bg Any valid color specification. The background color of the plot characters.
#' @param errCol Any valid color specification. The color of the error bars.
#' @param errLwd Numeric. The line weight of the error bars. \code{errLwd = 1} is a line 1/96 inch thick, \code{errLwd = 0.75} is a line 1 point thick.
#' @param errFeet Logical. Should horizontal lines ("feet") be drawn on the error bars? Default \code{FALSE} will draw vertical lines only.
#' @param main Character. Title for plot.
#' @param xCol Integer. The column of \code{data} containing the independent variable.
#' @param yCol Integer. The column of \code{data} containing observed values.
#' @param uncCol Integer. Optional. The column of \code{data} containing uncertainties in observed values.
#' @param ... Additional arguments to \code{plot} or \code{points}
#'
#' @return No return value, changes state of graphics device.
#"
#' @export
#'
#' @examples
#' #Plot with linear x-axis scale
#' plotBinding(TuSC_Spc110_binding)
#' #Plot with "feet" on the error bars the same color as the plot symbols
#' plotBinding(TuSC_Spc110_binding, errFeet=TRUE, pch=19, col='red', errCol='red')
#' #Plot with log-10 x-axis scale
#' plotBinding(TuSC_Spc110_binding, xlog=T)
#
plotBinding <- function(data, xlab="Concentration (nM)", ylab="FRET", ylim=NULL, xlim=NULL, pch=21, cex=1, xlog=F, logBase=10, col="black", bg=1, errCol="black", errLwd=1, errFeet=FALSE, main=NULL, xCol=1, yCol=2, uncCol=3, ...) {
	#Store the current graphic parameters for mar and las, which will be changed then reset
	mar <- par()$mar
	las <- par()$las
	#Set mar and las
	par(mar=c(4, 4, 1, 0.5), las=1)
	#Check whether data is output from fretFit and assign data to convenient variables
	if(all(c("Concentration", "FRET", "SD", "N") %in% names(data))) {
		xData <- data$Concentration
		yData <- data$FRET
		uncData <- data$SD
	} else { #Otherwise get the columns specified by user from data
		xData <- data[, xCol]
		yData <- data[, yCol]
		if(ncol(data) > 2) {
			uncData <- data[, uncCol]
		} else {
			#If there are no uncertainties, make them all NA
			uncData <- rep(NA, length(xData))
			warning("No uncertainties supplied. Error bars will not be plotted.")
		}
	}
	#Get ylim values taking into account error bars
	if(is.null(ylim)) {
		minIndex <- which(yData == min(yData))
		maxIndex <- which(yData == max(yData))
		ylim <- c(yData[minIndex] - ifelse(is.na(uncData[minIndex]), 0, uncData[minIndex]), yData[maxIndex] + ifelse(is.na(uncData[minIndex]), 0, uncData[maxIndex]))
	}
	#Deal with log-scale x-axis
	if(xlog == TRUE) {
		#Get the positions for x-axis ticks - will be logBase ^ integer, with integer depending on range(xData)
		at <- unique(floor(log(sort(xData[xData != 0]), logBase)))
		#Get labels for x-axis ticks
		labels <- as.character(round(logBase ^ at, 1))
		#Get minor x-axis tick positions (In log-10, these will be 10, 30, 30, etc., 100, 200, 300, etc...)
		minor.at <- c()
		for(i in 1 : (length(at)-1)) {
			minor.at <- c(minor.at, seq(from = logBase^at[i], to = logBase^at[i+1], by = logBase^at[i])[2:9])
		}
		last.few <- seq(from = logBase^max(at), to=round(max(xData)/logBase^max(at), 1) * logBase^max(at), by=logBase^max(at))
		minor.at <- c(minor.at, last.few)
		#Setup the the plot and add axes
		plot(log(xData, logBase), yData, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, main=main, ..., axes=F, type="n")
		axis(1, at=at, labels=labels, tcl=-0.5)
		axis(1, at=log(minor.at, logBase), labels=FALSE, tcl=-0.3)
		axis(2)
		box(bty='l')
		#Draw errorbars with feet if requested
		segments(x0=log(xData, logBase), y0=yData-uncData, y1=yData+uncData, col=errCol, lwd=errLwd)
		if(errFeet == TRUE) {
			footwidth <- .017 * abs(diff(par()$usr[1:2]))  * cex
			segments(x0=rep(log(xData, logBase)-footwidth*0.5, 2), x1=rep(log(xData, logBase)+footwidth*0.5, 2), y0=c(yData-uncData, yData+uncData), col=errCol, lwd=errLwd)
		}
		#Plot the data
		points(log(xData, logBase), yData, pch=pch, col=col, bg=bg, cex=cex, ...)
	} else { #For linear x-axis scaling
		#Setup the plot
		plot(xData, yData, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, main=main, ..., type="n", bty='l')
		#Draw the errorbars with feet if requested
		segments(x0=xData, y0=yData-uncData, y1=yData+uncData, col=errCol, lwd=errLwd)
		if(errFeet == TRUE) {
			footwidth <- .017 * abs(diff(par()$usr[1:2]))  * cex
			segments(x0=rep(xData-footwidth*0.5, 2), x1=rep(xData+footwidth*0.5, 2), y0=c(yData-uncData, yData+uncData), col=errCol, lwd=errLwd)
		}
		#Plot the data
		points(xData, yData, pch=pch, col=col, bg=bg, cex=cex, ...)
	}
	#Rest par to original settings
	par(mar = mar, las = las)
}

#' Add binding data to an existing plot
#'
#' Given an existing graphics device, adds points and error bars.
#'
#' @param data Dataframe with at least two columns, typically returned by \code{fretFit}. The two required columns must contain independent variable and one observations. Uncertainties in the observed values may be included as an optional third column.
#' @param pch Any valid plot character (numeric or character).
#' @param cex Plot character magnification relative to default size (e.g. \code{cex=2} will make the plot character two times larger than default).
#' @param xlog Logical. Should the x-axis values be logarithmically transformed?
#' @param logBase Numeric. If \code{xlog=TRUE}, the base of the logarithmic transform of the x-axis values.
#' @param col Any valid color specification. The foreground color of the plot characters.
#' @param bg Any valid color specification. The background color of the plot characters.
#' @param errCol Any valid color specification. The color of the error bars.
#' @param errLwd Numeric. The line weight of the error bars. \code{errLwd = 1} is a line 1/96 inch thick, \code{errLwd = 0.75} is a line 1 point thick.
#' @param errFeet Logical. Should horizontal lines ("feet") be drawn on the error bars? Default \code{FALSE} will draw vertical lines only.
#' @param xCol Integer. The column of \code{data} containing the independent variable.
#' @param yCol Integer. The column of \code{data} containing observed values.
#' @param uncCol Integer. Optional. The column of \code{data} containing uncertainties in observed values.
#' @param ... Additional arguments to \code{points}
#'
#' @return No return value, changes state of graphics device.
#"
#' @export
#'
#' @examples
#' #Plot some data
#' plotBinding(TuSC_Spc110_binding, pch=21, bg='red', errFeet=TRUE, errCol='red')
#' #Make a dataframe with the inverse curve
#' (just as an example, this doesn't have any practical purpose related to analyzing these data)
#' inverse <- TuSC_Spc110_binding
#' inverse$FRET <- max(inverse$FRET) - inverse$FRET + min(inverse$FRET)
#' #Add the inverse curve to the plot
#' pointsBinding(inverse, pch=21, bg='blue', errFeet=TRUE, errCol='blue')

pointsBinding <- function(data, pch=21, cex=1, xlog=F, logBase=10, col="black", bg=1, errCol="black", errLwd=1, errFeet=FALSE, xCol=1, yCol=2, uncCol=3, ...) {
	#Check whether data is output from fretFit and assign data to convenient variables
	if(all(c("Concentration", "FRET", "SD", "N") %in% names(data))) {
		xData <- data$Concentration
		yData <- data$FRET
		uncData <- data$SD
	} else { #Otherwise get the columns specified by user from data
		xData <- data[, xCol]
		yData <- data[, yCol]
		if(ncol(data) > 2) {
			uncData <- data[, uncCol]
		} else {
			#If there are no uncertainties, make them all NA
			uncData <- rep(NA, length(xData))
			warning("No uncertainties supplied. Error bars will not be plotted.")
		}
	}
	if(xlog == TRUE) {
		#Draw errorbars with feet if requested
		segments(x0=log(xData, logBase), y0=yData-uncData, y1=yData+uncData, col=errCol, lwd=errLwd)
		if(errFeet == TRUE) {
			footwidth <- .017 * abs(diff(par()$usr[1:2]))  * cex
			segments(x0=rep(log(xData, logBase)-footwidth*0.5, 2), x1=rep(log(xData, logBase)+footwidth*0.5, 2), y0=c(yData-uncData, yData+uncData), col=errCol, lwd=errLwd)
		}
		#Plot the data
		points(log(xData, logBase), yData, pch=pch, col=col, bg=bg, cex=cex, ...)
	} else { #For linear x-axis scaling
		#Draw the errorbars with feet if requested
		segments(x0=xData, y0=yData-uncData, y1=yData+uncData, col=errCol, lwd=errLwd)
		if(errFeet == TRUE) {
			footwidth <- .017 * abs(diff(par()$usr[1:2]))  * cex
			segments(x0=rep(xData-footwidth*0.5, 2), x1=rep(xData+footwidth*0.5, 2), y0=c(yData-uncData, yData+uncData), col=errCol, lwd=errLwd)
		}
		#Plot the data
		points(xData, yData, pch=pch, col=col, bg=bg, cex=cex, ...)
	}
}