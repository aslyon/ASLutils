#Functions for fitting and plotting binding data

#' Fit a binding model to data
#'
#' Given a dataframe with an independent variable (typically concentration) and observed values (typically FRET or fluorescence intensity) with optional uncertanties, this function fits one of three binding models using the Levenberg-Marquardt non-linear least squares optimizer. The three available models are
#' \itemize{
#'     \item \code{"s"} for simple, single-site saturable binding
#'     \item \code{"t"}  for tight single-site saturable binding where total and unbound concentrations of the titrated
#'           component cannot be assumed equal (i.e. when the non-titrated component concentration is similar to the
#'            dissociation constant)
#'     \item \code{"c"} for cooperative binding
#' }
#'
#' Binding model definitions are as follows, where \eqn{B_{max}} is the theoretical observed value as the independent variable approaches infinity, \eqn{B_{min}} is the observed value when the independent variable is zero, \eqn{K_{d}} is the equilibrium dissociation constant, \eqn{c} is the concentration of the non-titrated component, and \eqn{n} is the Hill coefficient. 
#'
#' Simple:
#' \deqn{y = (B_{max}-B_{min})\frac{x}{x + K_d}+B_{min}}
#'
#' Tight:
#' \deqn{y= (B_{max}-B_{min})\frac{(x + c + K_d) - \sqrt{(x + c + K_d)^2 - 4xc}}{2c}+B_{min}}
#'
#' Cooperative:
#' \deqn{y = (B_{max}-B_{min})\frac{x^n}{x^n + K_d^n}+B_{min}}
#'
#' The \code{kdGuess = "auto"} option finds two data points close to the half-maximal observed value, fits a line to them, then interpolates to find the x-axis value corresponding to the half-maximum point on the binding curve. It is likely to have problems with noisy data, data with outliers, or binding curves where saturation has not been reached.
#'
#' Weighted least-squares is performed with weights \eqn{1/\sigma^2}, so the user should supply uncertainties as standard deviations.
#'
#' @param data Dataframe with at least two columns, typically the output of \code{fretFit}. Columns must contain an independent variable, observed values. An optional column may contain uncertainties for the observed values.
#' @param kdGuess Either \code{"auto"} or numeric. If \code{"auto"}, an initial estimate for the dissociation constant is determined from the data. Alternatively, users may specify an initial estimate as a numeric value.
#' @param probeConcentration Numeric. For \code{model="t"}, the concentration of the non-titrated component in the binding experiment in the same units as the independent variable in \code{data}.
#' @param hillCoefficient Numeric. Initial estimate of the Hill coefficient for \code{model="c"}.
#' @param model One of \code{"s"}, \code{"t"}, or \code{"c"}. The binding model to fit: \code{"s"} for simple, \code{"t"} for tight, \code{"c"} for cooperative.
#' @param weight Logical. Perform weighted non-linear least squares fit using user-provided uncertainties?
#' @param xCol Integer. The column in \code{data} containing the independent variable.
#' @param yCol Integer. The column in \code{data} containing observed values.
#' @param uncCol Integer. The column in \code{data} containing uncertainties for the observed values. Required if \code{weight=TRUE}.
#' @return Object of class \code{"nls"} with optimized binding model parameters.
#' @export
#'
#' @examples
#' fitBinding(data=TuSC_Spc110_binding, model="s")
              
fitBinding <- function(data, kdGuess="auto", probeConcentration, hillCoefficient=2, model=c("s", "t", "c"), weight = TRUE, xCol=1, yCol=2, uncCol=3) {
	#Test if the data comes from fretFit
	if(all(names(data) %in% c("Concentration", "FRET", "SD", "N"))) {
		#Sort data by increasing concentration
		data <- data[order(data$Concentration), ]
		yData <- data$FRET
		xData <- data$Concentration
		if(weight == TRUE) {
			uncData <- data$SD
			if(any(is.na(uncData))) {
				warning("Missing uncertainty for one or more data points. Fitting unweighted data.")
				uncData <- rep(1, length(yData))
			}
		} else {
			uncData <- rep(1, length(yData))
		}
		# dfUse <- data.frame(xData=xData, yData=yData, uncData=uncData)
	} else {
		#Sort data by increasing x-axis value
		data <- data[order(data[, xCol]), ]
		yData <- data[, yCol]
		xData <- data[, xCol]
		if(weight == TRUE) {
			uncData <- data[, uncCol]
			if(any(is.na(uncData))) {
				warning("Missing uncertainty for one or more data points. Fitting unweighted data.")
				uncData <- rep(1, length(yData))
			}
		} else {
			uncData <- rep(1, length(yData))
		}
		# dfUse <- data.frame(xData=xData, yData=yData, uncData=uncData)
	}
	#Get initial parameter estimates
	Bmax <- max(yData)
	Bmin <- min(yData)
	if(kdGuess == "auto") {
		#Get half-maximal y-value
		halfMax <- (diff(range(yData)) / 2) + min(yData)
		#Find the index of the y-value closest to the half-max
		closestIndex <- which(abs(yData - halfMax) == min(abs(yData - halfMax)))
		#Of the two adjacent y-values, find the one closest to the half-max
		nextClosest <- c(closestIndex-1, closestIndex+1)[which(abs(yData[c(closestIndex-1, closestIndex+1)] - halfMax) == min(abs(yData[c(closestIndex-1, closestIndex+1)] - halfMax)))]
		#Perform linear fit at these two points
		obsMatrix <- matrix(yData[c(closestIndex, nextClosest)])
		basis <- matrix(c(1, 1, xData[c(closestIndex, nextClosest)]), ncol=2)
		coefs <- solve(t(basis) %*% basis) %*% (t(basis) %*% obsMatrix)
		#Linear interpolation to get x-axis value corresponding to half-max
		kdInit <- (halfMax - coefs[1]) / coefs[2]
	} else {
		kdInit <- kdGuess
	}
	model <- model[1]
	fit <- switch(model, 
				  s = minpack.lm::nlsLM(yData ~ (Bmax - Bmin) * xData / (xData + Kd) + Bmin, start = list(Bmax=Bmax, Bmin=Bmin, Kd=kdInit), weights=1/uncData^2, model=TRUE),
				  c = minpack.lm::nlsLM(yData ~ (Bmax - Bmin) * xData^n / (xData^n + Kd^n) + Bmin, start = list(Bmax=Bmax, Bmin=Bmin, Kd=kdInit, n=hillCoefficient), weights=1/uncData^2, model=TRUE),
				  t = minpack.lm::nlsLM(yData ~ (Bmax-Bmin) * ((probeConcentration + xData + Kd) - sqrt((probeConcentration + xData + Kd)^2 - 4*probeConcentration*xData))/(2*probeConcentration)+Bmin, start=list(Bmax=Bmax, Bmin=Bmin, Kd=kdInit), weights=1/uncData^2, model=TRUE)
				 )
	return(fit)
}

