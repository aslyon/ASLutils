#Functions for linear unmixing of FRET data

#' Fit spectra as a linear combination of basis spectra
#'
#' A spectrum recorded from a sample containing more than one fluorophore (or chromophore) is the linear combination of the fluorophores in isolation. Spectra from samples containing unknown contributions from multiple fluorophores can be "unmixed" by least-squares fitting as the sum of multiple basis spectra. This is particularly useful for Forster resonance energy transfer (FRET) experiments. This function performs the least-squares fit and provides several tools useful in FRET applications.
#'
#' @param spectra Matrix with spectral data arranged so that each column contains one spectrum.
#' @param concentrations Numeric vector of concentrations, for instance the concentrations of a titrated component in a binding experiment. The order of concentrations must correspond to the ordering of columns in the \code{spectra} matrix. Required if \code{fitted = FALSE}.
#' @param accConc Numeric. The concentration of the acceptor fluorophore. Required if \code{accCorrect = TRUE} to allow for correction for direct excitation of the acceptor fluorophore in FRET experiments.
#' @param accCorrect Logical. If \code{TRUE}, the contribution of the acceptor signal due to direct excitation (i.e., not FRET) is subtracted before further analysis.
#' @param normalize Logical. If \code{TRUE} FRET is calculated as Acceptor/(Donor + Acceptor). Otherwise FRET is simply Acceptor/Donor.
#' @param fitted Logical. If \code{TRUE}, this function returns a list with components \code{Fitted}, a matrix with fitted spectral values, and \code{Coefs}, containing the coefficients from the least-squares fit.
#' @param average Logical. If \code{TRUE} calculates mean FRET values and standard deviations. Otherwise unaveraged FRET values are returned.
#' @param basis Matrix where \code{nrow(basis) == nrow{spectra}}. Basis spectra for the least-squares fit. The default \code{cfp_yfp_ref} is basis spectra for CFP-YFP FRET experiments. See \code{?cfp_yfp_ref}.
#' @param donorCol Numeric integer. The column of \code{basis} corresponding to the donor fluorophore in FRET experiments.
#' @param accCol Numeric integer. The column of \code{basis} corresponding to the acceptor fluorophore in FRET experiments.
#' @param accCorr Numeric vector with components \code{slope} and \code{intercept}. Values for linear correction term for subtracting acceptor fluorophore signal due to direct excitation of the acceptor (i.e., not FRET). See details.
#' @return A dataframe with FRET values, unless \code{fitted=TRUE} in which a list with components \code{Fitted}, a matrix with fitted spectral values, and \code{Coefs}, containing the coefficients from the least-squares fit.
#' @export
#'
#' @examples
#' #Get the background spectrum
#' background_spectrum <- rowMeans(getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3, startRow=1,
#' endRow=2, startColumn=1, endColumn=3))
#' #Subtract background from experimental spectra
#' TuSC_Spc110_spectra_bgsub <- bgSub(getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3,
#' startRow=3, endRow=9, startColumn=1, endColumn=3), background_spectrum)
#' #Fit the spectra and get FRET values with statistics
#' TuSC_Spc110_binding <- fretFit(TuSC_Spc110_spectra_bgsub, concentrations=sort(rep(c(1500,
#' 750, 375, 187.5, 93.75, 46.875, 0), 3)), accConc=50, average=TRUE)
#' #Fit the spectra and get the fitted spectra
#' TuSC_Spc110_binding <- fretFit(TuSC_Spc110_spectra_bgsub, concentrations=sort(rep(c(1500,
#' 750, 375, 187.5, 93.75, 46.875, 0), 3)), accConc=50, fitted=TRUE)
#' #Fit the spectra and get unaveraged FRET values
#' TuSC_Spc110_binding <- fretFit(TuSC_Spc110_spectra_bgsub, concentrations=sort(rep(c(1500,
#' 750, 375, 187.5, 93.75, 46.875, 0), 3)), accConc=50, average=FALSE)

fretFit <- function(spectra, concentrations, accConc, accCorrect=TRUE, normalize=TRUE, fitted=FALSE, average=TRUE, basis=cfp_yfp_ref, donorCol=1, accCol=2, accCorr = c(slope=0.3986, intercept=-0.7486)) {
	#Add column of 1s for intercept term
	basis_intercept <- cbind(matrix(rep(1, 29), ncol=1), basis)
	dimnames(basis_intercept)[[2]][1] <- "Intercept"
	#Calculate least-squares fit
	fit_coefs <- solve(t(basis_intercept) %*% basis_intercept) %*% (t(basis_intercept) %*% spectra)
	
	#Calculate fitted spectra
	if(fitted == TRUE) {
		fitted_spectra <- basis_intercept %*% fit_coefs
		return(list(Fitted = fitted_spectra, Coefs=fit_coefs))
	}
	
	#Adjust acceptor values to correct for direct excitation of acceptor
	if(accCorrect == TRUE) {
		fit_coefs[accCol+1, ] <- fit_coefs[accCol+1, ] - (accConc*accCorr["slope"] + accCorr["intercept"])
	}
	
	#Calculate FRET
	if(normalize == TRUE) {
		#FRET = YFP / (CFP+YFP)
		FRET <- fit_coefs[accCol+1, ] / apply(fit_coefs[-1, ], 2, sum)
	} else {
		#FRET = YFP/CFP
		FRET <- fit_coefs[accCol+1, ] / fit_coefs[donorCol+1, ]
	}
	
	#Calculate statistics
	stdev <- aggregate(data.frame(Concentrations = concentrations, FRET=FRET), by=list(concentrations), sd)[, 3]
	N <- aggregate(data.frame(Concentrations = concentrations, FRET=FRET), by=list(concentrations), length)[, 3]
	if(average == TRUE) {
		average <- aggregate(data.frame(Concentrations = concentrations, FRET=FRET), by=list(concentrations), mean)
		fret_avg_df <- data.frame(Concentration = average[, 2], FRET=average[, 3], SD=stdev, N=N)
		return(fret_avg_df)
	} else {
		#Make long-format data.frame sorted by increasing concentration
		fret_df_long <- data.frame(Concentration = concentrations, FRET=FRET)[order(concentrations), ]
		uniq_conc <- unique(fret_df_long$Concentration)
		#Assign replicate number to each FRET value - have to do it this way in case replicate count is not equal across all condtions
		replicates <- c()
		for(i in uniq_conc) {
			n <- sum(concentrations == i)
			replicates <- c(replicates, 1:n)
		}
		fret_df_long$Replicate <- replicates
		fret_df_wide <- reshape(fret_df_long, v.names="FRET", idvar="Concentration", timevar="Replicate", direction="wide")
		fret_df_wide$SD <- stdev
		fret_df_wide$N <- N
		rownames(fret_df_wide) <- NULL
		return(fret_df_wide)
	}
}