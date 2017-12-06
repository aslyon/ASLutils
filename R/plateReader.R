#Functions for manipulating plate reader data

#' Open a file containing spectral data
#'
#' Simple wrapper around \code{read.table} for opening tab-delimited spectral data files.
#'
#' @param file Character string giving absolute or relative path to a tab-delimited data file.
#'
#' @return A matrix where each column is a spectrum.
#' @export 
openSpectraFile <- function(file="") {
	spectra_matrix <- as.matrix(read.table(file=file, header=TRUE, row.names=1))
	return(spectra_matrix)
}

#' Get a set of spectra corresponding to a group of wells from a microplate
#'
#' Microplate experiments often have a set of control samples and a set of experimental 
#' samples (for example) where it might be convenient to separate the spectra into 
#' different objects for processing and analysis. Given an input matrix containing 
#' spectral data from a microplate experiment, this function returns the spectra from 
#' a rectangular block of wells.
#'
#' The \code{spectra} matrix must be arranged so that the groups of adjacent columns correspond to wells within a row of a microplate. That is, for a microplate with rows labeled A, B, C,  etc., and columns labeled 1, 2, 3, 4, etc., the columns of the \code{spectra} matrix correspond to wells A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4, etc., in that order.
#'
#' @param spectra Matrix, typically returned by \code{openSpectraFile}. See details for required features of this matrix.
#' @param nRows Numeric or length 2 character vector. If numberic, the number of rows containing samples in the microplate. If character vector, the letter codes for the first and last rows containing sample.
#' @param nCols Numeric. The number of columns containing samples in the microplate.
#' @param startRow Numeric or character. The first row containing samples of interest. Can be specified as the ordinal number for the row of interest or the letter code.
#' @param endRow Numeric or character. The last row containing samples of interest. Can be specified as the ordinal number for the row of interest or the letter code.
#' @param startColumn Numeric. The first column containing samples of interest.
#' @param endRow Numeric. The last column containing samples of interest.
#'
#' @return A matrix where each column is a spectrum.
#'
#' @examples
#' #The matrix TuSC_110_spectra was recorded in wells A7-I9 of a microplate,
#' with buffer-only controls in wells A7-B9.
#' #Get the spectra for the buffer only controls
#' getSpectra(TuSC_110_spectra, nRows=9, nCols=3, startRow=1, endRow=2, startColumn=1,
#' endColumn=3)
#' #Get the experimental spectra using row numbers
#' getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3, startRow=3, endRow=9, startColumn=1,
#' endColumn=3)
#' #Get the experimental spectra using row letters
#' getSpectra(TuSC_Spc110_spectra, nRows=c("A", "I"), ncols=3, startRow="C", endRow="I",
#' startColumn=1, endColumn=3)
#'
#' @export
getSpectra<-function(spectra, nRows, nCols, startRow, endRow=startRow, startColumn, endColumn=startColumn) 	{
		char_num_df <- data.frame(Character = LETTERS, Number = 1:26)
		if(is.character(nRows)) {
			if(length(nRows) != 2) {
				stop("length(nRows) must be 2 if specifying letter code for rows.")
			}
			nRows <- toupper(nRows)
			row_range <- char_num_df[char_num_df$Character %in% nRows, "Number"]
			nRows <- diff(row_range) + 1
		}
		if(is.character(startRow)) {
			if(!is.character(endRow)) {
				stop("startRow and endRow must both be characters if specifiying letter code for rows.")
			}
			startRow <- char_num_df[char_num_df$Character %in% startRow, "Number"]
			endRow <- char_num_df[char_num_df$Character %in% endRow, "Number"]
		} else {
			if(!is.numeric(endRow)) {
				stop("startRow and endRow must both be numeric if specifying row numbers.")
			}
		}
		m<-matrix(1: (nRows*nCols), ncol=nCols, byrow=T)
		spectra.positions<-m[startRow:endRow,startColumn:endColumn]
		output<-spectra[, sort(as.numeric(spectra.positions))]
		return(output)
}

#' Plot a set of spectra
#'
#' Given a matrix containing spectral data, plot each spectrum on the same set of axes.
#'
#' @param spectra Matrix containing spectral data, typically returned by \code{openSpectraFile} or \code{getSpectra}.
#' @param wl Numeric vector of wavelengths.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param bg Any valid R color specification. Background color of plot.
#' @param ... Further arguments to \code{plot}.
#'
#' @return No return value, changes state of graphics device.
#'
#' @examples
#' blanks <- getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3, startRow=1, endRow=2,
#' startColumn=1, endColumn=3)
#' plotSpectra(blanks)
#'
#' @export
plotSpectra <- function(spectra, wl=seq(460, 600, by=5), xlab="Wavelength (nm)", ylab="Fluorescence intensity (AU)", bg='grey80', ...) {
	arguments <- list(...)
	n <- ncol(spectra)
	plot(-1e6, -1e6, xlab=xlab, ylab=ylab, xlim=range(wl), ylim=range(spectra), ...)
	bounds <- par('usr')
	rect(bounds[1], bounds[3], bounds[2], bounds[4], col=bg, border=NA)
	for(i in 1:n) {
		lines(wl, spectra[, i], col=rainbow(n, 5/6)[i])
	}
}

#' Subtract a background spectrum
#'
#' Given a matrix containing spectral data and a vector containing a background spectrum, subtract the background spectrum from each spectrum in the matrix.
#' @param spectra Matrix containing spectral data, typically returned by \code{openSpectraFile} or \code{getSpectra}.
#' @param bgSpectrum Numeric vector where \code{length(bgSpectrum) == nrow(spectra)}.
#'
#' @return Matrix of same dimensions as \code{spectra} containing background-subtracted spectra.
#'
#' @examples
#' blanks <- getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3, startRow=1, endRow=2,
#' startColumn=1, endColumn=3)
#' spectra <- getSpectra(TuSC_Spc110_Spectra, nRows=9, nCols=3, startRow=3, endRow=9,
#' startColumn=1, endColumn=3)
#' spectra_bgsub <- bgSub(spectra=spectra, bgSpectrum=rowMeans(blanks))
#'
#' @export
bgSub <- function(spectra, bgSpectrum) {
	bg.matrix <- matrix(rep(bgSpectrum, ncol(spectra)), ncol=ncol(spectra))
	return(spectra - bg.matrix)
}