#' Basis spectra for linear unmixing of CFP/YFP spectra
#'
#' Fluorescence spectra from gamma-tubulin small complex (gTuSC) containing 
#' Spc97-YFP or Spc98-CFP were recorded on a Molecular Devices SpectraMax
#' M5 plate reader. For gTuSC-CFP, the excitation wavelength was 420 nm with
#' emission recorded through a 455 nm longpass filter in 5 nm increments from 460
#' to 600 nm. Spectra were recorded in the same manner for gTuSC-YFP, but with
#' excitation at 475 nm, a 495 nm longpass filter, and spectra recorded from 495
#' to 600 nm. YFP emission from 460 to 490 nm is set to zero. Multiple spectra
#' were recorded, background subtracted, and averaged. The spectra were then 
#' scaled so that the maximum intensity is one.
#'
#' @format A matrix with 29 rows and 2 columns. The \code{dimnames} attributes
#' indicate CFP or YFP and the wavelengths at which the spectra were recorded 
#' in the columns and rows, respectively.
#'
#' @source Lyon et al. 2016. Molecular Biology of the Cell 27: 2245, figure S1A. \url{http://www.molbiolcell.org/content/27/14/2245.long}
"cfp_yfp_ref"

#' Fluorescence spectra of gamma-TuSC-CFP/YFP recorded at a variety of Spc110 concentrations
#'
#' 50 nm gTuSC-CFP/YFP and varying concentration so Spc110 were mixed to allow
#' assembly of gTuRCs, then fluorescence spectra were recorded on a Molecular Devices SpectraMax
#' M5 platereader. Columns 1-6 (wells A7-B9) are blank spectra recorded from buffer alone. Columns 7-27 are gTuSC+Spc110 spectra. The Spc110 concentrations are \code{sort(rep(c(1500, 750, 375, 187.5,
#' 93.75, 46.875, 0), 3), decreasing=TRUE)}, in nanomolar. These are raw spectra so background 
#' spectra from buffer-only controls
#' must be subtracted before further processing.
#'
#' @format 29 x 27 matrix.
#'
#' @source Lyon et al. 2016. Molecular Biology of the Cell 27: 2245, figure 2E (Dimer WT curve). \url{http://www.molbiolcell.org/content/27/14/2245.long}
"TuSC_Spc110_spectra"

#' Spc110-induced gTuRC assembly curve
#'
#' FRET data for gTuSC-CFP/YFP in the presence of varying concentrations of Spc110.
#'
#' @format 7 x 4 dataframe with columns Concentration, FRET, SD, and N.
#' 
#' @examples
#' #TuSC_Spc110_binding was prepared as follows:
#' background_spectrum <- rowMeans(getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3, startRow=1,
#' endRow=2, startColumn=1, endColumn=3))
#' TuSC_Spc110_spectra_bgsub <- bgSub(getSpectra(TuSC_Spc110_spectra, nRows=9, nCols=3,
#' startRow=3, endRow=9, startColumn=1, endColumn=3), background_spectrum)
#' TuSC_Spc110_binding <- fretFit(TuSC_Spc110_spectra_bgsub, concentrations=sort(rep(c(1500, 750,
#' 375, 187.5, 93.75, 46.875, 0), 3)), accConc=50)
#'
#' @source Lyon et al. 2016. Molecular Biology of the Cell 27: 2245, figure 2E (Dimer WT curve). \url{http://www.molbiolcell.org/content/27/14/2245.long}
"TuSC_Spc110_binding"