% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_fpconc.R
\name{calc_fpconc}
\alias{calc_fpconc}
\title{Calculate FP concentration in molar units}
\usage{
calc_fpconc(
  data_csv,
  timecourse = TRUE,
  flu_channels,
  flu_labels,
  remove_wells,
  get_rfu_vol = TRUE,
  get_mol_vol = FALSE,
  od_specific_total_volume = NULL,
  odmeasure = NULL,
  odmeasure_conversion = NULL,
  plate_type = 96,
  outfolder = "."
)
}
\arguments{
\item{data_csv}{path to a CSV file containing processed plate reader data}

\item{timecourse}{logical. Is the data timecourse/kinetic data and does it
include a variable called \code{time}?}

\item{flu_channels}{the column names for the NORMALISED fluorescence data}

\item{flu_labels}{the column names for the CALIBRATED fluorescence data}

\item{remove_wells}{list of coordinates of wells to be removed from analysis
(e.g. empty wells)}

\item{get_rfu_vol}{logical. if TRUE, uses \code{normalised_FP} and OD-specific
cell volume (\code{od_specific_total_volume}, specified in ul) to calculate FP
concentration as \code{FP/cellvolume} in relative fluorescence units/litre
(RFU/L).}

\item{get_mol_vol}{logical. if TRUE, uses calibrated_FP and OD-specific cell
volume (\code{od_specific_total_volume}, specified in ul) to calculate FP
concentration as moles/L in Molar (M) units.}

\item{od_specific_total_volume}{numeric. OD600-specific total cellular volume
in ul x OD-1 x cm, i.e. the total cellular volume represented by 1 OD600
unit (in 1 cm path length). Recommended value is 3.6 from Volkmer et al.,
2011.}

\item{odmeasure}{character. Which OD measurement is being used in the data?
Specifically, which measurement is represented by the 'normalised_OD_cm1'
column? e.g. "OD600" or "OD700". Purpose is to record this in the table.}

\item{odmeasure_conversion}{numeric. How to convert the measurement specified
by \code{odmeasure} to OD600? i.e. OD600 = OD used / x. Use '1' for OD600 (no
conversion) and 0.79 for OD700.}

\item{plate_type}{type of plate. numeric, i.e. '96' for 96-well plate.
Defines the rows and columns used for plotting figures. Defaults to '96'.}

\item{outfolder}{path to folder where output files should be saved. Defaults
to current working directory.}
}
\value{
a data.frame with columns for each FP/cell calculation
}
\description{
Takes as input plate reader data processed with \code{process_plate()} and uses
normalised/calibrated values to calculate FP concentration. Adds column(s)
for FP concentration as either: (a) \code{normalisedFP/cellvolume} (RFU/L), (b)
\code{calibratedFP/cellvolume} (moles/L, or M). Plots results and returns a
dataframe.
}
\examples{
\dontrun{
conc_data_mCherry <- calc_fpconc(
  data_csv = "mcherry_parsed_processed.csv",
  flu_channels = c("red1red1"), flu_labels = c("mCherry"),
  remove_wells = c("A11"),
  get_rfu_vol = TRUE, get_mol_vol = TRUE,
  od_specific_total_volume = 3.6, odmeasure = "OD700", odmeasure_conversion = 0.79,
  outfolder = file.path("plots"))
}
}
