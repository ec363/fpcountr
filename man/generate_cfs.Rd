% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cfs.R
\name{generate_cfs}
\alias{generate_cfs}
\title{Generate conversion factors for FP calibrations}
\usage{
generate_cfs(
  calibration_csv,
  more_csvs = FALSE,
  more_plots = FALSE,
  sensitivity_plots = FALSE,
  include_only = NULL,
  exclude = NULL,
  gain_fix = FALSE,
  rename_from = NULL,
  rename_to = NULL,
  subset_rows = FALSE,
  rows_to_keep = c("C", "D"),
  separator = "",
  complete_blank = FALSE,
  outfolder = "."
)
}
\arguments{
\item{calibration_csv}{character string. Path of the calibration data csv
file.}

\item{more_csvs}{logical. Optionally save all the intermediate tables in this
function. Defaults to FALSE.}

\item{more_plots}{logical. Optionally save more plots. Defaults to FALSE.}

\item{sensitivity_plots}{logical. Optionally adds a few extra columns to the
data such as the min/max normalised fluorescence and detectable molecules.
It also plots a few extra plots. Defaults to FALSE.}

\item{include_only}{character string. If specified, only includes the
measures (column names) specified here.}

\item{exclude}{character string. If specified, excludes any measures (column
names) specified here.}

\item{gain_fix}{logical. Optionally add "0" before 2-digit gain names ie "GFP
40" -> "GFP 040", or "blueblue 40" -> "blueblue 040"- which fixes ordering
of plots. Defaults to FALSE.}

\item{rename_from}{character string. Rename measures (column names)
containing character string \code{rename_from} to character string specified in
\code{rename_to}. Both \code{rename_from} and \code{rename_to} need to be completed to
trigger renaming.}

\item{rename_to}{character string. Rename measures (column names) containing
character string \code{rename_from} to character string specified in
\code{rename_to}. Both \code{rename_from} and \code{rename_to} need to be completed to
trigger renaming.}

\item{subset_rows}{logical. Should script take a subset of the rows (or whole
table)? Defaults to FALSE.}

\item{rows_to_keep}{character array. If \code{subset_rows} is TRUE, script will
choose rows to keep from this list. Defaults to \code{c("C","D")}.}

\item{separator}{character string that represents the separator between the
channel name and the gain value in the measures columns, eg. for "GFP 40"
it is " ", for "GFP40" it is "" and for "GFP_40" it is "_". Required for
plotting the gain vs conversion factors. Defaults to "".}

\item{complete_blank}{logical. Optionally adds "0" to the concentration_ngul
column for wells identified as blank (protein = "none"). Useful if data
layout is missing these values. Defaults to FALSE.}

\item{outfolder}{character string. Path to folder where output files should
be saved. Defaults to current working directory.}
}
\description{
Generate conversion factors for FP calibrations. Originally based on
\code{flopr::generate_cfs} but with many changes. ... The original function was
intended for fluorescein and microspheres, this one can be used for any
fluorescent calibrant, including proteins. Takes as input a parsed csv of
fluorescence data of a dilution series of FPs at one or more gains, that
contains both data and metadata columns (including: instrument, plate, seal,
channel_name, channel_ex, channel_em, media, calibrant, protein, replicate,
volume, molecular_weight_gmol, concentration_ngul, dilution, rev_dilution,
well, the columns for the fluorescence data, row, column). ... A number of
arguments allow the tweaking of the original data set. Following this, the
data is reshaped, normalised, trimmed of saturated points, summarised and
used to fit a model for the conversion factors from arbitrary to absolute
units. Optional extras: the \code{sensitivity_plots} argument extends this
analysis to identify the limits of detection and relative sensitivity and
dynamic range of each gain. Plots are saved to record the processed data at
every step, allowing for visual sanity checks and troubleshooting. A csv file
with the fitted conversion factors is saved (along with the processed data if
requested with \code{more_csvs}).
}
\examples{
fitvals <- generate_cfs(calibration_csv = "data_parsed.csv", subset_rows = TRUE, rows_to_keep = c("C","D"), outfolder = "cfs_mCherry")
}
