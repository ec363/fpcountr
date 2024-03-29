% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_kfactor.R
\name{get_kfactor}
\alias{get_kfactor}
\title{Get k-factor of a buffer}
\usage{
get_kfactor(
  buffer_used = "water",
  concentration_used = 0,
  temperature_used = 25
)
}
\arguments{
\item{buffer_used}{string corresponding to buffer. Must match a buffer entry
in dataset kfactors_buffers_data. Default is "water".}

\item{concentration_used}{numeric value of concentration of buffer used.
Default is 0 (pure water), so it needs changing if buffer isn't water. Take
care to ensure the units are correct (buffers specified in M require
concentrations in M not mM). This value is ignored if buffer = "water".}

\item{temperature_used}{numeric value of temperature in oC. Default is 25.}
}
\description{
Calculates k-factor of a given buffer/temperature combination according to
data from Thermo Scientific for different buffers and temperatures. The
k-factor is defined as the A975-A900 for a 1cm pathlength. The calculations
are made under the assumption that the fold changes for buffers from water,
and for temperatures from 25oC may be combined to get an approximate estimate
of the kfactor for a given condition.
}
\examples{
get_kfactor(buffer_used = "TBS", concentration_used = 0.005, temperature_used = 30)
}
