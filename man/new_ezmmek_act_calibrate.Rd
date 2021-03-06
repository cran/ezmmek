% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_ezmmek_act_calibrate.R
\name{new_ezmmek_act_calibrate}
\alias{new_ezmmek_act_calibrate}
\title{new_ezmmek_act_calibrate}
\usage{
new_ezmmek_act_calibrate(
  std.data.fn,
  act.data.fn,
  ...,
  method = NA,
  columns = NULL
)
}
\arguments{
\item{std.data.fn}{Standard data file as character string}

\item{act.data.fn}{Activity data file as character string}

\item{...}{User defined column names to join std.data.fn and act.data.fn}

\item{method}{Enzyme assay protocol. Must define method as 'isc' or 'ibc'}

\item{columns}{Column names carried over from parent functions if parent functions used}
}
\description{
Creates dataframe containing calibrated enzyme activity data of class 'new_ezmmek_calibrate'
}
\examples{
\dontrun{new_obj <- new_ezmmek_act_calibrate("data/tyson_std_04172020.csv",
  "data/tyson_sat_steen_04172020.csv",
  site_name,
  std_type,
  method = "isc",
  columns = NULL)
new_obj <- new_ezmmek_act_calibrate("data/tyson_std_04172020.csv",
  "data/tyson_sat_german_04172020.csv",
  site_name,
  std_type,
  method = "ibc",
  columns = NULL)}
}
