## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)
library(ezmmek)
library(dplyr)

## ---- echo = FALSE, out.width = "60%", fig.cap = "Figure 1: Standards and controls of each protocol."----
knitr::include_graphics(system.file("images", 
                                    "protocols.png", 
                                    package = "ezmmek"))

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("ezmmek")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)

## ---- eval = FALSE------------------------------------------------------------
#  install_github("ccook/ezmmek")

## ---- eval = FALSE------------------------------------------------------------
#  assertable
#  dplyr
#  magrittr
#  nls2
#  purrr
#  rlang
#  tidyr

## ---- eval = FALSE------------------------------------------------------------
#  library(ezmmek)

## ---- include = FALSE---------------------------------------------------------
std_data <- system.file("extdata", 
                        "tyson_std_04172020.csv", 
                        package = "ezmmek")

fsu_data_ibc <- system.file("extdata", 
                            "tyson_sat_german_04172020.csv", 
                            package = "ezmmek")

fsu_data_isc <- system.file("extdata", 
                            "tyson_sat_steen_04172020.csv", 
                            package = "ezmmek")

## -----------------------------------------------------------------------------
std_data <- read.csv(std_data)
dplyr::glimpse(std_data)

## ---- eval = FALSE------------------------------------------------------------
#  site_name: Descriptor column (Tyson Park, Knoxville, TN)
#  std_type: Descriptor column (7-Amino-4-methylcoumarin; AMC)
#  std_conc: Standard concentration ([micro]M)
#  homo_signal: Signal from standard in homogenate (fsu)
#  buffer_signal: Signal from standard in buffer (fsu)

## -----------------------------------------------------------------------------
raw_fsu_data_ibc <- read.csv(fsu_data_ibc)
dplyr::glimpse(raw_fsu_data_ibc)

## ---- eval = FALSE------------------------------------------------------------
#  site_name: Descriptor column (Tyson Park, Knoxville, TN)
#  std_type: Descriptor column (7-Amino-4-methylcourmarin; AMC)
#  time: Timepoints at which signal was collected (hours)
#  signal: Signal from assay (fsu)
#  substrate_conc: Substrate concentration, ([micro]M)
#  replicate: Sample replicate (integer)
#  homo_vol: Homogenate volume (L)
#  soil_mass: Soil mass, or homogenate volume for a water sample (L)
#  std_vol: Standard volume (L)
#  homo_control: Homogenate control (fsu)
#  substrate_control: Substrate control (fsu)

## -----------------------------------------------------------------------------
raw_fsu_data_isc <- read.csv(fsu_data_isc)
dplyr::glimpse(raw_fsu_data_isc)

## ---- eval = FALSE------------------------------------------------------------
#  site_name: Descriptor column (Tyson Park, Knoxville, TN)
#  std_type: Descriptor column (7-Amino-4-methylcourmarin; AMC)
#  time: Timepoints at which signal was collected (hours)
#  signal: Signal from assay (fsu)
#  substrate_conc: Substrate concentration ([micro]Molar)
#  replicate: Sample replicate (integer)
#  kill_control: Signal from autoclaved sample (fsu)

## ---- eval = FALSE------------------------------------------------------------
#  new_ezmmek_sat_fit
#  new_ezmmek_act_calibrate
#  new_ezmmek_act_group
#  new_ezmmek_std_group

## ---- eval = FALSE------------------------------------------------------------
#  new_obj <- new_ezmmek_sat_fit(std.data.fn,
#                                act.data.fn,
#                                ...,
#                                km = NULL,
#                                vmax = NULL,
#                                method = NA)

## ---- eval = FALSE------------------------------------------------------------
#  std.data.fn: Standard curve data file as character string
#  act.data.fn: Raw activity data file as character string
#  ...: User defined column names to join and group std.data.fn and act.data.fn
#  km: Starting value to estimate km. Default value is median of 'sub.conc' values
#  vmax: Starting value to estimate vmax. Default value is max calibrated activity
#  method: Enzyme assay protocol. Must define method as '"ibc"' or '"isc"'

## ---- eval = FALSE------------------------------------------------------------
#  new_obj <- new_ezmmek_act_calibrate(std.data.fn,
#                                      act.data.fn,
#                                      ...,
#                                      method = NA,
#                                      columns = NULL)

## ---- eval = FALSE------------------------------------------------------------
#  std.data.fn: Standard curve data file as character string
#  act.data.fn: Raw activity data file as character string
#  ...: User defined column names to join and group std.data.fn and act.data.fn
#  method: Enzyme assay protocol. Must define method as '"ibc"' or '"isc"'
#  columns: User defined column names to join and group std.data.fn and act.data.fn

## ---- eval = FALSE------------------------------------------------------------
#  new_obj <- new_ezmmek_act_group(act.data.fn,
#                                  ...,
#                                  method = NA,
#                                  columns = NULL)

## ---- eval = FALSE------------------------------------------------------------
#  act.data.fn: Raw activity data file as character string
#  ...: User defined column names to join and group std.data.fn and act.data.fn
#  method: Enzyme assay protocol. Must define method as '"ibc"' or '"isc"'
#  columns: User defined column names to join and group std.data.fn and act.data.fn

## ---- eval = FALSE------------------------------------------------------------
#  new_obj <- new_ezmmek_std_group(std.data.fn,
#                                  ...,
#                                  method = NA,
#                                  columns = NULL)

## ---- eval = FALSE------------------------------------------------------------
#  std.data.fn: Standard curve data file as character string
#  ...: User defined column names to join and group std.data.fn and act.data.fn
#  method: Enzyme assay protocol. Must define method as '"ibc"' or '"isc"'
#  columns: User defined column names to group std.data.fn

## ---- include = FALSE---------------------------------------------------------
std_data <- system.file("extdata", 
                        "tyson_std_04172020.csv", 
                        package = "ezmmek")

fsu_data_ibc <- system.file("extdata", 
                            "tyson_sat_german_04172020.csv", 
                            package = "ezmmek")

fsu_data_isc <- system.file("extdata", 
                            "tyson_sat_steen_04172020.csv", 
                            package = "ezmmek")

## ---- eval = FALSE------------------------------------------------------------
#  sat_fit_ibc <- new_ezmmek_sat_fit(std_data,
#                                    fsu_data_ibc,
#                                    site_name,
#                                    std_type,
#                                    km = NULL,
#                                    vmax = NULL,
#                                    method = "ibc")
#  
#  dplyr::glimpse(sat_fit_ibc)

## -----------------------------------------------------------------------------
sat_fit_isc <- new_ezmmek_sat_fit(std_data, 
                                  fsu_data_isc, 
                                  site_name, 
                                  std_type,
                                  km = NULL,
                                  vmax = NULL,
                                  method = "isc")

dplyr::glimpse(sat_fit_isc)

## -----------------------------------------------------------------------------
act_calibrate_ibc <- new_ezmmek_act_calibrate(std_data,
                                              fsu_data_ibc,
                                              site_name,
                                              std_type,
                                              method = "ibc",
                                              columns = NULL)

dplyr::glimpse(act_calibrate_ibc)

## -----------------------------------------------------------------------------
act_calibrate_isc <- new_ezmmek_act_calibrate(std_data,
                                              fsu_data_isc,
                                              site_name,
                                              std_type,
                                              method = "isc",
                                              columns = NULL)

dplyr::glimpse(act_calibrate_isc)

## -----------------------------------------------------------------------------
act_group_ibc <- new_ezmmek_act_group(fsu_data_ibc,
                                      site_name,
                                      std_type,
                                      method = "ibc",
                                      columns = NULL)

dplyr::glimpse(act_group_ibc)

## -----------------------------------------------------------------------------
act_group_isc <- new_ezmmek_act_group(fsu_data_isc,                                                                   site_name,
                                      std_type,
                                      method = "isc",
                                      columns = NULL)

dplyr::glimpse(act_group_isc)

## -----------------------------------------------------------------------------
std_group_ibc <- new_ezmmek_std_group(std_data,
                                      site_name,
                                      std_type,
                                      method = "ibc",
                                      columns = NULL)

dplyr::glimpse(std_group_ibc)

## -----------------------------------------------------------------------------
std_group_isc <- new_ezmmek_std_group(std_data,
                                      site_name,
                                      std_type,
                                      method = "isc",
                                      columns = NULL)
dplyr::glimpse(std_group_isc)

## ---- fig.width = 6, fig.asp = 0.618, fig.align = "center", warning = FALSE----
plot(sat_fit_isc, site_name, std_type, substrate_type)

## ---- fig.width = 6, fig.asp = 0.618, fig.align = "center", warning = FALSE----
plot(act_calibrate_ibc, site_name, std_type, substrate_type)
plot(act_calibrate_isc, site_name, std_type, substrate_type)

## ---- fig.width = 6, fig.align = "center"-------------------------------------
plot(act_group_ibc, site_name, std_type, substrate_type)
plot(act_group_isc, site_name, substrate_type, std_type, substrate_conc)

## ---- fig.width = 6, fig.asp = 0.618, fig.align = "center"--------------------
plot(std_group_ibc, site_name, std_type)
plot(std_group_isc, site_name, std_type)

