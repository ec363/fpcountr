#' Correct normalised fluorescence to compensate for cellular quenching effects
#'
#' Used by `process_plate` function for cell quench correction. The presence of
#' cells, particularly at high densities, can 'quench' fluorescence such that
#' fluorescence readings are lower than expected from the equivalent FP
#' concentration in vitro. This function compensates for this effect using
#' empirical data that relates the fold decrease in fluorescence expected, given
#' a certain cell density.
#'
#' @param pr_data a data.frame of parsed plate reader data
#' @param od_type Which OD-type was used? Required for quench correction.
#'   "OD600" or "OD700".
#' @param flu_channel fluorescent channel name
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return
correct_flu <- function(pr_data,
                        od_type, # OD600 or OD700
                        flu_channel
                        ){

  # Get cell quench correction model ---------------------------------------------------------------------------------------

  # Get quench data
  quench_data <- fpcountr::cell_quench_data
  quench_data

  # Rename OD value so that the model will include the correct column name
  # Which column to use for model?
  if(od_type == "OD600"){
    idx <- which(names(quench_data) == "cells_od600" )
    idx
  } else if(od_type == "OD700"){
    idx <- which(names(quench_data) == "cells_od700" )
    idx
  } else {
    message("Error: od_type needs to be either 'OD600' or 'OD700'.")
    return()
  }
  # Rename
  names(quench_data)[idx] <- "normalised_OD_cm1" # whether you used OD600 or OD700, "normalised_OD_cm1" should be the column title after processing
  names(quench_data)

  # Fit model
  # y = k / (x+b)^2 + a
  quench_model <- stats::nls(norm_fc_value ~ k / ((normalised_OD_cm1+b)^2) + a,
                             start = c(k = 1, a = 0.5, b = 2), data = quench_data)
  quench_model

  # Add predicted quenching into table ---------------------------------------------------------------------------------------

  pr_data <- pr_data %>%
    # dplyr::mutate(flu_quench = param1*(normalised_OD_cm1)^2 + param2*(normalised_OD_cm1) + param3) %>% # handwritten version
    dplyr::mutate(flu_quench = stats::predict(quench_model, .)) %>%
    dplyr::mutate(v1 = .data[[paste0("normalised_", flu_channel)]] / flu_quench)
  # flu_quench represents expected quenching that has already happened. so to correct for it we must divide by this
  # eg if quench = 0.90, that suggests we are seeing 90% of the real value. dividing by 0.9 would help us correct this.
  pr_data[1:5,]

  # Rename last column:
  names(pr_data)[ncol(pr_data)] <- paste0("corrected_normalised_", flu_channel)
  names(pr_data)

  return(pr_data)

}

