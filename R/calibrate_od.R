#' Convert arbitrary optical density units to calibrated units
#'
#' Used by `process_plate` function for optical density calibration. Function
#' adds calibrated_OD column to the data, which is returned. Originally based on
#' `flopr::calibrate_od`, except `instr` argument added to allow selection of
#' conversion factor from table that may include multiple instruments, and
#' includes error checks that report to the user if conversion factors are
#' missing.
#'
#' @param pr_data a dataframe of parsed plate reader data
#' @param od_name the column name for the optical density data
#' @param instr character string to represent instrument. If do_calibrate =
#'   TRUE, used for filtering `od_coeffs_csv` and `fluor_coeffs_csv` files for
#'   conversion factors of the relevant instrument.
#' @param conversion_factors_csv path of the CSV file containing
#'   conversion factors for optical density
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return an updated data.frame with an additional column for calibrated OD
#'
#' @keywords internal
#' @export
calibrate_od <- function(pr_data, od_name, instr, conversion_factors_csv) {

  conversion_factors <- utils::read.csv(conversion_factors_csv) # copy in conversion factors
  conversion_factors

  # Filter by each relevant factor
  # instrument # od_name like OD700
  od_cfs <- conversion_factors %>%
    dplyr::filter(.data$instrument == instr) %>%
    # dplyr::filter(.data$measure == od_name)
    dplyr::filter(tolower(.data$measure) == tolower(od_name)) # case insensitive
  od_cfs

  # Error checking
  if(nrow(od_cfs) == 0) {
    # if there is no calibration for the OD, return from this function

    message("Error: No conversion factor found for: ", instr, " - ", od_name)

    return(pr_data)

  } else if(nrow(od_cfs) > 1){

    message("Warning: More than one conversion factor found for: ", instr, " - ", od_name)
    message("Taking first.")

    od_cfs <- od_cfs[1,]
  }

  # Get conversion factor for OD --------------------------------------------

  # Get cfs from remaining table
  if(is.numeric(od_cfs$cf)){ # checks column existence and whether numeric

    this_cf <- od_cfs$cf
    message("Calibrating ", od_name, " channel with conversion factor ", signif(this_cf, digits = 3), "...")

  } else {

    message("Error: Numeric conversion factor column does not exist.") # throw error if it doesn't exist or isn't numeric
    return(pr_data)

  }

  ### OD conversion factor converts FROM particle number TO od value
  ### (hence tiny: from ~100,000,000 to ~0.1 gives 1e-9)
  ### To go from particle number to od value, MULTIPLY by this_cf
  ### To go from od value to particle number, DIVIDE by this_cf

  pr_data$calibrated_OD <- pr_data$normalised_OD / this_cf
  # head(pr_data[, "normalised_OD"])
  # head(pr_data[, "calibrated_OD"])

  return(pr_data)
}
