#' Normalise optical density against blank wells
#'
#' Used by `process_plate` function for OD normalisation. Remains virtually
#' unchanged from `flopr::od_norm`.
#'
#' @param pr_data a long data.frame containing your plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param od_name the column name for the optical density data
#' @param timecourse logical. Is the data timecourse/kinetic data and does it
#'   include a variable called 'time'?
#'
#' @return an updated data.frame with an additional column "normalised_OD"
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
od_norm <- function(pr_data, blank_well, od_name, timecourse) {

  pr_data$normalised_OD <- pr_data[, od_name] # copies OD into normalised OD column

  # remove background optical density signal -------------------------------------

  if(isFALSE(timecourse)){
    pr_data <- pr_data %>%
      dplyr::mutate(normalised_OD = .data$normalised_OD - mean(.data$normalised_OD[.data$well %in% blank_well]))
  } else if(isTRUE(timecourse)){
    pr_data <- pr_data %>%
      dplyr::group_by(.data$time) %>%
      dplyr::mutate(normalised_OD = .data$normalised_OD - mean(.data$normalised_OD[.data$well %in% blank_well])) %>%
      dplyr::ungroup()
    # for each timepoint, normalises OD to mean of blank wells
    # OD = (copied over OD value) - mean(ODs in blank wells)
  }

  return(as.data.frame(pr_data))
}
