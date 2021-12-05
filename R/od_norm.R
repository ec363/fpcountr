#' Normalise optical density against blank wells
#'
#' Used by `process_plate` function for OD normalisation. Remains virtually
#' unchanged from `flopr::od_norm`.
#'
#' @param pr_data a long data.frame containing your plate reader data
#' @param blank_well the well coordinates of one or more media blanks
#' @param od_name the column name for the optical density data
#'
#' @return an updated data.frame with an additional column "normalised_OD"
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
od_norm <- function(pr_data, blank_well, od_name) {

  pr_data$normalised_OD <- pr_data[, od_name] # copies OD into normalised OD column

  # remove background optical density signal -------------------------------------

  pr_data <- pr_data %>%
    dplyr::group_by(.data$time) %>%
    dplyr::mutate(normalised_OD = .data$normalised_OD - mean(.data$normalised_OD[.data$well %in% blank_well]))
  # for each timepoint, normalises OD to mean of blank wells
  # OD = (copied over OD value) - mean(ODs in blank wells)

  return(as.data.frame(pr_data))
}
