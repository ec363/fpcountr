#' Convert arbitrary fluorescence units to calibrated units
#'
#' Used by `process_plate` function for fluorescence calibration. Function adds
#' calibrated fluorescence column to the data, which is returned. Originally
#' based on `flopr::calibrate_flu`, but with multiple changes. A list of
#' arguments have been added to allow selection of required conversion factor
#' from table that may include conversion factors from multiple instruments,
#' FPs, etc, and function now includes error checks that report to the user if
#' conversion factors are missing.
#'
#' @param pr_data a data.frame of parsed plate reader data
#' @param od_name the column name for the optical density data
#' @param flu_instr instrument name
#' @param flu_channel fluorescent channel name
#' @param flu_gain gain
#' @param flu_slug name of fluorescent protein in FPbase slug format
#' @param conversion_factors_csv path of the csv file
#'   containing predicted conversion factors for the fluorescent channels
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return
calibrate_flu <- function(pr_data,
                          od_name,

                          flu_instr,
                          flu_channel,
                          flu_gain,
                          flu_slug,
                          flu_label,

                          do_quench_correction = do_quench_correction,

                          conversion_factors_csv
) {

  # Get conversion factors ---------------------------------------------------------------------------------------

  conversion_factors <- utils::read.csv(conversion_factors_csv)
  conversion_factors

  # Filter by each relevant factor ---------------------------------------------------------------------------------------

  # instrument # channel # gain # FP
  flu_cfs <- conversion_factors %>%
    dplyr::filter(.data$instrument == flu_instr) %>%
    dplyr::filter(.data$channel_name == flu_channel) %>%
    # dplyr::filter(.data$gain == flu_gain) %>% # do gain later
    dplyr::filter(.data$calibrant == flu_slug)
  flu_cfs

  # Error checking
  if(nrow(flu_cfs) == 0) {

    # if there is no calibration for the fluorescent protein
    # return from this function
    message("Error: No conversion factor found for: ", flu_instr, " - ", flu_channel, " - ", flu_slug)
    return(pr_data)

  }

  # Filter by gain ---------------------------------------------------------------------------------------

  if(flu_gain %in% flu_cfs$gain){

    # if flu_gain is listed in conversion table

    # select it
    flu_cfs <- flu_cfs %>%
      dplyr::filter(.data$gain == flu_gain)
    flu_cfs

    # Error checking
    if(nrow(flu_cfs) > 1){

        message("Warning: More than one conversion factor found for: ", flu_instr, " - ", flu_channel, " - ", flu_gain, " - ", flu_slug)
        message("Taking first.")

        flu_cfs <- flu_cfs[1,]
    }

    # Get cfs from remaining table
    if(is.numeric(flu_cfs$cf)){ # checks column existence and whether numeric

      this_cf <- flu_cfs$cf
      message("Calibrating ", flu_channel, " fluorescence channel with conversion factor ", signif(this_cf, digits = 3), "...")

    } else {

      message("Error: Numeric conversion factor column does not exist.") # throw error if it doesn't exist or isn't numeric
      return(pr_data)

    }

  } else {

    # if flu_gain is not listed in conversion table

    # message("No conversion factor found for gain = ", flu_gain, ".")

    ### Fit cf to gain relation to get cf for specific gain
    message("Fitting conversion factor ~ gain model to estimate cf at gain = ", flu_gain, ".")
    model <- stats::lm(log10(cf) ~ poly(gain, 2), data = flu_cfs)
    model
    this_cf <- 10 ^ stats::predict(model, data.frame(gain = flu_gain))
    this_cf
    message("Calibrating ", flu_channel, " fluorescence channel with conversion factor ", signif(this_cf, digits = 3), "...")

    # plot
    plt_gain <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = flu_cfs$gain,
                                      y = 10^stats::predict(model, flu_cfs))) +
      ggplot2::geom_point(ggplot2::aes(x = flu_cfs$gain,
                                       y = flu_cfs$cf)) +
      ggplot2::geom_vline(xintercept = flu_gain, linetype = 2) +
      ggplot2::geom_hline(yintercept = 10 ^ stats::predict(model, data.frame(gain = flu_gain)),
                          linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(x = flu_gain,
                                       y = 10 ^ stats::predict(model, data.frame(gain = flu_gain))),
                          colour = "red", shape = 1, size = 2) +
      ggplot2::scale_x_continuous("Gain") +
      ggplot2::scale_y_continuous("Conversion factor (rfu/molecules)", trans = "log10") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        panel.grid.minor = ggplot2::element_blank()
      )
    plt_gain
    # plotname = paste0(flu_label, "_gain_plot.pdf")
    # ggplot2::ggsave(filename = plotname,
    #                 plot = plt_gain,
    #                 height = 10, width = 10, units = "cm")

  } # if flu_gain is in conversion table

  ## Add new column called "calibrated_[fluorescence channel]" ---------------------------------------------------

  if(isFALSE(do_quench_correction)){
    # if no quench correction has been performed, the column we need to calibrate is called "normalised"

    # "calibrated_[fluorescence channel]" = normalised_[fluorescence channel] / this_cf

    # pr_data[,paste("calibrated_", flu_name, sep="")] <- (pr_data[,paste("normalised_", flu_name, sep="")] / this_cf)
    # pr_data[,paste("calibrated_", flu_channel, sep="")] <- (pr_data[,paste("normalised_", flu_channel, sep="")] / this_cf)
    pr_data[,paste("calibrated_", flu_label, sep="")] <- (pr_data[,paste("normalised_", flu_channel, sep="")] / this_cf) # better. allows multiple calibrations of same channel, eg. for different FPs.
    # head(pr_data)

  } else {
    # if quench correction has been performed, the column we need to calibrate is called "corrected_normalised_"

    pr_data[,paste("calibrated_", flu_label, sep="")] <- (pr_data[,paste("corrected_normalised_", flu_channel, sep="")] / this_cf)
    # head(pr_data)

  }

  return(pr_data)
}
