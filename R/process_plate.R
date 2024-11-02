#' Process timecourse plate reader data (normalise, correct and calibrate)
#'
#' Normalise and calibrate plate reader measurements to obtain molecular units
#' from the raw data. Originally based on `flopr::process_plate` but with many
#' changes. Takes as input timecourse plate reader data containing optical
#' density measurements and fluorescence measurements as CSV file. First, the
#' function normalises optical density readings and fluorescence readings to an
#' autofluorescence model based on the specified negative wells (see 2020
#' Fedorec et al ACS SynBio for details). The identity of the blank wells need
#' to be specified in `blank_wells`, the reading to be used for optical density
#' in `od_name`, and the fluorescence channels in `flu_channels`. The
#' autofluorescence model can be set with `af_model` on the `neg_well` wells,
#' alternatively, fluorescence can be normalised to blank wells by setting
#' `af_model` to `NULL`. Second, if `do_quench_correction` is TRUE, it
#' compensates for cellular quenching of fluorescence measurements according to
#' the cell density. This requires the specification of `od_type` as "OD600" or
#' "OD700". Third, if `do_calibrate` is TRUE, normalised values are calibrated,
#' where conversion factors are provided as `od_coeffs_csv` for OD and
#' `fluor_coeffs_csv` for fluorescence, and each calibration is specified by the
#' `flu_slugs` to represent the FP used, `flu_gains` to provide the gain used,
#' and `flu_labels` to specify how the relevant plots should be labelled, (e.g.
#' `flu_slugs = c("mcherry", "mtagbfp2")`, `flu_gains = c(60,80)`, `flu_labels =
#' c("RFP, BFP")`).
#'
#' @param data_csv path to a CSV file containing parsed plate reader data
#' @param blank_well the well coordinates of one or more media blanks. Defaults
#'   to "A1".
#' @param timecourse logical. Is the data timecourse/kinetic data and does it
#'   include a variable called 'time'?
#' @param od_name the column name for the optical density data. Defaults to
#'   "OD". If no OD measurements were taken, use `NULL`.
#' @param flu_channels the column names for the fluorescence data. Defaults to
#'   "green1green2".
#' @param flu_channels_rename if specified, what to change the flu_channels
#'   column names to, before computing normalisations and calibrations. Can be
#'   useful if the channel_name in the calibration file is named differently
#'   from the columns in the data file. Needs to be same length as flu_channels,
#'   if not all require changing, specify them anyway to allow positional
#'   replacement (first element in flu_channels_rename replaces first in
#'   flu_channels, etc.). Defaults to NULL.
#' @param af_model model used to fit negative control autofluorescence. For now
#'   these include "polynomial", "inverse_poly", "exponential", "spline" and
#'   "loess". If set to NULL, no model is used, and fluorescence is normalised
#'   akin to OD: by subtracting the value for the blanks. Defaults to "spline".
#' @param neg_well the well coordinates of a non-fluorescent control. Defaults
#'   to "A2".
#' @param do_quench_correction logical. Should function correct for anticipated
#'   quenching of fluorescence, depending on the cell density?
#' @param od_type Which OD-type was used? Required for quench correction.
#'   "OD600" or "OD700".
#' @param do_calibrate logical. Should function convert OD and fluorescence data
#'   to calibrated units? Defaults to FALSE.
#' @param instr character string to represent instrument. If do_calibrate =
#'   TRUE, used for filtering `od_coeffs_csv` and `fluor_coeffs_csv` files for
#'   conversion factors of the relevant instrument.
#' @param flu_slugs character array representing fluorescent proteins (format =
#'   FPbase slug). If do_calibrate = TRUE, used for filtering `fluor_coeffs_csv`
#'   for conversion factors of the relevant FP.
#' @param flu_gains numeric array representing gains of each fluorescent
#'   channel. If do_calibrate = TRUE, used for filtering `fluor_coeffs_csv` for
#'   conversion factors of the relevant gain.
#' @param flu_labels If do_calibrate = TRUE, the column names to be given to
#'   each calibration. May be identical to flu_slug or flu_channel, but
#'   recommended is to make it obvious which FP is being calibrated, e.g.
#'   "mCherry", as channel names may be non-specifically named e.g. "red1red1".
#'   Needs to be same length as flu_slugs and flu_gains.
#' @param od_coeffs_csv if do_calibrate = TRUE, path of the CSV file containing
#'   conversion factors for optical density
#' @param fluor_coeffs_csv if do_calibrate = TRUE, path of the CSV file
#'   containing predicted conversion factors for the fluorescent channels
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @return a data.frame with columns for raw plate reader data, normalised data
#'   and, if do_calibrate = TRUE, calibrated OD and FP data
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   processed_data <- process_plate(
#'     data_csv = "mcherry_parsed.csv",
#'     blank_well = c("A11"),
#'     od_name = "OD600", flu_channels = c("red1"), flu_channels_rename = c("red1red1"),
#'     af_model = NULL,
#'     do_quench_correction = TRUE, od_type = "OD700",
#'     do_calibrate = TRUE, instr = "spark1",
#'     flu_slugs = c("mcherry"), flu_gains = c(80), flu_labels = c("mcherry"),
#'     od_coeffs_csv = "od_coeffs.csv", fluor_coeffs_csv = "flu_coeffs.csv",
#'     outfolder = file.path("data_processed")
#'   )
#' }
process_plate <- function(
    data_csv, blank_well = "A1",

    # timecourse
    timecourse = TRUE,

    # od
    od_name = "OD",

    # fluorescence
    flu_channels = c("green1green2"), # column names in data
    flu_channels_rename = NULL, # if not NULL, for every flu_channels specify new name
    # rename data columns to make sure they match conversion factor table entries

    # autofluorescence model
    af_model = "spline", # options: NULL, "spline", "loess"
    neg_well = "A2", # if af_model = NULL this can be hashed

    # cell quench correction
    do_quench_correction = FALSE,
    od_type, # "OD600" or "OD700"

    # calibrations
    do_calibrate = FALSE,
    instr, # instrument name for the calibration(s)
    # flu_channels above specifies channel(s) for the calibration(s)
    flu_slugs = c(), # FP name in FPbase's slug format for the calibration(s)
    # OR whatever matches the format used in the conversion factor tables
    flu_gains = c(), # gain for calibration(s)
    # each of the above need to be specified for every flu_channels being calibrated
    flu_labels = c(), # how to display the FP/channel combination in plots and file names
    # important if slug is unhelpful or multiple gains used for same FP

    # conversion factors
    od_coeffs_csv, # microsphere conversion factors for od600, od700
    fluor_coeffs_csv, # FP conversion factors

    outfolder = ".") {

  # Get parsed data --------------------------------------------------------------------------------------------------

  pr_data <- utils::read.csv(data_csv, check.names = FALSE)
  pr_data

  # Rename flu channels if requested ---------------------------------------------------------------------------------

  if(!is.null(flu_channels_rename)){

    ## Check there are replacements for each flu channel
    if(!isTRUE(all.equal(length(flu_channels), length(flu_channels_rename), check.attributes = FALSE))){
      message("Error: number of flu_channels must match number of flu_channels_rename.
              To remove requirement for renaming, remove flu_channels_rename or set it to NULL.")
      return()
    }

    ## Rename columns
    for (flu_name_idx in seq_len(length(flu_channels))){
      names(pr_data)
      col_idx <- which(names(pr_data) == flu_channels[flu_name_idx])
      names(pr_data)[col_idx] <- flu_channels_rename[flu_name_idx]
      names(pr_data)
    }

    # Rest of script thinks flu_channels are the final column names needed.
    # So flu_channels needs to be overwritten with flu_channels_rename to prevent downstream errors
    flu_channels <- flu_channels_rename
  }

  # Location for saved plots -----------------------------------------------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # OD processing sections --------------------------------------------------------------------

  if(is.null(od_name)){

    # skip OD processing...

    # rename data
    od_norm_pr_data <- pr_data # otherwise done at od_norm() step # NB. this won't have normalised_OD or normalised_OD_cm1 columns

  } else {

    # proceed with OD processing...

    # Normalise OD (with external od_norm function) --------------------------------------------------------------------

    od_norm_pr_data <- fpcountr::od_norm(pr_data, blank_well, od_name, timecourse = timecourse)
    # adds one extra column - normalised_OD. for each timepoint, normalises OD to mean of blank wells

    # Plot OD data in plate format

    # plot
    if(isFALSE(timecourse)){

      # heatmap1 - raw OD
      max_value <- max(od_norm_pr_data[[od_name]], na.rm = TRUE)
      plt_od <- ggplot2::ggplot(data = od_norm_pr_data,
                                # ggplot2::aes(x = row, y = column, fill = .data[[od_name]])) +
                                ggplot2::aes(x = .data$column, y = row, fill = .data[[od_name]])) +

        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete("", position = "top", # put axis labels at top
                                  limits = factor(unique(od_norm_pr_data$column))) + # put first element at top # factor to make it discrete rather than continuous
        ggplot2::scale_y_discrete("", limits = rev(unique(od_norm_pr_data$row))) + # put first element at top
        viridis::scale_fill_viridis("raw OD",
                                    discrete = FALSE, limits = c(0, max_value),
                                    alpha = 0.4,
                                    na.value = "white") +

        ggplot2::geom_text(ggplot2::aes(label = round(.data[[od_name]], 2)), na.rm = TRUE, size = 2.5) +

        # ggplot2::coord_fixed(ratio = 1) +
        ggplot2::theme_bw() + # base_size = 8
        ggplot2::theme(
          aspect.ratio = 8/12,
          panel.grid = ggplot2::element_blank()
        )
      plt_od
      # plotname <- gsub(".csv", "_OD1.pdf", basename(data_csv))
      plotname <- gsub(".csv", "_OD_1a_raw.pdf", basename(data_csv)) # plot numbering: no-timecourse plot1a
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_od,
                      height = 16, width = 24, units = "cm")

      # normalised OD
      max_value <- max(od_norm_pr_data$normalised_OD, na.rm = TRUE)
      plt_od <- ggplot2::ggplot(data = od_norm_pr_data,
                                ggplot2::aes(x = .data$column, y = row, fill = .data$normalised_OD)) +

        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(od_norm_pr_data$column))) +
        ggplot2::scale_y_discrete("", limits = rev(unique(od_norm_pr_data$row))) +
        viridis::scale_fill_viridis("normalised OD",
                                    discrete = FALSE, limits = c(0, max_value),
                                    alpha = 0.4,
                                    na.value = "white") +

        ggplot2::geom_text(ggplot2::aes(label = round(.data$normalised_OD, 2)), na.rm = TRUE, size = 2.5) +

        ggplot2::theme_bw() + # base_size = 8
        ggplot2::theme(
          aspect.ratio = 8/12,
          panel.grid = ggplot2::element_blank()
        )
      plt_od
      # plotname <- gsub(".csv", "_OD2.pdf", basename(data_csv))
      plotname <- gsub(".csv", "_OD_1b_normalised.pdf", basename(data_csv)) # plot numbering: no-timecourse plot1b
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_od,
                      height = 16, width = 24, units = "cm")

    } else if(isTRUE(timecourse)){

      # plot with caption
      plt_od <- ggplot2::ggplot(od_norm_pr_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[od_name]]),
                           colour = "black", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$normalised_OD),
                           colour = "red", linewidth = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::labs(caption = "black: raw, red: normalised") +
        ggplot2::scale_colour_discrete("") +
        ggplot2::facet_grid(row~column) +
        ggplot2::theme_bw(base_size = 8) +
        ggplot2::theme(
          aspect.ratio = 1,
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
          panel.grid.minor = ggplot2::element_blank()
        )
      plt_od
      # plotname <- gsub(".csv", "_OD.pdf", basename(data_csv))
      plotname <- gsub(".csv", "_OD_1_raw_normalised.pdf", basename(data_csv)) # plot numbering: timecourse plot1
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_od,
                      height = 16, width = 24, units = "cm")
    }

    # Convert to OD (cm-1) --------------------------------------------------------------------

    if("volume" %in% names(od_norm_pr_data)){ # only run if data contains the column 'volume'

      od_norm_pr_data <- od_norm_pr_data %>%
        dplyr::mutate(pathlength = fpcountr::get_pathlength(.data$volume)) %>% # get (estimated) pathlength from volume
        dplyr::mutate(normalised_OD_cm1 = .data$normalised_OD/.data$pathlength) # get OD (cm-1)
      od_norm_pr_data[1,]

      # Plot in plate format
      if(isFALSE(timecourse)){

        # heatmap1 - raw OD
        max_value <- max(od_norm_pr_data$normalised_OD_cm1, na.rm = TRUE)
        plt_od <- ggplot2::ggplot(data = od_norm_pr_data,
                                  ggplot2::aes(x = .data$column, y = row, fill = .data$normalised_OD_cm1)) +

          ggplot2::geom_tile() +
          ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(od_norm_pr_data$column))) +
          ggplot2::scale_y_discrete("", limits = rev(unique(od_norm_pr_data$row))) +
          viridis::scale_fill_viridis(paste0("normalised ", od_name, " (cm-1)"),
                                      discrete = FALSE, limits = c(0, max_value),
                                      alpha = 0.4,
                                      na.value = "white") +

          ggplot2::geom_text(ggplot2::aes(label = round(.data$normalised_OD_cm1, 2)), na.rm = TRUE, size = 2.5) +

          ggplot2::theme_bw() + # base_size = 8
          ggplot2::theme(
            aspect.ratio = 8/12,
            panel.grid = ggplot2::element_blank()
          )
        plt_od

      } else if(isTRUE(timecourse)){

        plt_od <- ggplot2::ggplot(od_norm_pr_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$normalised_OD_cm1,
                                          colour = "normalised_cm1"), linewidth = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0("normalised ", od_name, " (cm-1)")) + # "normalised OD (cm-1)"
          ggplot2::labs(caption = "") +
          ggplot2::scale_colour_discrete("") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8) +
          ggplot2::theme(
            aspect.ratio = 1,
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
          )
        plt_od

      }
      # plotname <- gsub(".csv", "_normODcm1.pdf", basename(data_csv))
      plotname <- gsub(".csv", "_OD_2_pathlength-normalised.pdf", basename(data_csv)) # plot numbering: both plot2
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_od,
                      height = 16, width = 24, units = "cm")

    }

  } # od processing sections

  # Normalise fluorescence data --------------------------------------------------------------------------------------------

  flu_norm_pr_data <- od_norm_pr_data

  ## for each fluorescence channel...
  for (flu_idx in seq_len(length(flu_channels))) {

    ## for each fluorescent protein, run fluorescent normalisation function..

    # Force od_name = NULL data to be af_model = NULL (if by accident it isn't)
    if(is.null(od_name) & !is.null(af_model)){
      message("Autofluorescence model function cannot be run without OD data. Setting af_model to NULL and normalising fluorescence to blank_wells...")
      af_model <- NULL
    }
    # Force timecourse = FALSE data to be af_model = NULL (if by accident it isn't)
    if(isFALSE(timecourse) & !is.null(af_model)){
      message("Autofluorescence model function cannot be run without timecourse data. Setting af_model to NULL and normalising fluorescence to blank_wells...")
      af_model <- NULL
    }

    flu_norm_pr_data <- fpcountr::flu_norm(pr_data = flu_norm_pr_data, neg_well = neg_well, blank_well = blank_well,
                                           flu_name = flu_channels[flu_idx], af_model = af_model, data_csv = data_csv,
                                           timecourse = timecourse, outfolder = outfolder)
    # adds one column per FP of (eg) normalised_[GFP] to the table

    # Plot in plate format

    if(isFALSE(timecourse)){

      # heatmap1 - raw fluor
      max_value <- max(flu_norm_pr_data[[flu_channels[flu_idx]]], na.rm = TRUE)
      plt_flu <- ggplot2::ggplot(data = flu_norm_pr_data,
                                 ggplot2::aes(x = .data$column, y = row, fill = .data[[flu_channels[flu_idx]]])) +

        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(flu_norm_pr_data$column))) +
        ggplot2::scale_y_discrete("", limits = rev(unique(flu_norm_pr_data$row))) +
        viridis::scale_fill_viridis(paste0(flu_labels[flu_idx], " (rfu)"),
                                    discrete = FALSE, limits = c(0, max_value),
                                    alpha = 0.4,
                                    na.value = "white") +

        ggplot2::geom_text(ggplot2::aes(label = round(.data[[flu_channels[flu_idx]]], 2)), na.rm = TRUE, size = 2.5) +

        ggplot2::theme_bw() + # base_size = 8
        ggplot2::theme(
          aspect.ratio = 8/12,
          panel.grid = ggplot2::element_blank()
        )
      plt_flu
      # plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "1.pdf", sep = ""), basename(data_csv))
      plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_1a_raw.pdf", sep = ""), basename(data_csv)) # plot numbering: no-timecourse plot3a
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")

      # heatmap2 - normalised fluor
      max_value <- max(flu_norm_pr_data[[paste0("normalised_", flu_channels[flu_idx])]], na.rm = TRUE)
      plt_flu <- ggplot2::ggplot(data = flu_norm_pr_data,
                                 ggplot2::aes(x = .data$column, y = row, fill = .data[[paste("normalised_", flu_channels[flu_idx], sep = "")]])) +

        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete("", position = "top", # put axis labels at top
                                  limits = factor(unique(flu_norm_pr_data$column))) + # put first element at top # factor to make it discrete rather than continuous
        ggplot2::scale_y_discrete("", limits = rev(unique(flu_norm_pr_data$row))) + # put first element at top
        viridis::scale_fill_viridis(paste0("normalised ", flu_labels[flu_idx], " (rfu)"),
                                    discrete = FALSE, limits = c(0, max_value),
                                    alpha = 0.4,
                                    na.value = "white") +

        ggplot2::geom_text(ggplot2::aes(label = round(.data[[paste("normalised_", flu_channels[flu_idx], sep = "")]], 2)),
                           na.rm = TRUE, size = 2.5) +

        # ggplot2::coord_fixed(ratio = 1) +
        ggplot2::theme_bw() + # base_size = 8
        ggplot2::theme(
          aspect.ratio = 8/12,
          panel.grid = ggplot2::element_blank()
        )
      plt_flu
      # plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "2.pdf", sep = ""), basename(data_csv))
      plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_1b_normalised.pdf", sep = ""), basename(data_csv)) # plot numbering: no-timecourse plot3b
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")

    } else if(isTRUE(timecourse)){

      plt_flu <- ggplot2::ggplot(flu_norm_pr_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[flu_channels[flu_idx]]]),
                           colour = "black", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                        y = .data[[paste("normalised_", flu_channels[flu_idx], sep = "")]]),
                           colour = "red", linewidth = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], " (rfu)")) + # FP name
        ggplot2::labs(caption = "black: raw, red: normalised") +
        ggplot2::scale_colour_discrete("") +
        ggplot2::facet_grid(row~column) +
        ggplot2::theme_bw(base_size = 8) +
        ggplot2::theme(
          aspect.ratio = 1,
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
          panel.grid.minor = ggplot2::element_blank()
        )
      plt_flu
      # plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], ".pdf", sep = ""), basename(data_csv))
      plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_1_raw_normalised.pdf", sep = ""), basename(data_csv)) # plot numbering: timecourse plot3
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")

    }

  }

  out_data <- flu_norm_pr_data # copy normalised data to out_data df that will be returned if do_calibrate = FALSE

  # Correct for cell quenching ----------------------------------------------------------------------

  # Force od_name = NULL data to be do_quench_correction = FALSE (if by accident it isn't)
  if(is.null(od_name) & isTRUE(do_quench_correction)){
    message("Quench correction cannot be run without OD data. Setting do_quench_correction to FALSE and skipping this step...")
    do_quench_correction <- FALSE
  }

  if (do_quench_correction) {

    ### Correct for cell quenching

    ## for each fluorescence channel...
    for (flu_idx in seq_len(length(flu_channels))) {

      out_data <- fpcountr::correct_flu(pr_data = out_data,
                                        od_type = od_type,
                                        flu_channel = flu_channels[flu_idx])

      # plot calibrated fluorescence data
      if(isFALSE(timecourse)){

        # heatmap - corrected
        max_value <- max(out_data[[paste0("corrected_normalised_", flu_channels[flu_idx])]], na.rm = TRUE)
        plt_flu <- ggplot2::ggplot(data = out_data,
                                   ggplot2::aes(x = .data$column, y = row, fill = .data[[paste("corrected_normalised_", flu_channels[flu_idx], sep = "")]])) +

          ggplot2::geom_tile() +
          ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(out_data$column))) +
          ggplot2::scale_y_discrete("", limits = rev(unique(out_data$row))) +
          viridis::scale_fill_viridis(paste0("corrected ", flu_labels[flu_idx], " (rfu)"),
                                      discrete = FALSE, limits = c(0, max_value),
                                      alpha = 0.4,
                                      na.value = "white") +

          ggplot2::geom_text(ggplot2::aes(label = round(.data[[paste("corrected_normalised_", flu_channels[flu_idx], sep = "")]], 2)),
                             na.rm = TRUE, size = 2.5) +

          ggplot2::theme_bw() + # base_size = 8
          ggplot2::theme(
            aspect.ratio = 8/12,
            panel.grid = ggplot2::element_blank()
          )
        plt_flu

      } else if(isTRUE(timecourse)){

        plt_flu <- ggplot2::ggplot(out_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                          y = .data[[paste("normalised_", flu_channels[flu_idx], sep = "")]],
                                          colour = "normalised"),
                             colour = "black", linewidth = 0.5) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                          y = .data[[paste("corrected_normalised_", flu_channels[flu_idx], sep = "")]],
                                          colour = "corrected"),
                             colour = "red", linewidth = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0("corrected ", flu_labels[flu_idx], " (rfu)")) + # FP name
          ggplot2::labs(caption = "black: normalised, red: corrected") +
          ggplot2::scale_colour_discrete("") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8) +
          ggplot2::theme(
            aspect.ratio = 1,
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
          )
        plt_flu

      }
      # plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_corrected.pdf", sep = ""), basename(data_csv))
      plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_2_quench-corrected.pdf", sep = ""), basename(data_csv)) # both plot4
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")

    } # for each fluorescence channel

  } # if do correct

  # Calibrate data ----------------------------------------------------------------------

  if (do_calibrate) {

    if(is.null(od_name)){

      # skip OD calibration...

      # rename data done at calibrate_od() step - n/a
      # NB. data here won't have calibrated_OD column

    } else {

      # proceed w OD calibration...

      ### Calibrate OD to Particle Number ----------------------------------------------------

      out_data <- fpcountr::calibrate_od(pr_data = out_data,
                                         od_name = od_name,
                                         instr = instr,
                                         conversion_factors_csv = od_coeffs_csv)
      utils::head(out_data)
      # adds calibrated_OD column

      # plot
      if(isFALSE(timecourse)){

        # heatmap - calibrated OD
        max_value <- max(out_data$calibrated_OD, na.rm = TRUE)
        plt_od_calib <- ggplot2::ggplot(data = out_data,
                                        ggplot2::aes(x = .data$column, y = row, fill = .data$calibrated_OD)) +

          ggplot2::geom_tile() +
          ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(out_data$column))) +
          ggplot2::scale_y_discrete("", limits = rev(unique(out_data$row))) +
          viridis::scale_fill_viridis(paste0("calibrated ", od_name, " (particles)"),
                                      discrete = FALSE, limits = c(0, max_value),
                                      alpha = 0.4,
                                      na.value = "white") +

          ggplot2::geom_text(ggplot2::aes(label = scales::scientific(.data$calibrated_OD)), na.rm = TRUE, size = 2.5) +

          ggplot2::theme_bw() + # base_size = 8
          ggplot2::theme(
            aspect.ratio = 8/12,
            panel.grid = ggplot2::element_blank()
          )
        plt_od_calib

      } else if(isTRUE(timecourse)){

        plt_od_calib <- ggplot2::ggplot(out_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$calibrated_OD,
                                          colour = "calibrated"), linewidth = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0("calibrated ", od_name, " (particles)")) + # OD700 (particles)
          ggplot2::labs(caption = "") +
          ggplot2::scale_colour_discrete("") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8) +
          ggplot2::theme(
            aspect.ratio = 1,
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
          )
        plt_od_calib

      }
      # plotname <- gsub(".csv", "_ODcalib.pdf", basename(data_csv))
      plotname <- gsub(".csv", "_OD_3_calibrated.pdf", basename(data_csv)) # plot numbering: both plot5
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_od_calib,
                      height = 16, width = 24, units = "cm")

    } # od calibration

    ### Calibrate RFU to Molecule number ----------------------------------------------------

    ## for each fluorescence channel...
    for (flu_idx in seq_len(length(flu_channels))) {

      # original: 'as long as the current index is within the few that you've decided to calibrate....'
      # removing this bc here, all fluorescence channels ought to be calibratable!
      # if(length(flu_gains) >= flu_idx){

      out_data <- fpcountr::calibrate_flu(pr_data = out_data,
                                          flu_instr = instr,
                                          flu_channel = flu_channels[flu_idx],
                                          flu_gain = flu_gains[flu_idx],
                                          flu_slug = flu_slugs[flu_idx],
                                          flu_label = flu_labels[flu_idx],
                                          do_quench_correction = do_quench_correction,
                                          conversion_factors_csv = fluor_coeffs_csv)

      # plot calibrated fluorescence data
      if(isFALSE(timecourse)){

        # heatmap - calibrated OD
        max_value <- max(out_data[[paste0("calibrated_", flu_labels[flu_idx])]], na.rm = TRUE)
        plt_flu_calib <- ggplot2::ggplot(data = out_data,
                                         ggplot2::aes(x = .data$column, y = row, fill = .data[[paste("calibrated_", flu_labels[flu_idx], sep = "")]])) +

          ggplot2::geom_tile() +
          ggplot2::scale_x_discrete("", position = "top", limits = factor(unique(out_data$column))) +
          ggplot2::scale_y_discrete("", limits = rev(unique(out_data$row))) +
          viridis::scale_fill_viridis(paste0(flu_labels[flu_idx], " (molecules)"),
                                      discrete = FALSE, limits = c(0, max_value),
                                      alpha = 0.4,
                                      na.value = "white") +

          ggplot2::geom_text(ggplot2::aes(label = scales::scientific(.data[[paste("calibrated_", flu_labels[flu_idx], sep = "")]])), na.rm = TRUE, size = 2.5) +

          ggplot2::theme_bw() + # base_size = 8
          ggplot2::theme(
            aspect.ratio = 8/12,
            panel.grid = ggplot2::element_blank()
          )
        plt_flu_calib

      } else if(isTRUE(timecourse)){

        plt_flu_calib <- ggplot2::ggplot(out_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                          y = .data[[paste("calibrated_", flu_labels[flu_idx], sep = "")]],
                                          colour = "calibrated"),
                             linewidth = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], " (molecules)")) + # FP name
          ggplot2::labs(caption = "") +
          ggplot2::scale_colour_discrete("") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8) +
          ggplot2::theme(
            aspect.ratio = 1,
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
          )
        plt_flu_calib

      }
      # plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "calib.pdf", sep = ""), basename(data_csv))
      plotname <- gsub(".csv", paste("_", flu_labels[flu_idx], "_3_calibrated.pdf", sep = ""), basename(data_csv)) # plot numbering: both plot6
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu_calib,
                      height = 16, width = 24, units = "cm")

      # original: 'as long as the current index is within the few that you've decided to calibrate....'
      # removing this bc here, all fluorescence channels ought to be calibratable!
      #   } # if(length(flu_gains) >= flu_idx){
      #   else {break} # acts to break if loop once number exceeds number of flu_channels
      # } # if else

    } # for each fluorescence channel

  } # if do calibrate

  # Save ------------------------------------------------------------------------------------------------

  filename <- gsub(".csv", "_processed.csv", basename(data_csv))
  utils::write.csv(x = out_data, file = file.path(outfolder, filename), row.names = FALSE)

  return(out_data)
}
