#' Calculate FP/cell values
#'
#' Takes as input timecourse plate reader data processed with `process_plate`
#' and uses normalised/calibrated values to calculate per-cell values. Adds
#' column(s) for FP/cell as either: (a) normalisedFP/normalisedOD (rfu/od), (b)
#' calibratedFP/calibratedOD (molecules/cell). Plots results and returns df.
#' Note that technically, units of molecules are 'molecules of equivalent FP'
#' and cells are 'particles of equivalent microspheres'.
#'
#' @param data_csv path to a csv file containing processed plate reader data
#' @param flu_channels the column names for the NORMALISED fluorescence data
#' @param flu_labels the column names for the CALIBRATED fluorescence data
#' @param remove_wells list of coordinates of wells to be removed from analysis
#'   (eg. empty wells)
#' @param get_rfu_od logical. if TRUE, uses normalised_FP and normalised_OD to
#'   calculate FP per cell as FP/OD in relative fluorescence units/relative OD
#'   (rfu/od).
#' @param get_mol_cell logical. if TRUE, uses calibrated_FP and calibrated_OD to
#'   calculate FP per cell as FP/OD in relative fluorescence units/relative OD
#'   (molecules/cell).
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @return a data.frame with columns for each FP/cell calculation
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#'
#' @examples
#' pc_data_mCherry <- calc_fppercell(data_csv = "mcherry_parsed_processed.csv")), flu_channels = c("red1red1"), flu_labels = c("mCherry"), remove_wells = c("A11"), get_rfu_od = TRUE, get_mol_cell = TRUE, outfolder = file.path("plots"))
calc_fppercell <- function(data_csv,
                           flu_channels, # colnames of normalised values (rfu), eg c("red1red1", "red1red2")
                           flu_labels, # colnames of calibrated values (molecules), eg c("mCherry", "mScarlet")
                           remove_wells,
                           get_rfu_od = TRUE, # get_rfu_pems = FALSE, # pointless
                           get_mol_cell = FALSE, # used to be get_mefl_pems
                           outfolder = "."){

  # Get parsed data -------------------------------------------------

  pr_data <- utils::read.csv(data_csv, check.names = FALSE)
  pr_data

  percell_data <- pr_data # copy over data
  percell_data

  # Remove wells eg blank wells -------------------------------------------------

  # Blank wells tend to distort scales as FP/OD leads to division by close to zero values = huge values
  if(!is.null(remove_wells)){

    ## to mutate blank wells to NA would be ideal, but complex

    ## to remove blank wells, do this
    percell_data <- percell_data %>%
      dplyr::filter(!.data$well %in% remove_wells)
    percell_data

  }

  # Location for saved plots -------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # RFU/OD calculation -------------------------------------------------

  if (get_rfu_od){
    # If get_rfu_od = TRUE, calculate FP(relative fluor units)/OD...

    # normalised_OD required
    # normalised_FPs required
    # This will be normalised_FP/normalised_OD (rfu/od).

    for (flu_idx in seq_len(length(flu_channels))) {
      # For each FP...

      # Calculate FP/OD -----------------------------

      # cols to use
      cellmeasure <- "normalised_OD"
      if(!any(names(percell_data) %in% cellmeasure)){
        # if normalised_OD column doesn't exist
        message("Normalised_OD needed to calculate FP/OD in rfu/od.
                Measure missing, skipping all rfu/od calculations.")
        break
      }
      flumeasure <- paste0("normalised_", flu_channels[flu_idx])
      if(!any(names(percell_data) %in% flumeasure)){
        # if normalised_FP column doesn't exist
        message(flumeasure, " needed to calculate FP/OD in rfu/od.
                Measure missing, skipping FP.")
        next
      }

      # Calc
      percell_data <- percell_data %>%
        dplyr::group_by(.data$time) %>%
        dplyr::mutate(v1 = (.data[[flumeasure]]) / (.data[[cellmeasure]]) )
      # for each timepoint, divides norm GFP by norm OD

      # Rename last column:
      names(percell_data)[ncol(percell_data)] <- paste0("normalised", flu_channels[flu_idx], "_perOD")
      names(percell_data)

      # Plot FP - raw, normalised -----------------------------------------

      if(any(names(percell_data) %in% flu_channels[flu_idx]) &
         any(names(percell_data) %in% paste0("normalised_", flu_channels[flu_idx]))
      ){
        # doesn't run if FP/normalised_FP column do not exist

        # plot with legend
        # plt_flu_calib <- ggplot2::ggplot(percell_data) +
        #   ggplot2::geom_line(ggplot2::aes(x = .data$time,
        #                                   y = .data[[flu_channels[flu_idx]]],
        #                                   colour = "raw"), size = 0.5) +
        #   ggplot2::geom_line(ggplot2::aes(x = .data$time,
        #                                   y = .data[[paste0("normalised_", flu_channels[flu_idx])]],
        #                                   colour = "normalised"), size = 0.5) +
        #   ggplot2::scale_x_continuous("time") +
        #   ggplot2::scale_y_continuous(name = paste0(flu_channels[flu_idx], " (rfu)"),
        #                               labels = scales::label_scientific()) +
        #   ggplot2::labs(caption = "") +
        #   ggplot2::scale_colour_discrete("") +
        #   ggplot2::facet_grid(row~column) +
        #   ggplot2::theme_bw(base_size = 8) +
        #   ggplot2::theme(
        #     aspect.ratio = 1,
        #     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
        #     panel.grid.minor = ggplot2::element_blank()
        #   )
        # plt_flu_calib

        # plot with caption
        plt_flu_calib <- ggplot2::ggplot(percell_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[flu_channels[flu_idx]]]),
                             colour = "black", size = 0.5) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data[[paste0("normalised_", flu_channels[flu_idx])]]),
                             colour = "red", size = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0(flu_channels[flu_idx], " (rfu)"),
                                      labels = scales::label_scientific(digits = 2)) +
          ggplot2::labs(caption = "black: raw, red: normalised") +
          ggplot2::facet_grid(row~column) +
          ggplot2::theme_bw(base_size = 8) +
          ggplot2::theme(
            aspect.ratio = 1,
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
            panel.grid.minor = ggplot2::element_blank()
          )
        plt_flu_calib
        plotname <- paste0("normalised", flu_channels[flu_idx], "_total.pdf")
        ggplot2::ggsave(file.path(outfolder, plotname),
                        plot = plt_flu_calib,
                        height = 16, width = 24, units = "cm")
      }

      # Plot FP/OD ---------------------------------------------------------

      # plot
      plt_flu <- ggplot2::ggplot(percell_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                        y = .data[[paste0("normalised", flu_channels[flu_idx], "_perOD")]]),
                           size = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::scale_y_continuous(name = paste0(flu_channels[flu_idx], "/OD (rfu/od)"),
                                    labels = scales::label_scientific(digits = 2)) +
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
      plt_flu
      plotname <- paste0("normalised", flu_channels[flu_idx], "_perOD.pdf")
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")

    } # then repeat for each FP

  } # get_rfu_od

  # FP/Cell in RFU/CELL --------------------------------------------------

  # pointless

  # if(get_rfu_pems){
  #   # If get_rfu_pems = TRUE, calculate FP(rfu)/Cell(pems)...
  #
  #   # calibrated_OD required
  #   # normalised_FPs required
  #   # This will be normalised_FP/calibrated_OD (rfu/pems).
  #
  #   for (flu_idx in seq_len(length(flu_channels))) {
  #     # For each FP...
  #
  #     # Calculate FP/Cell -----------------------------
  #
  #     # cols to use
  #     cellmeasure <- "calibrated_OD"
  #     if(!any(names(percell_data) %in% cellmeasure)){
  #       # if calibrated_OD column doesn't exist
  #       message("Calibrated_OD needed to calculate FP/cell in rfu/pems.
  #               Measure missing, skipping all rfu/pems calculations.")
  #       break
  #     }
  #     flumeasure <- paste0("normalised_", flu_channels[flu_idx])
  #     if(!any(names(percell_data) %in% flumeasure)){
  #       # if normalised_FP column doesn't exist
  #       message(flumeasure, " needed to calculate FP/cell in rfu/pems.
  #               Measure missing, skipping FP.")
  #       next
  #     }
  #
  #     percell_data <- percell_data %>%
  #       dplyr::group_by(.data$time) %>%
  #       dplyr::mutate(v1 = (.data[[flumeasure]]) / (.data[[cellmeasure]]) )
  #     # for each timepoint, divides norm GFP by calib OD
  #
  #     # Rename last column:
  #     names(percell_data)[ncol(percell_data)] <- paste0("normalised", flu_channels[flu_idx], "_perCell")
  #     names(percell_data)
  #
  #     # Plot FP/Cell -----------------------------
  #
  #     # plot
  #     plt_flu <- ggplot2::ggplot(percell_data) +
  #       ggplot2::geom_line(ggplot2::aes(x = .data$time,
  #                                       y = .data[[paste0("normalised", flu_channels[flu_idx], "_perCell")]]),
  #                          size = 0.5) +
  #       ggplot2::scale_x_continuous("time") +
  #       ggplot2::scale_y_continuous(name = paste0(flu_channels[flu_idx], "/cell (rfu/pems)"),
  #                                   labels = scales::label_scientific()) +
  #       # pems = particles of equiv microspheres
  #       ggplot2::labs(caption = "") +
  #       ggplot2::scale_colour_discrete("") +
  #       ggplot2::facet_grid(row~column) +
  #       ggplot2::theme_bw(base_size = 8) +
  #       ggplot2::theme(
  #         aspect.ratio = 1,
  #         legend.position = "none",
  #         axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
  #         panel.grid.minor = ggplot2::element_blank()
  #       )
  #     plt_flu
  #     plotname <- paste0("normalised", flu_channels[flu_idx], "_perCell.pdf")
  #     ggplot2::ggsave(file.path(outfolder, plotname),
  #                     plot = plt_flu,
  #                     height = 16, width = 24, units = "cm")
  #   } # then repeat for each FP
  #
  # } # get_rfu_pems

  # FP/Cell in molecules/cell --------------------------------------------------

  if(get_mol_cell){
    # If get_mol_cell = TRUE, calculate FP(molecules)/Cell...

    # calibrated_OD required
    # calibrated_FPs required
    # This will be calibrated_GFP/calibrated_OD (molecules/cell).

    for (flu_idx in seq_len(length(flu_labels))) {
      # For each FP...

      # Calculate FP/Cell -----------------------------

      # cols to use
      cellmeasure <- "calibrated_OD"
      if(!any(names(percell_data) %in% cellmeasure)){
        # if calibrated_OD column doesn't exist
        message("Calibrated_OD needed to calculate FP/cell in molecules/cell.
                Measure missing, skipping all molecules/cell calculations.")
        break
      }
      flumeasure <- paste0("calibrated_", flu_labels[flu_idx])
      if(!any(names(percell_data) %in% flumeasure)){
        # if calibrated_FP column doesn't exist
        message(flumeasure, " needed to calculate FP/cell in molecules/cell.
                Measure missing, skipping FP.")
        next
      }

      percell_data <- percell_data %>%
        dplyr::group_by(.data$time) %>%
        dplyr::mutate(v1 = (.data[[flumeasure]]) / (.data[[cellmeasure]]) )
      # for each timepoint, divides norm/calib GFP to norm/calib OD

      # Rename last column:
      names(percell_data)[ncol(percell_data)] <- paste0("calibrated", flu_labels[flu_idx], "_perCell")
      names(percell_data)

      # Plot FP - raw, normalised, calibrated ----------------------------------

      if(any(names(percell_data) %in% paste0("calibrated_", flu_labels[flu_idx]))
      ){
        # doesn't run if calibrated_FP column does not exist

        plt_flu_calib <- ggplot2::ggplot(percell_data) +
          ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                          y = .data[[paste0("calibrated_", flu_labels[flu_idx])]],
                                          colour = "calibrated"), size = 0.5) +
          ggplot2::scale_x_continuous("time") +
          ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], " (molecules)"),
                                      labels = scales::label_scientific()) +
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
        plotname <- paste0("calibrated", flu_labels[flu_idx], "_total.pdf")
        ggplot2::ggsave(file.path(outfolder, plotname),
                        plot = plt_flu_calib,
                        height = 16, width = 24, units = "cm")
      }

      # Plot FP/Cell -----------------------------

      # plot
      plt_flu <- ggplot2::ggplot(percell_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                        y = .data[[paste0("calibrated", flu_labels[flu_idx], "_perCell")]]),
                           size = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], "/cell (molecules/cell)"),
                                    labels = scales::label_scientific()) +
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
      plt_flu
      plotname <- paste0("calibrated", flu_labels[flu_idx], "_perCell.pdf")
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu,
                      height = 16, width = 24, units = "cm")
    } # then repeat for each FP

  } # get_mol_cell

  # Save CSV --------------------------------------------------

  filename <- gsub(".csv", "_pc.csv", basename(data_csv))
  utils::write.csv(x = percell_data,
                   file = file.path(outfolder, filename),
                   row.names = FALSE)

  return(percell_data)
}
