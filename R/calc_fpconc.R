#' Calculate FP concentration in molar units
#'
#' Takes as input timecourse plate reader data processed with `process_plate`
#' and uses normalised/calibrated values to calculate FP concentration. Adds
#' column(s) for FP concentration as either: (a) normalisedFP/cellvolume
#' (rfu/L), (b) calibratedFP/cellvolume (moles/L, or M). Plots results and
#' returns df.
#'
#' @param data_csv path to a csv file containing processed plate reader data
#' @param flu_channels the column names for the NORMALISED fluorescence data
#' @param flu_labels the column names for the CALIBRATED fluorescence data
#' @param remove_wells list of coordinates of wells to be removed from analysis
#'   (eg. empty wells)
#' @param get_rfu_vol logical. if TRUE, uses normalised_FP and OD-specific cell
#'   volume (`od_specific_total_volume`, specified in ul) to calculate FP
#'   concentration as FP/cellvolume in relative fluorescence units/litre
#'   (rfu/L).
#' @param get_mol_vol logical. if TRUE, uses calibrated_FP and OD-specific cell
#'   volume (`od_specific_total_volume`, specified in ul) to calculate FP
#'   concentration as moles/L in Molar (M) units.
#' @param od_specific_total_volume numeric. OD600-specific total cellular volume
#'   in ul*OD-1*cm, ie. the total cellular volume represented by 1 OD600 unit
#'   (in 1 cm path length). Recommended value is 3.6 from Volkmer et al., 2011.
#' @param odmeasure character. Which OD measurement is being used in the data?
#'   Specifically, which measurement is represented by the 'normalised_OD_cm1'
#'   column? eg. "OD600" or "OD700". Purpose is to record this in the table.
#' @param odmeasure_conversion numeric. How to convert the measurement specified
#'   by odmeasure to OD600? ie. OD600 = ODused / x. Use '1' for OD600 (no
#'   conversion) and 0.79 for OD700.
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @return a data.frame with columns for each FP/cell calculation
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#'
#' @examples
#' conc_data_mCherry <- calc_fpconc(data_csv = "mcherry_parsed_processed.csv", flu_channels = c("red1red1"), flu_labels = c("mCherry"), remove_wells = c("A11"), get_rfu_vol = TRUE, get_mol_vol = TRUE, od_specific_total_volume = 3.6, odmeasure = "OD700", odmeasure_conversion = 0.79, outfolder = file.path("plots"))
calc_fpconc <- function(data_csv,
                        flu_channels, # colnames of normalised values (rfu), eg c("red1red1", "red1red2")
                        flu_labels, # colnames of calibrated values (molecules), eg c("mCherry", "mScarlet")
                        remove_wells,
                        get_rfu_vol = TRUE, # rfu per volume (rfu/L)
                        get_mol_vol = FALSE, # molecules per volume (Molar)
                        od_specific_total_volume = NULL, # 3.6 ul od-1cm-1 from 2011 Volkmer.
                        odmeasure = NULL, # "OD600" or "OD700" # which OD measurement is being used in the data?
                        odmeasure_conversion = NULL, # how to convert odmeasure to OD600 # typically 1 for OD600, 0.79 for OD700
                        outfolder = "."){

  # Messages -------------------------------------------------

  # error checking
  if(is.null(od_specific_total_volume) | is.null(odmeasure) | is.null(odmeasure_conversion)){
    message("Error: specify arguments `od_specific_total_volume`, `odmeasure` and `odmeasure_conversion`")
    return()
  }

  message("Note that the default 'OD-specific total volume' is 3.6 ul per (OD600 cm-1) - and requires measurement in OD600 or a conversion to estimate OD600 from the OD used.")
  message(paste0("Using OD-specific total cell volume: ", od_specific_total_volume, "ul per (OD600 cm-1)."))
  message("The empirical ratio between E. coli absorbance at OD700/OD600 is typically 0.79.")
  message(paste0("Using OD: ", odmeasure, "."))
  message(paste0("Using conversion: ", odmeasure_conversion, "."))

  # Get parsed data -------------------------------------------------

  pr_data <- utils::read.csv(data_csv, check.names = FALSE)
  pr_data

  percell_data <- pr_data # copy over data
  percell_data

  # Remove wells eg blank wells -------------------------------------------------

  # Blank wells tend to distort scales as FP/cell volume leads to division by close to zero values = huge values
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

  # Correction of OD value to use if not using OD600 -------------------------------------------------

  ## record odmeasure used and correct it if necessary
  percell_data <- percell_data %>%
    dplyr::mutate(odmeasure_used = odmeasure) %>% # record od measure used
    dplyr::mutate(odmeasure_conversion_to_od600 = odmeasure_conversion) %>% # and whether to correct it to od600

    # OD-specific cell volume is for an OD (cm-1) type value
    # this should already exist in the data as 'normalised_OD_cm1' provided that the volume was specified, as this is calculated by process_plate()
    dplyr::mutate(normalised_OD600_cm1 = .data$normalised_OD_cm1 / odmeasure_conversion) # corrected OD: eg. OD600_cm1 = OD600_cm1 / 1 .. eg. OD600_cm1 = OD700_cm1 / 0.79
  percell_data[1,]

  # Calculation of total cellular volume -------------------------------------------------

  ## od_specific_total_volume should be specified in arguments.
  # default is 3.6 ul od-1cm-1 from Volkmer et al., 2011.
  # this is basically the cellular volume 'per OD measured'.
  # od_specific_total_volume # units = ul / (od/cm) = ul/ (odcm-1) = ul*od-1*cm
  # od_specific_total_volume/1e6 # L # convert to Litres. units = L*od-1*cm

  # add to table
  percell_data <- percell_data %>%
    dplyr::mutate(od600_specific_total_volume_L_OD1_cm = od_specific_total_volume/1e6 ) # units = L*od-1*cm
  percell_data[1,]

  ## calculate total cell volume (L)
  percell_data <- percell_data %>%
    dplyr::mutate(total_cell_volume_L = .data$normalised_OD600_cm1 * .data$od600_specific_total_volume_L_OD1_cm) # cell volume = OD600*(volume per OD600)
  percell_data[1,]

  # RFU/cell volume calculation -------------------------------------------------

  if (get_rfu_vol){
    # If get_rfu_vol = TRUE, calculate FP(relative fluor units)/cellvolume...

    # normalised_FPs required
    # This will be normalised_FP/cellvolume (rfu/cellvolume).

    for (flu_idx in seq_len(length(flu_channels))) {
      # For each FP...

      # Calculate FP concentration (rfu/vol) -----------------------------

      # cols to use

      # fluorescence:
      flumeasure <- paste0("normalised_", flu_channels[flu_idx])
      if(!any(names(percell_data) %in% flumeasure)){
        # if normalised_FP column doesn't exist
        message(flumeasure, " needed to calculate FP/vol in rfu/vol.
                Measure missing, skipping FP.")
        next
      }

      # Calc
      percell_data <- percell_data %>%
        dplyr::group_by(.data$time) %>%
        dplyr::mutate(v1 = (.data[[flumeasure]]) / total_cell_volume_L )
      # for each timepoint, divides norm GFP by cell volume
      as.data.frame(percell_data)[1,]

      # Rename last column:
      names(percell_data)[ncol(percell_data)] <- paste0("normalised_", flu_labels[flu_idx], "_percellvolume")
      names(percell_data)

      # Plot FP - normalised ----------------------------------

      plt_flu_calib <- ggplot2::ggplot(percell_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                        y = .data[[paste0("normalised_", flu_labels[flu_idx], "_percellvolume")]],
                                        colour = "normalised"), size = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], " concentration (rfu/cellvolume)")) +
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
      plotname <- paste0("normalised_", flu_labels[flu_idx], "_concentration.pdf")
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu_calib,
                      height = 16, width = 24, units = "cm")

    } # then repeat for each FP

  } # get_rfu_vol

  # FP/cell volume in molar units --------------------------------------------------

  if(get_mol_vol){
    # If get_mol_vol = TRUE, calculate FP concentration (M)

    # calibrated_FPs required
    # This will be calibrated_GFP moles/L (moles/L or M).

    for (flu_idx in seq_len(length(flu_labels))) {
      # For each FP...

      # Calculate FP concentration (M) -----------------------------

      # cols to use

      # fluorescence:
      # 1 Molar = 1 mole/L.
      # #molecules we know
      # #moles = molecules / (molecules/mol) = molecules/avogradro's constant
      # #molar: 1M = 1mol/L
      flumeasure <- paste0("calibrated_", flu_labels[flu_idx])
      if(!any(names(percell_data) %in% flumeasure)){
        # if calibrated_FP column doesn't exist
        message(flumeasure, " needed to calculate FP/cell in molecules/cell.
                Measure missing, skipping FP.")
        next
      }

      # add col for moles
      # #moles = molecules / (molecules/mol) = molecules/avogradro's constant
      percell_data <- percell_data %>%
        dplyr::mutate(v1 = (.data[[flumeasure]]) / (6.022*10^23) )
      percell_data
      # Rename last column:
      names(percell_data)[ncol(percell_data)] <- paste0(flumeasure, "_moles")
      names(percell_data)
      moles_columnname <- names(percell_data)[ncol(percell_data)] # save for calc
      moles_columnname

      # Calc
      # #molar: 1M = 1mol/L
      percell_data <- percell_data %>%
        dplyr::group_by(.data$time) %>%
        dplyr::mutate(v1 = (.data[[moles_columnname]]) / total_cell_volume_L )
      # for each timepoint, divides calib GFP by cell volume
      as.data.frame(percell_data)[1,]

      # Rename last column:
      names(percell_data)[ncol(percell_data)] <- paste0("calibrated_", flu_labels[flu_idx], "_Molar")
      names(percell_data)

      # Plot FP - calibrated ----------------------------------

      plt_flu_calib <- ggplot2::ggplot(percell_data) +
        ggplot2::geom_line(ggplot2::aes(x = .data$time,
                                        y = .data[[paste0("calibrated_", flu_labels[flu_idx], "_Molar")]]*1e6,
                                        colour = "calibrated"), size = 0.5) +
        ggplot2::scale_x_continuous("time") +
        ggplot2::scale_y_continuous(name = paste0(flu_labels[flu_idx], " concentration (uM)")) +
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
      plotname <- paste0("calibrated_", flu_labels[flu_idx], "_concentration.pdf")
      ggplot2::ggsave(file.path(outfolder, plotname),
                      plot = plt_flu_calib,
                      height = 16, width = 24, units = "cm")

    } # then repeat for each FP

  } # get_mol_vol

  # Save CSV --------------------------------------------------

  filename <- gsub(".csv", "_pc.csv", basename(data_csv))
  utils::write.csv(x = percell_data,
                   file = file.path(outfolder, filename),
                   row.names = FALSE)

  return(percell_data)
}
