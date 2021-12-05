#' Plot absorbance spectrum
#'
#' Processes *parsed absorbance spectrum data*, collected as a *dilution series*
#' measured with an *A200-1000 absorbance spectrum*. Corrects raw data to path
#' length of 1cm by a user-defined method, and normalises to the blanks. Plots
#' spectra and returns processed data.
#'
#' @param spectrum_csv path of a .csv file of your spectrum data
#' @param subset_rows logical. should script take a subset of the rows (or whole
#'   table)? Defaults to FALSE.
#' @param rows_to_keep character array. If `subset_rows` is TRUE, script will
#'   choose rows to keep from this list. Defaults to `c("C","D")`.
#' @param columns_to_keep numeric array. If `subset_rows` is TRUE, script will
#'   choose rows to keep from this list. Defaults to `c(1,12)`. Note that blanks
#'   are required for normalisation, so if blanks are in column12, it is
#'   necessary to include 12 in this list.
#' @param xrange numerical lower and upper bounds of wavelengths, between which
#'   the is subset for plotting. This can be useful for clear plates, which have
#'   high background <300nm, so can set the xrange as c(300,800) or similar.
#' @param pl_method string denoting which method of path length normalisation to
#'   use. Options are `calc_each`, `calc_blanks` and `volume`. All three are
#'   always calculated (and compared in one of the output plots), but only the
#'   chosen method is used. The `calc_each` method uses the path length
#'   calculated from each well separately using the A900 and A975 measurements
#'   and the k-factor of the buffer used. The `calc_blanks` method uses the path
#'   length calculated from the blanks for all wells. Finally, the `volume`
#'   method uses the path length expected from the given volume, using internal
#'   path length data measured with water via the internal
#'   `fpcountR::get_pathlength()` function. Defaults to `calc_blanks`.
#' @param buffer_used string corresponding to buffer used for assay, for use in
#'   path length determination. Must be in table used by
#'   `fpcountR::get_kfactor()`. To look at table, use
#'   `fpcountR::view_kfactors()`. Defaults to "water".
#' @param concentration_used numeric value corresponding to concentration of
#'   buffer, for use in path length determination. Must be in the same units as
#'   used by table in `fpcountR::get_kfactor()`. To look at table, use
#'   `fpcountR::view_kfactors()`. Defaults to 0, which effectively means
#'   water. Needs changing if using different buffer!
#' @param temperature_used numeric value corresponding to temperature of assay,
#'   for use in path length determination. Defaults to 25.
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#'
#' @examples
#' spectrum <- plot_absorbance_spectrum(spectrum_csv = "spectrum_parsed.csv", subset_rows = TRUE, rows_to_keep = c("C","D"), columns_to_keep = c(1:12), pl_method = "calc_blanks", buffer_used = "TBS", concentration_used = 0.005, temperature_used = 30, xrange = c(250,1000), outfolder = "spectrum")

plot_absorbance_spectrum <- function(spectrum_csv,

                                     subset_rows = FALSE, rows_to_keep = c("C","D"), columns_to_keep = c(1,12),
                                     xrange = c(200,1000),

                                     # for path length calcs:
                                     pl_method = "calc_blanks", # options: "calc_each", "calc_blanks", "volume"
                                     buffer_used = "water", concentration_used = 0, temperature_used = 25,

                                     outfolder = "."){

  # Get data --------------------------------------------------------------------------------------------------

  spectrum_data <- utils::read.csv(spectrum_csv, header = TRUE, check.names = FALSE)

  # Location for saved outputs --------------------------------------------------------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # Get types of measure -----------------------------------------------------------------------------------------------------

  well_idx <- which(names(spectrum_data) == "well")
  row_idx <- which(names(spectrum_data) == "row")
  measures <- names(spectrum_data)[(well_idx+1):(row_idx-1)]
  ### columns of measurements are looked for BETWEEN "well" column and "row" column

  # Remove rows if requested --------------------------------------------------------------------------------------------------

  # always remove rows without protein annotated (empty wells)
  spectrum_data <- spectrum_data %>%
    dplyr::filter(.data$protein != "")

  # remove rows if requested
  if(subset_rows == TRUE){
    spectrum_data <- spectrum_data %>%
      dplyr::filter(.data$row %in% rows_to_keep)

    spectrum_data <- spectrum_data %>%
      dplyr::filter(.data$column %in% columns_to_keep)
  }

  # error checking - one calibrant only
  if(length(unique(spectrum_data$calibrant)) > 1){
    message("Error: data contains multiple calibrants.")
    return()
  }

  # Data transformation steps -------------------------------------------------------

  # 1. Transform data into long format
  # 2. Adjust path length to 1cm
  # 3. Normalise to buffer wells
  # 4. Summarise (take means)

  # 1. Transform data into long format  -------------------------------------------------------

  raw_values <- spectrum_data %>%
    tidyr::pivot_longer(tidyselect::all_of(measures),
                        names_to = "measure",
                        values_to = "raw_value")
  raw_values
  raw_values$measure <- as.numeric(raw_values$measure) # make measure column (of entries like "200", "201", .. "800"), numeric

  # Plots
  # Raw data full spectrum
  data.to.plot <- raw_values
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  ymax <- max(data.to.plot$raw_value)*1.1 # set y axis max
  plot1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data.to.plot, ggplot2::aes(x = measure, y = raw_value, colour = as.factor(replicate)), size = 0.5) +
    ggplot2::scale_x_continuous("wavelength (nm)") +
    ggplot2::scale_y_continuous("absorbance", limits = c(0, ymax)) +
    ggplot2::scale_colour_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ .) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot1a_raw.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  # Raw data full spectrum, blanks only
  data.to.plot <- data.to.plot %>%
    dplyr::filter(is.na(dilution))
  ymax <- max(data.to.plot$raw_value)*1.1 # set y axis max
  plot1 <- ggplot2::ggplot(data = data.to.plot) +
    ggplot2::geom_point(ggplot2::aes(x = measure, y = raw_value, colour = as.factor(replicate)), size = 1) +
    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(200,1000)) +
    ggplot2::scale_y_continuous("absorbance", limits = c(0, ymax)) +
    ggplot2::labs(title = "raw absorbance of blanks") +
    ggplot2::scale_colour_discrete("replicate") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio = 1,
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
                   strip.background = ggplot2::element_blank(),
                   # strip.text.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   panel.grid.minor = ggplot2::element_blank())
  plot1
  plotname <- "plot1b_raw_blanks.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  # Raw data 250nm+
  data.to.plot <- raw_values %>%
    dplyr::filter(measure > 250) %>%
    dplyr::filter(measure >= xrange[1] & measure <= xrange[2])
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  ymax <- max(data.to.plot$raw_value)*1.1 # set y axis max
  ggplot2::ggplot() +
    ggplot2::geom_point(data = data.to.plot, ggplot2::aes(x = measure, y = raw_value, colour = as.factor(replicate)), size = 1) +
    ggplot2::scale_x_continuous("wavelength (nm)") +
    ggplot2::scale_y_continuous("absorbance", limits = c(0, ymax)) +
    ggplot2::scale_colour_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ .) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )


  # 2. Adjust path length to 1cm ----------------

  # to convert to pathlength = 1cm, need to know pathlength
  # if data contains a900-1000, and specifically measure == 900 and measure == 975 exist, can calculate this for each well
  # this may differ for each well so calculation needs doing before wells are compared (eg. during buffer normalisation)

  # Find kfactor for buffer
  raw_values <- raw_values %>%
    dplyr::mutate(kfactor_1cm = fpcountr::get_kfactor(buffer_used = buffer_used, concentration_used = concentration_used,
                                                              temperature_used = temperature_used))
  raw_values

  # Create new df
  pl_values <- c()
  skip_calc <- FALSE

  if( (!900 %in% raw_values$measure) | (!975 %in% raw_values$measure) ){
    message("\nWarning: one or both wavelengths required for path length calculation is missing.")
    message("Skipping calculation and setting pathlength to 'volume'.")
    pl_method <- "volume"
    skip_calc <- TRUE
  }

  ## 2A. Path length calcs for each well
  for(well_i in unique(raw_values$well)){

    # create temp df. required so we can add in the extra columns
    temp_well_values <- raw_values %>%
      dplyr::filter(.data$well == well_i)
    temp_well_values

    if(skip_calc){

      # where a900 or a975 don't exist in data
      # add empty columns
      temp_well_values <- temp_well_values %>%
        dplyr::mutate(kfactor_well = NA) %>%
        dplyr::mutate(pathlength_each = NA)
      temp_well_values

    } else {

      # find !RAW! value of 900 and 975
      a900 <- temp_well_values %>%
        dplyr::filter(.data$measure == 900) %>%
        dplyr::select(.data$raw_value) %>%
        as.numeric()
      a900
      a975 <- temp_well_values %>%
        dplyr::filter(.data$measure == 975) %>%
        dplyr::select(.data$raw_value) %>%
        as.numeric()
      a975

      # find kfactor_well (a975-a900) and pathlength
      temp_well_values <- temp_well_values %>%
        dplyr::mutate(kfactor_well = a975-a900) %>%
        dplyr::mutate(pathlength_each = .data$kfactor_well/.data$kfactor_1cm)
      temp_well_values

    }

    # rbind pathlength-containing values
    pl_values <- rbind(pl_values, temp_well_values)
    pl_values

  } # for each well

  pl_values

  ## 2B. Path length calcs from blanks only
  pl_by_blanks <- pl_values %>%
    dplyr::mutate(dilution = ifelse(is.na(dilution), 0, dilution)) %>% # make blanks 'dilution = 0'
    dplyr::filter(dilution == 0) %>%
    dplyr::select(pathlength_each) %>%
    dplyr::distinct() %>% # collapse identical rows due to all the wavelengths = 2 rows
    dplyr::summarise(pathlength_each = mean(pathlength_each)) %>% # take mean of duplicate blanks
    as.numeric()
  pl_by_blanks
  pl_values <- pl_values %>%
    dplyr::mutate(pathlength_blanks = pl_by_blanks)

  ## 2C. Path length calcs from volume only
  pl_values <- pl_values %>%
    dplyr::mutate(pathlength_volume = fpcountr::get_pathlength(test_volume = .data$volume))

  # Plots
  if(!skip_calc){

    # Absorbance at 900-1000nm
    data.to.plot <- pl_values %>%
      dplyr::filter(measure > 875) %>%
      dplyr::filter(measure >= xrange[1] & measure <= xrange[2])
    data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
    newlist <- levels(data.to.plot$dilution)
    data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
    ymax <- max(data.to.plot$raw_value)*1.1 # set y axis max
    plot1 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = data.to.plot, ggplot2::aes(x = measure, y = raw_value, colour = as.factor(replicate)), size = 1) +
      ggplot2::geom_vline(xintercept = 900, colour = "black") +
      ggplot2::geom_vline(xintercept = 975, colour = "black") +
      ggplot2::scale_x_continuous("wavelength (nm)") +
      ggplot2::scale_y_continuous("absorbance", limits = c(0, ymax)) +
      ggplot2::scale_colour_discrete("replicate") +
      ggplot2::facet_wrap(dilution ~ .) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", hjust = 0),
        panel.grid.minor = ggplot2::element_blank()
      )
    plot1
    plotname <- "plot2a_a9001000.pdf"
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = plot1,
                    width = 18, height = 18, units = "cm")

    # path lengths cm
    unique(pl_values$dilution)
    data.to.plot <- pl_values %>%
      dplyr::mutate(dilution = ifelse(is.na(dilution), 0, dilution)) # include dilution = 0 on plot
    unique(data.to.plot$dilution)
    plot1 <- ggplot2::ggplot() +

      # path length using data
      ggplot2::geom_point(data = data.to.plot %>%
                            dplyr::filter(dilution > 0),
                          ggplot2::aes(x = dilution, y = pathlength_each, colour = as.factor(replicate)),
                          size = 1) +
      ggplot2::geom_point(data = data.to.plot %>%
                            dplyr::filter(dilution == 0),
                          ggplot2::aes(x = dilution, y = pathlength_each, colour = as.factor(replicate)),
                          size = 1, shape = 1) +

      # path length using mean of blanks
      ggplot2::geom_hline(yintercept = pl_by_blanks, linetype = "dashed") +

      # path length according to volume
      ggplot2::geom_hline(yintercept = fpcountr::get_pathlength(test_volume = data.to.plot$volume[1]), colour = "grey", linetype = "dashed") +

      ggplot2::scale_x_continuous("dilution") +
      ggplot2::scale_y_continuous("path length (cm)") +
      ggplot2::scale_colour_discrete("replicate") +
      ggplot2::labs(title = "path length calculations",
                    caption = "black line: pl acc to blanks\ngrey line: pl acc to volume") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", hjust = 0),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
    plot1
    plotname <- "plot2b_pathlengths.pdf"
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = plot1,
                    width = 12, height = 12, units = "cm")

  } else {

    # If a900-1000 not available, don't plot anything, just print the path length acc to the volume.
    message(paste0("\nPath length for ", pl_values$volume[1], "ul is ", round(fpcountr::get_pathlength(test_volume = pl_values$volume[1]), 3), "cm."))

  }

  ##

  ## 2D. Path length normalisation decision & calculation

  message(paste0("\nCalculating path lengths using chosen method: ", pl_method, "."))
  if(pl_method == "calc_each"){
    message("Path length will calculated per well from the data.")
    pl_values <- pl_values %>%
      dplyr::mutate(pathlength_method = pl_method) %>% # record method
      dplyr::mutate(pathlength = pathlength_each) # assign chosen pathlength calc as "pathlength" column
  }
  if(pl_method == "calc_blanks"){
    message("Path length will calculated from the blanks data.")
    pl_values <- pl_values %>%
      dplyr::mutate(pathlength_method = pl_method) %>% # record method
      dplyr::mutate(pathlength = pathlength_blanks) # assign chosen pathlength calc as "pathlength" column
  }
  if(pl_method == "volume"){
    message("Path length will be based on volume.")
    message(paste0("Volume used = ", pl_values$volume[1], " ul."))
    pl_values <- pl_values %>%
      dplyr::mutate(pathlength_method = pl_method) %>% # record method
      dplyr::mutate(pathlength = pathlength_volume) # assign chosen pathlength calc as "pathlength" column
  }
  pl_values <- pl_values %>%
    dplyr::mutate(raw_cm1_value = raw_value/pathlength) # calculate cm-1 values

  # Plots
  # Raw data cm_1 250nm+
  data.to.plot <- pl_values  %>%
    dplyr::filter(measure > 250) %>%
    dplyr::filter(measure >= xrange[1] & measure <= xrange[2])
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  ymax <- max(data.to.plot$raw_cm1_value)*1.1 # set y axis max
  plot1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data.to.plot,
                        ggplot2::aes(x = measure, y = raw_cm1_value, colour = as.factor(replicate)),
                        size = 1) +
    ggplot2::scale_x_continuous("wavelength (nm)") +
    ggplot2::scale_y_continuous("absorbance (cm-1)", limits = c(0, ymax)) +
    ggplot2::scale_colour_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ .) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot2c_rawcm1.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  # 3. Normalise to buffer wells -------------------------------------------------------

  norm_values <- c() # to be populated by rbind of gain-by-gain filters of prev df

  for(meas in measures){

    # create temp df. required so we can add in the extra columns
    temp_meas_values <- pl_values %>%
      dplyr::filter(.data$measure == meas) # filter rows by measure/gain
    temp_meas_values

    # find mean of blanks
    blanks <- temp_meas_values %>%
      dplyr::filter(.data$protein == "none") # find blanks
    blanks
    mean_raw_cm1_blanks <- mean(blanks$raw_cm1_value, na.rm = TRUE) # calc mean
    mean_raw_cm1_blanks

    # add new columns
    temp_meas_values <- temp_meas_values %>%
      dplyr::mutate(raw_cm1_blanks = mean_raw_cm1_blanks) %>%
      dplyr::mutate(normalised_cm1_value = .data$raw_cm1_value - raw_cm1_blanks)
    temp_meas_values

    # rbind normalised values
    norm_values <- rbind(norm_values, temp_meas_values)
    norm_values

  } # measures

  # Plots
  # Normalised cm_1 250nm+
  data.to.plot <- norm_values  %>%
    dplyr::filter(measure > 250 & measure < 800) %>%
    dplyr::filter(measure >= xrange[1] & measure <= xrange[2])
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  ymax <- max(data.to.plot$normalised_cm1_value)*1.1 # set y axis max

  plot1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data.to.plot,
                        ggplot2::aes(x = measure, y = normalised_cm1_value, colour = as.factor(replicate)),
                        size = 1) +
    ggplot2::scale_x_continuous("wavelength (nm)") +
    ggplot2::scale_y_continuous("absorbance (norm, cm-1)") +
    ggplot2::scale_colour_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ .) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot3_normcm1.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  # 4. Summarise (take means) -------------------------------------------------------

  # Take mean of duplicate readings for: raw_value, normalised_value, normalised_cm1_value
  names(norm_values)
  summ_values <- norm_values %>%
    dplyr::group_by(.data$instrument, .data$plate, .data$seal,
                    # .data$channel_name, .data$channel_ex, .data$channel_em, # removed
                    .data$media, .data$calibrant, .data$protein,
                    # replicate will vary
                    .data$mw_gmol1, .data$concentration_ngul, .data$dilution, .data$rev_dilution, .data$volume,
                    # well, row, column will vary
                    .data$measure,
                    ## raw_value # TAKING MEAN
                    .data$kfactor_1cm,
                    # .data$kfactor_well, .data$pathlength_each, # will vary
                    # .data$pathlength_blanks, # remove as removing pathlength_each
                    # .data$pathlength_volume, # remove as removing pathlength_each
                    .data$pathlength_method,
                    # .data$pathlength, # might vary
                    ## raw_cm1_value, # TAKING MEAN
                    .data$raw_cm1_blanks,
                    ## normalised_cm1_value # TAKING MEAN
                    .drop = FALSE) %>%
    dplyr::summarise(dplyr::across(dplyr::ends_with("_value"), ~mean(.x, na.rm = TRUE)))
  summ_values

  # Normalised cm_1 250nm+
  data.to.plot <- summ_values %>%
    dplyr::filter(measure > 250 & measure < 800) %>%
    dplyr::filter(measure >= xrange[1] & measure <= xrange[2])
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  ymax <- max(data.to.plot$normalised_cm1_value)*1.1 # set y axis max
  plot1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data.to.plot,
                        ggplot2::aes(x = measure, y = normalised_cm1_value),
                        size = 1) +
    ggplot2::scale_x_continuous("wavelength (nm)") +
    ggplot2::scale_y_continuous("absorbance (norm, cm-1)", limits = c(0, ymax)) +
    # ggplot2::scale_colour_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ .) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot4_mean_normcm1.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  # 5. Save CSVs and Return -------------------------------------------------------

  # Save complete df norm_values as "_parsed_processed.csv" for use in get_conc_ functions
  csvname <- gsub(".csv", paste0("_processed.csv"), basename(spectrum_csv))
  utils::write.csv(x = norm_values, file = file.path(outfolder, csvname), row.names = FALSE)

  # Save summarised df summ_values as "_parsed_processed_summ.csv"
  csvname <- gsub(".csv", paste0("_processed_summ.csv"), basename(spectrum_csv))
  utils::write.csv(x = summ_values, file = file.path(outfolder, csvname), row.names = FALSE)

  return(summ_values)

}
