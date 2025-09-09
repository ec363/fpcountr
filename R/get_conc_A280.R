#' Get FP concentrations using A280 method
#'
#' Get protein's concentration from a dilution series measured with an
#' absorbance spectrum. Expects 'processed' data such as that produced by
#' `process_absorbance_spectrum()`, with a file name ending `_processed.csv`,
#' which contains values corrected for path length and normalised to blanks as a
#' column called `normalised_cm1_value`, but retains replicate data containing
#' positional (well) information required for exporting predicted concentrations
#' at the end of this function. Uses `get_extcoeff_a280()` to get EC in M-1cm-1
#' and wavelength, and converts it to an EC mass extinction coefficient in
#' `(mgml)-1cm-1` using the MW (worked out from `protein_seq` and
#' `fpcountr::get_mw`). Then the function uses the `EC_A280_mgml` to work out
#' the concentration of protein in each well, using three correction methods.
#' Instead of using the normalised data directly, the values used are based on a
#' LOESS fit through the absorption spectra to minimise fluctuations due to
#' noise. Finally, linear models are fitted to each concentration prediction
#' method, and a dataframe is built, returned and saved, containing predicted
#' concentrations according to the user's chosen correction method. Plots
#' showing each of the analytical steps are saved concurrently.
#'
#' @param protein_slug character string of protein name in 'slug' form to match
#'   slug of FPbase entry.
#' @param protein_seq character string of protein sequence using 1-letter code.
#'   Required for MW calculation.
#' @param buffer character string of buffer. Optional. Defaults to "".
#' @param processed_spectrum_csv Path to CSV file of a processed absorbance
#'   spectrum. Processing should be done with `process_absorbance_spectrum()`,
#'   which corrects for path lengths and normalises to blank wells.
#' @param wells_to_remove list of wells to remove before analysis. Defaults to
#'   NULL.
#' @param disulphides required for calculation of A280 extinction coefficient.
#'   logical. Does protein have disulphides? Defaults to FALSE.
#' @param showWarnings required for calculation of A280 extinction coefficient.
#'   logical. Should function show warnings? Defaults to TRUE.
#' @param showMessages required for calculation of A280 extinction coefficient.
#'   logical. Should function show messages? Defaults to TRUE.
#' @param corr_method string corresponding to type of correction method to use
#'   for the data to remove contribution of light scatter. Options are `none`,
#'   `baseline` and `scatter`. Method `none` applies no correction. Method
#'   `baseline` subtracts the absorbance value of the wavelength supplied in
#'   `wav_to_use1`. Method `scatter` subtracts a fraction of the absorbance
#'   value of the wavelength supplied in `wav_to_use2` according to scatter
#'   theory (details in section).
#' @param wav_to_use1 numerical value of wavelength (nm) to use for `baseline`
#'   correction. Defaults to 340nm.
#' @param wav_to_use2 numerical value of wavelength (nm) to use for `scatter`
#'   correction. Defaults to 333nm.
#' @param outfolder path to folder where output files should be saved.
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   a280_concs <- get_conc_a280(
#'     protein_slug = "mcherry", protein_seq = protein_seq, buffer = "T5N15_pi",
#'     processed_spectrum_csv = "abs_parsed_processed.csv",
#'     corr_method = "scatter", wav_to_use1 = 340, wav_to_use2 = 315,
#'     outfolder = "protquant_a280/mCherry_T5N15pi"
#'   )
#' }
get_conc_a280 <- function(protein_slug, protein_seq, buffer = "",
                          processed_spectrum_csv, wells_to_remove = NULL,
                          disulphides = FALSE, showWarnings = TRUE, showMessages = FALSE,
                          corr_method = "none", # "none", "baseline", "scatter"
                          wav_to_use1 = 340, wav_to_use2 = 333,
                          outfolder = ""
){

  # Get data -------------------------------------------------

  spectrum_data <- utils::read.csv(processed_spectrum_csv, header = TRUE, check.names = FALSE)

  # Where protein is none, dilution is usually NA, so these rows are lost in later steps.
  # To keep them, they need a value.
  spectrum_data <- spectrum_data %>%
    dplyr::mutate(dilution = ifelse(.data$protein == "none", 0, .data$dilution))
  spectrum_data
  unique(spectrum_data$dilution)

  # Remove wells if requested -------------------------------------------------

  if(!is.null(wells_to_remove)){
    spectrum_data_subset <- spectrum_data %>%
      dplyr::filter(!.data$well %in% wells_to_remove)
    # unique(spectrum_data_subset$well)
  } else {
    spectrum_data_subset <- spectrum_data
  }

  # Location for saved outputs -------------------------------------------------

  # check if parent directory exists
  parent_folder <- dirname(outfolder)
  if(!dir.exists(parent_folder)){
    message("Error: Please specify a valid path for the location 'outfolder' where files should be saved.")
    return()
  }

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # 1. Replicates data -------------------------------------------------

  # Check input data

  # Plot absorbance spectra (250-800nm) of replicates separately as sanity check for reproducibility
  data.to.plot <- spectrum_data_subset %>%
    dplyr::filter(.data$measure > 250 & .data$measure < 800)
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value,
                                     colour = as.factor(replicate)), # where you want to plot >1 replicate
                        size = 0.5) +
    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(250,800)) +
    ggplot2::scale_y_continuous("absorbance (cm-1)") +
    ggplot2::scale_color_discrete("replicate") +
    ggplot2::facet_wrap(dilution ~ ., scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
  plot1
  plotname <- "plot1_abs_spectra_replicates.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  # 2. Averaged data -------------------------------------------------

  # Take mean of duplicate readings for: raw_value, normalised_value, normalised_cm1_value
  names(spectrum_data_subset)
  summ_data <- spectrum_data_subset %>%
    dplyr::group_by(.data$instrument, .data$plate, .data$seal,
                    # .data$channel_name, .data$channel_ex, .data$channel_em, # removed
                    .data$media, .data$calibrant, .data$protein,
                    # replicate will vary

                    .data$mw_gmol1, .data$concentration_ngul, .data$dilution, .data$rev_dilution,

                    # well, row, column will vary
                    .data$measure,

                    ## raw_value # TAKING MEAN
                    .data$kfactor_1cm,
                    # .data$kfactor_well, .data$pathlength_each, # will vary
                    # .data$pathlength_blanks, # remove as removing pathlength_each
                    .data$volume,
                    # .data$pathlength_volume, # remove as removing pathlength_each
                    .data$pathlength_method,
                    # .data$pathlength, # might vary
                    ## raw_cm1_value, # TAKING MEAN
                    .data$raw_cm1_blanks,
                    ## normalised_cm1_value # TAKING MEAN

                    .drop = FALSE) %>%
    dplyr::summarise(dplyr::across(dplyr::ends_with("_value"), ~mean(.x, na.rm = TRUE)))
  names(summ_data)
  summ_data # half the rows bc going from duplicates to averaged rows

  ##

  # 3. Smooth data -------------------------------------------------

  # Instead of taking the raw absorbance values, fit a loess model through the points to make the quantifications less sensitive to noise

  # Plot absorbance spectrum loess with 'geom_smooth' (250-800nm)
  data.to.plot <- summ_data %>%
    dplyr::filter(.data$measure > 250 & .data$measure < 800)
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value),
                        colour = "lightblue") +
    ggplot2::geom_smooth(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value),
                         colour = "black",
                         span = 0.5/5.5 ### could leave as default, us 0.5, or add as argument in function # 0.5 used for 250-350
    ) +
    ggplot2::geom_vline(xintercept = 280, colour = "red") +

    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(250,800)) +
    ggplot2::scale_y_continuous("absorbance (cm-1)") +
    ggplot2::facet_wrap(dilution ~ ., scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
  plot1
  plotname <- "plot3a_abs_spectra_geomsmooth.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  # Plot absorbance spectrum loess with 'geom_smooth' (250-350nm)
  data.to.plot <- summ_data %>%
    dplyr::filter(.data$measure > 250 & .data$measure < 350)
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value),
                        colour = "lightblue") +
    ggplot2::geom_smooth(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value),
                         colour = "black",
                         span = 0.5 # 0.5 used for 250-350
    ) +
    ggplot2::geom_vline(xintercept = 280, colour = "red") +

    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(250,350)) +
    ggplot2::scale_y_continuous("absorbance (cm-1)") +
    ggplot2::facet_wrap(dilution ~ ., scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
  plot1
  plotname <- "plot3a2_uv_spectra_geomsmooth.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  ## Fit loess model for absorbance ~ wavelength (separately for each dilution)

  # method slightly complex. using both:
  # https://r4ds.had.co.nz/many-models.html
  # https://stackoverflow.com/questions/50163106/loess-regression-on-each-group-with-dplyrgroup-by

  # subset
  spectrum_data_subset <- summ_data %>%
    dplyr::filter(.data$measure > 250) %>%
    dplyr::filter(.data$measure < 800)
  # nest
  spectrum_data_nested <- spectrum_data_subset %>%
    dplyr::group_by(.data$instrument, .data$media, .data$calibrant, .data$dilution) %>%
    tidyr::nest()
  spectrum_data_nested
  spectrum_data_nested$data
  # purrr map to make model
  spectrum_data_nested <- spectrum_data_nested %>%
    dplyr::mutate(model = purrr::map(.data$data, stats::loess, formula = normalised_cm1_value ~ measure, span = 0.5/5.5)) # 0.5 for 250-350, 0.5/5.5 for 250-800
  spectrum_data_nested
  # extract fitted values!!
  spectrum_data_nested <- spectrum_data_nested %>%
    dplyr::mutate(fitted_cm1_value = purrr::map(.data$model, `[[`, "fitted"))
  spectrum_data_nested
  # remove model column (awkward and no longer necessary) and unnest
  spectrum_data_model <- spectrum_data_nested %>%
    dplyr::select(-.data$model) %>%
    tidyr::unnest(cols = c(.data$data, .data$fitted_cm1_value))
  spectrum_data_model
  spectrum_data_model$fitted_cm1_value

  # Model check, entire spectrum:
  # Plot loess with geom_line - effectively make line plot of the loess fit using fitted_values column
  # to check my loess fits are the same as what geom_smooth suggests:
  data.to.plot <- spectrum_data_model
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value), colour = "lightblue") +
    # geom_line of extracted loess fitted values
    ggplot2::geom_line(ggplot2::aes(x = .data$measure, y = .data$fitted_cm1_value)) +
    ggplot2::geom_vline(xintercept = 280, colour = "red") +
    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(250,800)) +
    ggplot2::scale_y_continuous("absorbance (cm-1)") +
    ggplot2::facet_wrap(dilution ~ ., scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
  plot1
  plotname <- "plot3b_abs_spectra_modelcheck.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  # Model check, uv:
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = .data$measure, y = .data$normalised_cm1_value), colour = "lightblue") +
    # geom_line of extracted loess fitted values
    ggplot2::geom_line(ggplot2::aes(x = .data$measure, y = .data$fitted_cm1_value)) +
    ggplot2::geom_vline(xintercept = 280, colour = "red") +
    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(250,350)) +
    ggplot2::scale_y_continuous("absorbance (cm-1)") +
    ggplot2::facet_wrap(dilution ~ ., scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
  plot1
  plotname <- "plot3b2_uv_spectra_modelcheck.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 18, height = 18, units = "cm")

  ##

  # 4. Get extinction coefficient ------------------------------------------------------------

  ## 4a. Get protein sequence
  protein_seq

  ## 4b. Get MW
  protein_mw <- fpcountr::get_mw(protein = protein_seq)
  protein_mw

  ## 4c. Get mgml extinction coefficient
  ectable <- get_extcoeff_a280(protein = protein_seq, disulphides = disulphides,
                               showWarnings = showWarnings, showMessages = showMessages,
                               protein_name = protein_slug, buffer = buffer, mol_weight = protein_mw,
                               outfolder = outfolder)
  ectable

  ## 4d. Add 'EC of 1mg/ml of protein' to table
  # EC of 0.1% solution = 0.1g/100ml = 0.001g/ml = 1mg/ml
  # EC/MW = (M-1cm-1)/(g/mol) = (mol/L)-1 * (cm)-1 / (g/mol) = (L/mol)*(1/cm)*(g/mol) = L/g * 1/cm = (g/L)-1(cm)-1 = (mg/ml)-1(cm)-1
  # = EC (mg/ml)-1(cm)-1

  df_2 <- spectrum_data_model %>%
    dplyr::mutate(EC_A280_mgml = ectable$E_0.1pc)
  utils::head(df_2)

  ##

  # 5. Calculate protein concentration ------------------------------------------------------------

  # 5a. Correction method - none ------------------------------------------------------------

  # Standard method with no correction

  df_3 <- c()

  for(diln in unique(df_2$dilution)){

    # Subset
    temp_diln_values <- df_2 %>%
      dplyr::filter(.data$dilution == diln)

    # Abs
    abs_280 <- temp_diln_values %>%
      dplyr::filter(.data$measure == 280) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_280

    # Add this Abs value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(a280_std = abs_280)
    temp_diln_values

    # Get conc from this Abs and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_a280_std = .data$a280_std/.data$EC_A280_mgml) # quicker
    temp_diln_values

    # Rbind
    df_3 <- rbind(df_3, temp_diln_values)

  }

  df_3

  # Plot processing steps
  plot1 <- ggplot2::ggplot(data = df_3 %>%
                             dplyr::filter(.data$measure == 280)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution, y = .data$raw_value), colour = "grey") +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution, y = .data$raw_cm1_value), colour = "black") +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution, y = .data$normalised_cm1_value), colour = "red") +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution, y = .data$fitted_cm1_value), colour = "blue") +
    ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
    ggplot2::scale_y_continuous("absorbance") +
    ggplot2::labs(subtitle = "raw data and processing steps",
                  caption = "raw: grey, raw_cm1: black,\nnormalised_cm1: red, fitted: blue") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot5a_a280.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  # Subset data to get rid of negatives:
  df_3_subset <- subset(df_3, df_3$conc_a280_std >= 0)
  # Fit model as (Y ~ X):
  model1 <- stats::lm(conc_a280_std ~ dilution + 0, data = df_3_subset) # force through zero to compare points to ideal
  model1
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_3,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_std*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_3_subset,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_std*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_3, # to extend line to zero
                       ggplot2::aes(x = .data$dilution,
                                    y = (stats::predict(model1, df_3))*1000),
                       colour = "red") +
    ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
    ggplot2::scale_y_continuous("concentration (ng/ul)") +
    ggplot2::labs(subtitle = "predicted conc (correction = none)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot5b_a280_stdmethod.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  ##

  # 5b. Correction method1 - baseline ------------------------------------------------------------------------

  # Correction method1 to remove scatter
  # Baseline aka Nanodrop method: normalise to A340
  # So use Abs = abs_a280 - abs_340

  # ND-1000 Spectrophotometer V3.5 User’s Manual > Section 8- Protein A280 Spectrum Normalization
  # The baseline is automatically set to the absorbance value of the sample at 340 nm,
  # which should be very nearly zero absorbance. All spectra are referenced off of this zero.

  # Edit: wavelength used here defined by wav_to_use1, with default 340nm.

  df_4 <- c()

  for(diln in unique(df_3$dilution)){

    # Subset
    temp_diln_values <- df_3 %>%
      dplyr::filter(.data$dilution == diln)

    # max
    abs_280 <- temp_diln_values %>%
      dplyr::filter(.data$measure == 280) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_280
    abs_340 <- temp_diln_values %>%
      dplyr::filter(.data$measure == wav_to_use1) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_340

    # Add corrected value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(wav_corr1 = wav_to_use1) %>% # record wav_to_use1
      dplyr::mutate(abs280_corr1 = abs_280 - abs_340) # normalise
    temp_diln_values

    # Get conc from this Abs and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_a280_corr1 = .data$abs280_corr1/.data$EC_A280_mgml)
    temp_diln_values

    # Rbind
    df_4 <- rbind(df_4, temp_diln_values)

  }

  df_4

  # Subset data to get rid of negatives:
  df_4_subset <- subset(df_4, df_4$conc_a280_corr1 >= 0)
  # Fit model as (Y ~ X):
  model2 <- stats::lm(conc_a280_corr1 ~ dilution + 0, data = df_4_subset) # force through zero to compare points to ideal
  model2
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_4,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_corr1*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_4_subset,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_corr1*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_4, # to extend line to zero
                       ggplot2::aes(x = .data$dilution,
                                    y = (stats::predict(model2, df_4))*1000),
                       colour = "red") +
    ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
    ggplot2::scale_y_continuous("concentration (ng/ul)") +
    ggplot2::labs(subtitle = "predicted conc (correction = baseline)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot5c_a280_baselinenorm.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  ##

  # 5c. Correction method2 - scatter ------------------------------------------------------------

  # Correction method2 to remove scatter
  # Scatter method:
  # For A280 measurements, A333 gives 1/2 the scatter as A280, so use Abs = abs_280 - 2*abs_333
  # For other wavelengths, calculate scatter_ratio between A333 and A280, then use Abs = abs_280 - scatter_ratio*abs_333
  # Edit: wavelength used here defined by wav_to_use2, with default 333nm.

  # ## Correcting for A280 - how and why?
  # # the optical density due to light scatter varies with 1/wavelength^4
  # # y = 1/x^4
  # eq = function(x){1/x^4}
  # ggplot2::ggplot(data = data.frame(x = c(200,800)),
  #                 ggplot2::aes(x = x)) +
  #   ggplot2::stat_function(fun = eq) +
  #   ggplot2::scale_x_continuous("wavelength (nm)") +
  #   ggplot2::scale_y_continuous("light scatter (au)") +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(
  #     aspect.ratio = 1,
  #     panel.grid.minor = ggplot2::element_blank()
  #   )
  # # light scatter at a280 is not same level as light scatter elsewhere
  # scatter_a280 <- eq(280)
  # scatter_a280 # 1.626926e-10
  # # to find light scatter at a280 when there's a peak in the way,
  # # need to find it by looking for light scatter at a different wavelength, then doing a transformation
  # # don't want to use wavelengths lower than A280 bc proteins and plastics absorb highly here
  # # using higher wavelengths with lower scatter works better
  # # eg. where is light scatter 2*less than at a280?
  # scatter_to_find <- scatter_a280/2
  # scatter_to_find # 8.134631e-11
  # # if y = 1/x^4
  # # x^4 = 1/y
  # # x = 4√(1/y) = (1/y)^(1/4)
  # wavelength_to_use <- (1/scatter_to_find)^(1/4)
  # wavelength_to_use # 332.978
  # # use A333! # this is good because it isn't so high that it would contribute to absorbance of an FP

  # Work out scatter_ratio.
  # What is ratio between scatter at A280 and scatter at the user-defined wav_to_use2 wavelength?
  eq = function(x){1/x^4}
  scatter_a333 <- eq(wav_to_use2) # 333 by default
  scatter_a333
  scatter_a280 <- eq(280)
  scatter_a280
  scatter_ratio <- scatter_a280/scatter_a333
  scatter_ratio # the scatter ratio should be >1 for A280.

  df_5 <- c()

  for(diln in unique(df_4$dilution)){

    # Subset
    temp_diln_values <- df_4 %>%
      dplyr::filter(.data$dilution == diln)

    # max
    abs_280 <- temp_diln_values %>%
      dplyr::filter(.data$measure == 280) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_280
    abs_333 <- temp_diln_values %>%
      dplyr::filter(.data$measure == wav_to_use2) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_333

    # Add corrected value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(wav_corr2 = wav_to_use2) %>% # record wav_to_use2
      dplyr::mutate(abs280_corr2 = abs_280 - scatter_ratio*abs_333) # normalise
    temp_diln_values

    # Get conc from this value and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_a280_corr2 = .data$abs280_corr2/.data$EC_A280_mgml)
    temp_diln_values

    # Rbind
    df_5 <- rbind(df_5, temp_diln_values)

  }

  df_5

  # Subset data to get rid of negatives:
  df_5_subset <- subset(df_5, df_5$conc_a280_corr2 >= 0)
  # Fit model as (Y ~ X):
  model3 <- stats::lm(conc_a280_corr2 ~ dilution + 0, data = df_5_subset) # force through zero to compare points to ideal
  model3

  # Plot
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_5,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_corr2*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_5_subset,
                        ggplot2::aes(x = .data$dilution, y = .data$conc_a280_corr2*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_5, # to extend line to zero
                       ggplot2::aes(x = .data$dilution,
                                    y = (stats::predict(model3, df_5))*1000),
                       colour = "red") +
    ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
    ggplot2::scale_y_continuous("concentration (ng/ul)") +
    ggplot2::labs(subtitle = "predicted conc (correction = scatter)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1
  plotname <- "plot5d_a280_scatternorm.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  ##

  # 5d. Compare methods

  # Redundant unless fits fail.

  # # Plot
  # df_6 <- df_5 %>%
  #   # dplyr::ungroup() %>%
  #   dplyr::select(calibrant, protein, dilution,
  #                 measure,
  #                 conc_a280_std, conc_a280_corr1, conc_a280_corr2) %>%
  #   dplyr::filter(measure == 280)
  # df_6 <- df_6 %>%
  #   tidyr::pivot_longer(cols = c(
  #     conc_a280_std, conc_a280_corr1, conc_a280_corr2),
  #     names_to = "type", values_to = "value")
  # df_6
  #
  # df_6$type <- factor(df_6$type, levels = c("conc_a280_std", "conc_a280_corr1", "conc_a280_corr2"))
  # type.labs <- c("standard", "baseline", "scatter")
  # names(type.labs) <- c("conc_a280_std", "conc_a280_corr1", "conc_a280_corr2")
  # plot1 <- ggplot2::ggplot(data = df_6) +
  #
  #   ggplot2::geom_point(ggplot2::aes(x = dilution, y = value*1000)) +
  #   ggplot2::geom_hline(yintercept = 0) +
  #
  #   ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
  #   ggplot2::scale_y_continuous("concentration (ng/ul)") +
  #   ggplot2::facet_grid(. ~ type,
  #                       labeller = ggplot2::labeller(type = type.labs)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(
  #     aspect.ratio = 1,
  #     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
  #     panel.grid.minor = ggplot2::element_blank(),
  #     strip.background = ggplot2::element_blank(),
  #     strip.text = ggplot2::element_text(face = "bold", hjust = 0)
  #   )
  # plot1
  # plotname <- "plot5e1_a280_all.pdf"
  # ggplot2::ggsave(file.path(outfolder, plotname),
  #                 plot = plot1,
  #                 width = 12, height = 12, units = "cm")
  #
  # plot1 <- ggplot2::ggplot(data = df_6 %>%
  #                            dplyr::filter(dilution != 0)) +
  #   ggplot2::geom_point(ggplot2::aes(x = dilution, y = value*1000)) +
  #   ggplot2::scale_x_log10("dilution") +
  #   ggplot2::scale_y_log10("concentration (ng/ul)") +
  #   ggplot2::facet_grid(. ~ type,
  #                       labeller = ggplot2::labeller(type = type.labs)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(
  #     aspect.ratio = 1,
  #     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
  #     panel.grid.minor = ggplot2::element_blank(),
  #     strip.background = ggplot2::element_blank(),
  #     strip.text = ggplot2::element_text(face = "bold", hjust = 0),
  #   )
  # plot1
  # plotname <- "plot5e2_a280_all_logplot.pdf"
  # ggplot2::ggsave(file.path(outfolder, plotname),
  #                 plot = plot1,
  #                 width = 12, height = 12, units = "cm")

  ##

  # 6. Fit models (conc ~ dilution) ----------------------------------------------------------------

  ## Tidy data
  df_6 <- df_5 %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$media, .data$calibrant, .data$protein, .data$dilution,
                  .data$measure,
                  .data$conc_a280_std, .data$conc_a280_corr1, .data$conc_a280_corr2) %>%
    dplyr::filter(.data$measure == 280)
  df_6 <- df_6 %>%
    tidyr::pivot_longer(cols = c(
      .data$conc_a280_std, .data$conc_a280_corr1, .data$conc_a280_corr2),
      names_to = "type", values_to = "value")
  df_6
  df_6$type <- factor(df_6$type, levels = c("conc_a280_std", "conc_a280_corr1", "conc_a280_corr2"))

  ## Fit model for conc ~ dilution ---------

  # remove blanks
  df_7 <- df_6 %>%
    dplyr::filter(.data$protein != "none")

  # nest to get to df with three rows, one per model required
  # nesting methods from https://r4ds.had.co.nz/many-models.html
  df_7_nested <- df_7 %>%
    dplyr::group_by(.data$media, .data$calibrant, .data$protein, .data$measure, .data$type) %>%
    tidyr::nest() # basically nesting dilution and value cols only
  df_7_nested
  df_7_nested$data

  # purrr map to make models
  df_7_nested <- df_7_nested %>%
    dplyr::mutate(model = purrr::map(.data$data, stats::lm, formula = value ~ dilution + 0))
  df_7_nested
  df_7_nested$model

  # extract fitted values
  df_7_nested <- df_7_nested %>%
    dplyr::mutate(tidied = purrr::map(.data$model, broom::tidy)) %>% # tidy fn from broom # https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    tidyr::unnest(.data$tidied)
  df_7_nested

  # save model coefficients
  df_7_model_tosave <- df_7_nested %>%
    dplyr::select(-.data$model, -.data$data)
  df_7_model_tosave
  csvname <- "a280_coeffs.csv"
  utils::write.csv(x = df_7_model_tosave, file = file.path(outfolder, csvname), row.names = FALSE)

  ### Plot points and model fits

  df_9 <- df_7_nested %>%
    tidyr::unnest(.data$data) %>%
    dplyr::mutate(pred_conc = .data$estimate*.data$dilution)
  df_9
  type.labs <- c("none", "baseline", "scatter")
  names(type.labs) <- c("conc_a280_std", "conc_a280_corr1", "conc_a280_corr2")

  plot1 <- ggplot2::ggplot() +
    # data
    ggplot2::geom_point(data = df_6, ggplot2::aes(x = .data$dilution, y = .data$value*1000)) +
    # model
    ggplot2::geom_line(data = df_9, ggplot2::aes(x = .data$dilution, y = .data$pred_conc*1000)) +
    ggplot2::geom_hline(yintercept = 0) +

    ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
    ggplot2::scale_y_continuous("concentration (ng/ul)") +
    ggplot2::facet_grid(. ~ type,
                        labeller = ggplot2::labeller(type = type.labs)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0)
    )
  plot1
  plotname <- "plot6a_a280_models_all.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  plot1 <- ggplot2::ggplot() +
    # data
    ggplot2::geom_point(data = df_6 %>%
                          dplyr::filter(.data$dilution != 0),
                        ggplot2::aes(x = .data$dilution, y = .data$value*1000)) +
    # model
    ggplot2::geom_line(data = df_9 %>%
                         dplyr::filter(.data$dilution != 0),
                       ggplot2::aes(x = .data$dilution, y = .data$pred_conc*1000)) +
    ggplot2::scale_x_log10("dilution") +
    ggplot2::scale_y_log10("concentration (ng/ul)") +
    ggplot2::facet_grid(. ~ type,
                        labeller = ggplot2::labeller(type = type.labs)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    )
  plot1
  plotname <- "plot6b_a280_models_all_logplot.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  # 7. Export concentration predictions ---------------

  ### Tidy dfs for export
  # NB Need to use a base df that contains columns like well, row, column - use input df

  ## Prep base df
  names(spectrum_data)
  spectrum_data_tidy <- spectrum_data %>%
    dplyr::select(c( .data$media:.data$well )) %>% # take all cols between A:B - https://suzan.rbind.io/2018/01/dplyr-tutorial-1/
    dplyr::select(c( -.data$volume )) # remove volume as downstream assays may have used different volumes
  names(spectrum_data_tidy)
  # loads of duplicated rows here now that measure (wavelength) and absorption value columns have been deleted
  nrow(spectrum_data_tidy) # eg. 19,224
  spectrum_data_tidy <- spectrum_data_tidy %>%
    dplyr::distinct()
  nrow(spectrum_data_tidy) # eg. 24 (so there were 801 wavelengths)
  spectrum_data_tidy

  ## Build new df
  conc_data <- spectrum_data_tidy
  conc_data

  # complete MW column
  protein_seq
  conc_data$mw_gmol1 <- fpcountr::get_mw(protein = protein_seq)

  # Pick correction method
  df_7_nested
  if(corr_method == "none"){
    conc_fit_coef <- df_7_nested %>%
      dplyr::filter(.data$type == "conc_a280_std") %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$estimate) %>%
      as.numeric()
    conc_fit_coef
  }
  if(corr_method == "baseline"){
    conc_fit_coef <- df_7_nested %>%
      dplyr::filter(.data$type == "conc_a280_corr1") %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$estimate) %>%
      as.numeric()
    conc_fit_coef
  }
  if(corr_method == "scatter"){
    conc_fit_coef <- df_7_nested %>%
      dplyr::filter(.data$type == "conc_a280_corr2") %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$estimate) %>%
      as.numeric()
    conc_fit_coef
  }

  # add concs
  conc_data <- conc_data %>%
    # Concentration = predicted concentration for the highest dilution, multiplied down for each dilution
    dplyr::mutate(concentration_ngul = conc_fit_coef * .data$dilution * 1000) # *1000 to convert from ug/ul to ng/ul
  conc_data

  # fill in protein column where dilution exists (for ease of future joining)
  conc_data <- conc_data %>%
    dplyr::mutate(protein = ifelse(.data$protein == "" & !is.na(.data$dilution), .data$calibr, .data$protein))
  # if protein column is empty AND dilution is not NA, add the calibr into the protein column, else, leave protein entry as it is
  conc_data

  ##

  # Save and Return -----------

  # Everything
  csvname <- gsub(".csv", paste0("_a280.csv"), basename(processed_spectrum_csv))
  utils::write.csv(x = df_5, file = file.path(outfolder, csvname), row.names = FALSE)

  # Protein concentrations
  csvname <- "protein_concs_a280.csv"
  utils::write.csv(x = conc_data, file = file.path(outfolder, csvname), row.names = FALSE)

  return(conc_data)

}
