#' Get FP concentrations using ECmax method
#'
#' Get protein's concentration from a *dilution series* measured with an
#' *A200-1000 absorbance spectrum*. Based on `get_conc_A280`, but uses "ECmax",
#' the FPbase-stated extinction coefficient for the FPbase-stated maximal
#' excitation wavelength (which usually corresponds to its maximal absorbance
#' wavelength). ... Function expects an input such as the csv file exported from
#' `plot_absorbance_spectrum` called '_processed.csv', which contains values
#' corrected for path length and normalised to blanks as a column called
#' `normalised_cm1_value`, but retains replicate data containing positional
#' (well) information required for exporting predicted concentrations at the end
#' of this function. ... Function uses `fpcountR::get_properties` to get FPbase
#' EC in M-1cm-1 and wavelength, and converts it to an ECmax mass extinction
#' coefficient in (mgml)-1cm-1 using the MW (worked out from `protein_seq` and
#' `fpcountR::get_mw`). Then the function uses the EC_max_mgml to work out the
#' concentration of protein in each well, using three correction methods.
#' Instead of using the normalised data directly, the values used are based on a
#' LOESS fit through the absorption spectra to minimise fluctuations due to
#' noise. ... Finally, linear models are fitted to each concentration prediction
#' method, and a csv file is built containing predicted concentrations according
#' to the user's chosen correction method. Plots showing each of the analytical
#' steps are saved concurrently. Troubleshooting: for 'incompatible lengths'
#' errors, adjust xrange to avoid noisy wavelengths.
#'
#' @param protein_slug character string of protein name in 'slug' form to match
#'   slug of FPbase entry.
#' @param protein_seq character string of protein sequence using 1-letter code.
#'   Required for MW calculation.
#' @param processed_spectrum_csv Path to csv file of a processed absorbance
#'   spectrum. Processing should be done with `plot_absorbance_spectrum3`, which
#'   corrects for path lengths and normalises to blank wells.
#' @param wells_to_remove list of wells to remove before analysis. Defaults to
#'   NULL.
#' @param xrange list of two numerical values corresponding to the wavelength
#'   range to keep when fitting the loess model across the absorbance spectrum.
#'   By default these values are 250nm and 800nm but where the data at the UV
#'   range is noisy, adjusting the xrange can prevent errors in the fitting.
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
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#' @param filename filename of csv file from fpbase
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @examples
#' ecmax_concs <- get_conc_ECmax(protein_slug = "mcherry", protein_seq = protein_seq, processed_spectrum_csv = "abs_parsed_processed.csv", corr_method = "scatter", wav_to_use1 = 700, wav_to_use2 = 315, outfolder = "protquant_ecmax/mCherry_T5N15pi", filename = "mCherry_properties.csv")

get_conc_ECmax <- function(protein_slug, protein_seq,
                           processed_spectrum_csv, wells_to_remove = NULL,
                           xrange = c(250,800),
                           corr_method = "none", # "none", "baseline", "scatter"
                           wav_to_use1 = 340, wav_to_use2 = 333,
                           outfolder, filename
){

  # Get data -------------------------------------------------

  spectrum_data <- utils::read.csv(processed_spectrum_csv, header = TRUE, check.names = FALSE)

  # Where protein is none, dilution is usually NA, so these rows are lost in later steps.
  # To keep them, they need a value.
  spectrum_data <- spectrum_data %>%
    dplyr::mutate(dilution = ifelse(protein == "none", 0, dilution))
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

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # 1. Replicates data -------------------------------------------------

  # Check input data

  # Plot absorbance spectra (250-800nm) of replicates separately as sanity check for reproducibility
  data.to.plot <- spectrum_data_subset %>%
    dplyr::filter(measure > 250 & measure < 800)
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = measure, y = normalised_cm1_value,
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
    dplyr::filter(measure > 250 & measure < 800)
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order

  # Get ECmax wavelength of your FP for plotting
  fp_properties <- fpcountr::get_properties(slug = protein_slug, verbose = TRUE, outfolder = outfolder, filename = filename)
  fp_properties
  # fp_properties$ex_max # is excitation max
  # fp_properties$ext_coeff # is EC

  # If fp properties are missing, stop here
  if(nrow(fp_properties) == 0){
    message("Error: The FP properties table for ", protein_slug, " is empty.")
    print(fp_properties)
    message("Stopping..")
    return()
  }
  # If key fp properties are missing, throw an error here (stop later)
  if(is.na(fp_properties$ex_max) | is.na(fp_properties$ext_coeff)){
    message("Error: The FP properties table for ", protein_slug, " does not have all the properties required for calculating protein concentration via the ECmax method.\nBoth excitation maximum (ex_max) and extinction coefficient (ext_coeff) are required.")
    print(fp_properties)
    # if fp_properties$ext_coeff = NA: section 4 will lead to EC_max_mgml = NA, and section 5 will fail.
    # if fp_properties$ex_max = NA: section 3 plots ok with as.numeric(), then will fail at section 5.
    # stop both after section 3.
  }

  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = measure, y = normalised_cm1_value),
                        colour = "lightblue") +
    ggplot2::geom_smooth(ggplot2::aes(x = measure, y = normalised_cm1_value),
                         colour = "black",
                         span = 0.5/5.5 ### could leave as default, us 0.5, or add as argument in function # 0.5 used for 250-350
    ) +
    ggplot2::geom_vline(xintercept = as.numeric(fp_properties$ex_max), colour = "red") + # as.numeric to handle missing (NA) values without error

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

  ##

  ## Fit loess model for absorbance ~ wavelength (separately for each dilution)

  # method slightly complex. using both:
  # https://r4ds.had.co.nz/many-models.html
  # https://stackoverflow.com/questions/50163106/loess-regression-on-each-group-with-dplyrgroup-by

  # subset
  spectrum_data_subset <- summ_data %>%
    dplyr::filter(measure > xrange[1]) %>% # default 250
    dplyr::filter(measure < xrange[2]) # default 800
  # nest
  spectrum_data_nested <- spectrum_data_subset %>%
    dplyr::group_by(instrument, media, calibrant, dilution) %>%
    tidyr::nest()
  spectrum_data_nested
  spectrum_data_nested$data
  # purrr map to make model
  spectrum_data_nested <- spectrum_data_nested %>%
    dplyr::mutate(model = purrr::map(data, stats::loess, formula = normalised_cm1_value ~ measure, span = 0.5/5.5)) # 0.5 for 250-350, 0.5/5.5 for 250-800
  spectrum_data_nested
  # extract fitted values!!
  spectrum_data_nested <- spectrum_data_nested %>%
    dplyr::mutate(fitted_cm1_value = purrr::map(model, `[[`, "fitted"))
  spectrum_data_nested
  # remove model column (awkward and no longer necessary) and unnest
  spectrum_data_model <- spectrum_data_nested %>%
    dplyr::select(-model) %>%
    tidyr::unnest(cols = c(data, fitted_cm1_value))
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
    ggplot2::geom_point(ggplot2::aes(x = measure, y = normalised_cm1_value), colour = "lightblue") +
    # geom_line of extracted loess fitted values
    ggplot2::geom_line(ggplot2::aes(x = measure, y = fitted_cm1_value)) +
    ggplot2::geom_vline(xintercept = as.numeric(fp_properties$ex_max), colour = "red") + # as.numeric to handle missing (NA) values without error
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

  ##

  # If key properties are missing, stop here
  if(is.na(fp_properties$ex_max) | is.na(fp_properties$ext_coeff)){
    message("Stopping..")
    return()
  }

  # 4. Get extinction coefficient ------------------------------------------------------------

  ## 4a. Get protein sequence
  protein_seq

  ## 4b. Get MW
  protein_mw <- fpcountr::get_mw(protein = protein_seq)
  protein_mw

  ## 4c. Get mgml extinction coefficient (get_properties fn call moved up to allow indication of ECmax in previous plots)
  # fp_properties <- fpcountr::get_properties(slug = protein_slug, verbose = TRUE, outfolder = outfolder, filename = filename)
  # fp_properties
  # fp_properties$ex_max # is excitation max
  # fp_properties$ext_coeff # is EC

  ## 4d. Add 'ECmax of 1mg/ml of protein' to table
  # EC of 0.1% solution = 0.1g/100ml = 0.001g/ml = 1mg/ml
  # EC/MW = (M-1cm-1)/(g/mol) = (mol/L)-1 * (cm)-1 / (g/mol) = (L/mol)*(1/cm)*(g/mol) = L/g * 1/cm = (g/L)-1(cm)-1 = (mg/ml)-1(cm)-1
  # = EC (mg/ml)-1(cm)-1

  df_2 <- spectrum_data_model %>%
    dplyr::mutate(EC_max_mgml = fp_properties$ext_coeff / protein_mw)
  head(df_2)

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
    abs_max <- temp_diln_values %>%
      dplyr::filter(measure == fp_properties$ex_max) %>%
      dplyr::ungroup() %>%
      dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_max

    # Add this Abs value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(absmax_std = abs_max)
    temp_diln_values

    # Get conc from this Abs and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_max_std = absmax_std/EC_max_mgml)
    temp_diln_values

    # Rbind
    df_3 <- rbind(df_3, temp_diln_values)

  }

  df_3

  # Plot processing steps
  plot1 <- ggplot2::ggplot(data = df_3 %>%
                             dplyr::filter(measure == fp_properties$ex_max)) +
    ggplot2::geom_point(ggplot2::aes(x = dilution, y = raw_value), colour = "grey") +
    ggplot2::geom_point(ggplot2::aes(x = dilution, y = raw_cm1_value), colour = "black") +
    ggplot2::geom_point(ggplot2::aes(x = dilution, y = normalised_cm1_value), colour = "red") +
    ggplot2::geom_point(ggplot2::aes(x = dilution, y = fitted_cm1_value), colour = "blue") +
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
  plotname <- "plot5a_ecmax.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  # Subset data to get rid of negatives:
  df_3_subset <- subset(df_3, conc_max_std >= 0)
  # Fit model as (Y ~ X):
  model1 <- stats::lm(conc_max_std ~ dilution + 0, data = df_3_subset) # force through zero to compare points to ideal
  model1
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_3,
                        ggplot2::aes(x = dilution, y = conc_max_std*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_3_subset,
                        ggplot2::aes(x = dilution, y = conc_max_std*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_3, # to extend line to zero
                       ggplot2::aes(x = dilution,
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
  plotname <- "plot5b_ecmax_stdmethod.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  ##

  # 5b. Correction method1 - baseline ------------------------------------------------------------------------

  # Correction method1 to remove scatter
  # Baseline aka Nanodrop method: normalise to A340
  # So use Abs = abs_ecmax - abs_340

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
    abs_max <- temp_diln_values %>%
      dplyr::filter(measure == fp_properties$ex_max) %>%
      dplyr::ungroup() %>%
      dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_max
    abs_340 <- temp_diln_values %>%
      dplyr::filter(measure == wav_to_use1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_340

    # Add corrected value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(wav_corr1 = wav_to_use1) %>% # record wav_to_use1
      dplyr::mutate(absmax_corr1 = abs_max - abs_340) # normalise
    temp_diln_values

    # Get conc from this Abs and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_max_corr1 = absmax_corr1/EC_max_mgml)
    temp_diln_values

    # Rbind
    df_4 <- rbind(df_4, temp_diln_values)

  }

  df_4

  # Baseline fit check, entire spectrum:
  data.to.plot <- df_4
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  # use top dilution
  data.to.plot <- subset(data.to.plot, dilution == levels(data.to.plot$dilution)[1])
  # fit baseline to data at chosen wavelength
  baselinefit <- data.to.plot %>%
    dplyr::filter(measure == wav_to_use1) %>%
    dplyr::ungroup() %>%
    dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
    as.numeric()
  baselinefit
  #
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = measure, y = normalised_cm1_value), colour = "lightblue") +
    # geom_line of extracted loess fitted values
    ggplot2::geom_line(ggplot2::aes(x = measure, y = fitted_cm1_value)) +
    ggplot2::geom_vline(xintercept = fp_properties$ex_max, colour = "red") +
    # baseline fit
    ggplot2::geom_hline(yintercept = baselinefit, colour = "blue") +

    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(xrange[1],xrange[2])) +
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
  plotname <- "plot5c_ecmax_baselinenorm_baselinecheck.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  # Subset data to get rid of negatives:
  df_4_subset <- subset(df_4, conc_max_corr1 >= 0)
  # Fit model as (Y ~ X):
  model2 <- stats::lm(conc_max_corr1 ~ dilution + 0, data = df_4_subset) # force through zero to compare points to ideal
  model2
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_4,
                        ggplot2::aes(x = dilution, y = conc_max_corr1*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_4_subset,
                        ggplot2::aes(x = dilution, y = conc_max_corr1*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_4, # to extend line to zero
                       ggplot2::aes(x = dilution,
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
  plotname <- "plot5c_ecmax_baselinenorm.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  ##

  # 5c. Correction method2 - scatter ------------------------------------------------------------

  # Correction method2 to remove scatter
  # Scatter method:
  # For A280 measurements, A333 gives 1/2 the scatter as A280, so use Abs = abs_280 - 2*abs_333
  # For ECmax measurements, calculate scatter_ratio between A333 and ECmax wavelength, then use Abs = abs_ecmax - scatter_ratio*abs_333
  # Edit: wavelength used here defined by wav_to_use2, with default 333nm.

  # ## Correcting for ECmax - how and why?
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
  # # light scatter at ecmax (eg 500nm) is not same level as light scatter elsewhere
  # scatter_a500 <- eq(500)
  # scatter_a500 # 1.626926e-10
  # # to find light scatter at ecmax when there's a peak in the way,
  # # need to find it by looking for light scatter at a different wavelength, then doing a transformation
  # # use A333 or a wavelength where FP absorbs minimally, and find ratio
  # # eg. A333
  # scatter_a333 <- eq(333)
  # scatter_ecmax <- eq(587) # for eg
  # scatter_ratio <- scatter_ecmax/scatter_a333
  # scatter_ratio
  # # so scatter at ecmax should be ~0.1x scatter at a333...

  # Work out scatter_ratio.
  # What is ratio between scatter at ECmax wavelength and scatter at the user-defined wav_to_use2 wavelength?
  eq = function(x){1/x^4}
  scatter_a333 <- eq(wav_to_use2) # 333 by default
  scatter_a333
  scatter_ecmax <- eq(fp_properties$ex_max)
  scatter_ecmax
  scatter_ratio <- scatter_ecmax/scatter_a333
  scatter_ratio # the scatter ratio should be <1.

  df_5 <- c()

  for(diln in unique(df_4$dilution)){

    # Subset
    temp_diln_values <- df_4 %>%
      dplyr::filter(.data$dilution == diln)

    # max
    abs_max <- temp_diln_values %>%
      dplyr::filter(measure == fp_properties$ex_max) %>%
      dplyr::ungroup() %>%
      dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_max
    abs_333 <- temp_diln_values %>%
      dplyr::filter(measure == wav_to_use2) %>%
      dplyr::ungroup() %>%
      dplyr::select(fitted_cm1_value) %>% # so pathlength = 1cm
      as.numeric()
    abs_333

    # Add corrected value as column
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(wav_corr2 = wav_to_use2) %>% # record wav_to_use2
      dplyr::mutate(abs_corr2 = abs_333) %>% # record abs of wav_to_use2 (needed for plot5d scatter check)
      dplyr::mutate(absmax_corr2 = abs_max - scatter_ratio*abs_333) # normalise
    temp_diln_values

    # Get conc from this value and add it into df
    temp_diln_values <- temp_diln_values %>%
      dplyr::mutate(conc_max_corr2 = absmax_corr2/EC_max_mgml)
    temp_diln_values

    # Rbind
    df_5 <- rbind(df_5, temp_diln_values)

  }

  df_5

  # Scatter fit check, entire spectrum:
  data.to.plot <- df_5
  data.to.plot$dilution <- as.factor(data.to.plot$dilution) # make dilution a factor
  newlist <- levels(data.to.plot$dilution)
  data.to.plot$dilution <- factor(data.to.plot$dilution, levels = rev(newlist)) # reverse the order
  # use only one dilution (top dilution), as stat_function doesn't play well with faceted plots
  data.to.plot <- subset(data.to.plot, dilution == levels(data.to.plot$dilution)[1])
  ## fit scatter eqn to data at chosen wavelength
  eq = function(x){1/x^4}
  # eq(wav_to_use2) # value at scatter wavelength for 1/x^4 eqn
  # unique(data.to.plot$abs_corr2) # abs at scatter wavelength
  # coefficient to supply to coefficient*(1/x^4) eqn to fit scatter curve to data:
  scatterfit_coefficient <- unique(data.to.plot$abs_corr2)/eq(wav_to_use2)
  scatterfit_coefficient
  #
  plot1 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # bottom
    ggplot2::geom_point(ggplot2::aes(x = measure, y = normalised_cm1_value), colour = "lightblue") +
    # geom_line of extracted loess fitted values
    ggplot2::geom_line(ggplot2::aes(x = measure, y = fitted_cm1_value)) +
    ggplot2::geom_vline(xintercept = fp_properties$ex_max, colour = "red") +
    # scatter fit
    ggplot2::stat_function(fun = function(x) scatterfit_coefficient/x^4, colour = "blue") +

    ggplot2::scale_x_continuous("wavelength (nm)", limits = c(xrange[1],xrange[2])) +
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
  plotname <- "plot5d_ecmax_scatternorm_scattercheck.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 8, height = 8, units = "cm")

  # Subset data to get rid of negatives:
  df_5_subset <- subset(df_5, conc_max_corr2 >= 0)
  # Fit model as (Y ~ X):
  model3 <- stats::lm(conc_max_corr2 ~ dilution + 0, data = df_5_subset) # force through zero to compare points to ideal
  model3

  # Plot
  plot1 <- ggplot2::ggplot() +
    # all points
    ggplot2::geom_point(data = df_5,
                        ggplot2::aes(x = dilution, y = conc_max_corr2*1000)) +
    # points used in model
    ggplot2::geom_point(data = df_5_subset,
                        ggplot2::aes(x = dilution, y = conc_max_corr2*1000),
                        colour = "red") +
    ggplot2::geom_line(data = df_5, # to extend line to zero
                       ggplot2::aes(x = dilution,
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
  plotname <- "plot5d_ecmax_scatternorm.pdf"
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
  #                 conc_max_std, conc_max_corr1, conc_max_corr2) %>%
  #   dplyr::filter(measure == fp_properties$ex_max)
  # df_6 <- df_6 %>%
  #   tidyr::pivot_longer(cols = c(
  #     conc_max_std, conc_max_corr1, conc_max_corr2),
  #     names_to = "type", values_to = "value")
  # df_6
  #
  # df_6$type <- factor(df_6$type, levels = c("conc_max_std", "conc_max_corr1", "conc_max_corr2"))
  # type.labs <- c("standard", "baseline", "scatter")
  # names(type.labs) <- c("conc_max_std", "conc_max_corr1", "conc_max_corr2")
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
  # plotname <- "plot5e1_ecmax_all.pdf"
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
  # plotname <- "plot5e2_ecmax_all_logplot.pdf"
  # ggplot2::ggsave(file.path(outfolder, plotname),
  #                 plot = plot1,
  #                 width = 12, height = 12, units = "cm")

  ##

  # 6. Fit models (conc ~ dilution) ----------------------------------------------------------------

  ## Tidy data
  df_6 <- df_5 %>%
    dplyr::ungroup() %>%
    dplyr::select(media, calibrant, protein, dilution,
                  measure,
                  conc_max_std, conc_max_corr1, conc_max_corr2) %>%
    dplyr::filter(measure == fp_properties$ex_max)
  df_6 <- df_6 %>%
    tidyr::pivot_longer(cols = c(
      conc_max_std, conc_max_corr1, conc_max_corr2),
      names_to = "type", values_to = "value")
  df_6
  df_6$type <- factor(df_6$type, levels = c("conc_max_std", "conc_max_corr1", "conc_max_corr2"))

  ## Fit model for conc ~ dilution ---------

  # remove blanks
  df_7 <- df_6 %>%
    dplyr::filter(protein != "none")

  # nest to get to df with three rows, one per model required
  # nesting methods from https://r4ds.had.co.nz/many-models.html
  df_7_nested <- df_7 %>%
    dplyr::group_by(media, calibrant, protein, measure, type) %>%
    tidyr::nest() # basically nesting dilution and value cols only
  df_7_nested
  df_7_nested$data

  # purrr map to make models
  df_7_nested <- df_7_nested %>%
    dplyr::mutate(model = purrr::map(data, stats::lm, formula = value ~ dilution + 0))
  df_7_nested
  df_7_nested$model

  # extract fitted values
  df_7_nested <- df_7_nested %>%
    dplyr::mutate(tidied = purrr::map(model, broom::tidy)) %>% # tidy fn from broom # https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    tidyr::unnest(tidied)
  df_7_nested

  # save model coefficients
  df_7_model_tosave <- df_7_nested %>%
    dplyr::select(-model, -data)
  df_7_model_tosave
  csvname <- "ecmax_coeffs.csv"
  utils::write.csv(x = df_7_model_tosave, file = file.path(outfolder, csvname), row.names = FALSE)

  ### Plot points and model fits

  df_9 <- df_7_nested %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(pred_conc = estimate*dilution)
  df_9
  type.labs <- c("none", "baseline", "scatter")
  names(type.labs) <- c("conc_max_std", "conc_max_corr1", "conc_max_corr2")

  plot1 <- ggplot2::ggplot() +
    # data
    ggplot2::geom_point(data = df_6, ggplot2::aes(x = dilution, y = value*1000)) +
    # model
    ggplot2::geom_line(data = df_9, ggplot2::aes(x = dilution, y = pred_conc*1000)) +
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
  plotname <- "plot6a_ecmax_models_all.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  plot1 <- ggplot2::ggplot() +
    # data
    ggplot2::geom_point(data = df_6 %>%
                          dplyr::filter(dilution != 0),
                        ggplot2::aes(x = dilution, y = value*1000)) +
    # model
    ggplot2::geom_line(data = df_9 %>%
                          dplyr::filter(dilution != 0),
                        ggplot2::aes(x = dilution, y = pred_conc*1000)) +
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
  plotname <- "plot6b_ecmax_models_all_logplot.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot1,
                  width = 12, height = 12, units = "cm")

  # 7. Export concentration predictions ---------------

  ### Tidy dfs for export
  # NB Need to use a base df that contains columns like well, row, column - use input df

  ## Prep base df
  names(spectrum_data)
  spectrum_data_tidy <- spectrum_data %>%
    dplyr::select(c( media:well )) %>% # take all cols between A:B - https://suzan.rbind.io/2018/01/dplyr-tutorial-1/
    dplyr::select(c( -volume )) # remove volume as downstream assays may have used different volumes
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
      dplyr::filter(type == "conc_max_std") %>%
      dplyr::ungroup() %>%
      dplyr::select(estimate) %>%
      as.numeric()
    conc_fit_coef
  }
  if(corr_method == "baseline"){
    conc_fit_coef <- df_7_nested %>%
      dplyr::filter(type == "conc_max_corr1") %>%
      dplyr::ungroup() %>%
      dplyr::select(estimate) %>%
      as.numeric()
    conc_fit_coef
  }
  if(corr_method == "scatter"){
    conc_fit_coef <- df_7_nested %>%
      dplyr::filter(type == "conc_max_corr2") %>%
      dplyr::ungroup() %>%
      dplyr::select(estimate) %>%
      as.numeric()
    conc_fit_coef
  }

  # add concs
  conc_data <- conc_data %>%
    # Concentration = predicted concentration for the highest dilution, multiplied down for each dilution
    dplyr::mutate(concentration_ngul = conc_fit_coef * dilution * 1000) # *1000 to convert from ug/ul to ng/ul
  conc_data

  # fill in protein column where dilution exists (for ease of future joining)
  conc_data <- conc_data %>%
    dplyr::mutate(protein = ifelse(.data$protein == "" & !is.na(.data$dilution), calibr, protein))
  # if protein column is empty AND dilution is not NA, add the calibr into the protein column, else, leave protein entry as it is
  conc_data

  ##

  # Save and Return -----------

  # Everything
  csvname <- gsub(".csv", paste0("_ecmax.csv"), basename(processed_spectrum_csv))
  utils::write.csv(x = df_5, file = file.path(outfolder, csvname), row.names = FALSE)

  # Protein concentrations
  csvname <- "protein_concs_ecmax.csv"
  utils::write.csv(x = conc_data, file = file.path(outfolder, csvname), row.names = FALSE)

  return(conc_data)

}

