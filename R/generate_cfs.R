#' Generate conversion factors for FP calibrations
#'
#' Generate conversion factors for fluorescent protein (FP) calibrations.
#' Originally based on `flopr::generate_cfs()` but with numerous changes. The
#' original function was intended for fluorescein and microsphere calibrations,
#' whereas this one can be used for any fluorescent calibrant, including
#' proteins. Takes as input a parsed CSV of fluorescence data of a dilution
#' series of FPs at one or more gains, that contains both data and metadata
#' columns (including: `instrument`, `plate`, `seal`, `channel_name`,
#' `channel_ex`, `channel_em`, `media`, `calibrant`, `protein`, `replicate`,
#' `volume`, `mw_gmol1`, `concentration_ngul`, `well`, the columns for the
#' fluorescence data, `row`, `column`). A number of arguments allow the tweaking
#' of the original data set. Following this, the data is reshaped, normalised,
#' trimmed of saturated points, summarised and used to fit a model for the
#' conversion factors from arbitrary to absolute units. Optional extras: the
#' `sensitivity_plots` argument extends this analysis to identify the limits of
#' detection and relative sensitivity and dynamic range of each gain. Plots are
#' saved to record the processed data at every step, allowing for visual sanity
#' checks and troubleshooting. A CSV file with the fitted conversion factors is
#' saved (along with the processed data if requested with `more_csvs`).
#'
#' @param calibration_csv character string. Path of the calibration data CSV
#'   file.
#' @param more_csvs logical. Optionally save all the intermediate tables in this
#'   function. Defaults to FALSE.
#' @param more_plots logical. Optionally save more plots. Defaults to FALSE.
#' @param sensitivity_plots logical. Optionally adds a few extra columns to the
#'   data such as the min/max normalised fluorescence and detectable molecules.
#'   It also plots a few extra plots. Defaults to FALSE.
#' @param include_only character string. If specified, only includes the
#'   measures (column names) specified here.
#' @param exclude character string. If specified, excludes any measures (column
#'   names) specified here.
#' @param gain_fix logical. Optionally add "0" before 2-digit gain names i.e.
#'   'GFP 40' -> 'GFP 040', or `blueblue 40` -> `blueblue 040`- which fixes
#'   ordering of plots. Defaults to FALSE.
#' @param rename_from character string. Rename measures (column names)
#'   containing character string `rename_from` to character string specified in
#'   `rename_to`. Both `rename_from` and `rename_to` need to be completed to
#'   trigger renaming.
#' @param rename_to character string. Rename measures (column names) containing
#'   character string `rename_from` to character string specified in
#'   `rename_to`. Both `rename_from` and `rename_to` need to be completed to
#'   trigger renaming.
#' @param subset_rows logical. Should script take a subset of the rows (or whole
#'   table)? Defaults to FALSE.
#' @param rows_to_keep character array. If `subset_rows` is TRUE, script will
#'   choose rows to keep from this list. Defaults to `c("C","D")`.
#' @param separator character string that represents the separator between the
#'   channel name and the gain value in the measures columns, e.g. for "GFP 40"
#'   it is " ", for "GFP40" it is "" and for "GFP_40" it is "_". Required for
#'   plotting the gain vs conversion factors. Defaults to "".
#' @param complete_blank logical. Optionally adds "0" to the
#'   `concentration_ngul` column for wells identified as blank (`protein` =
#'   "none"). Useful if metadata is missing these values. Defaults to FALSE.
#' @param outfolder character string. Path to folder where output files should
#'   be saved.
#'
#' @export
#' @return data frame of conversion factors
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#'
#' @examples
#' \dontrun{
#'   fitvals <- generate_cfs(
#'     calibration_csv = "data_parsed.csv",
#'     subset_rows = TRUE, rows_to_keep = c("C","D"),
#'     outfolder = "cfs_mCherry"
#'   )
#' }
generate_cfs <- function(calibration_csv,
                         more_csvs = FALSE,
                         more_plots = FALSE,
                         sensitivity_plots = FALSE,
                         include_only = NULL, exclude = NULL,
                         gain_fix = FALSE,
                         rename_from = NULL, rename_to = NULL,
                         subset_rows = FALSE, rows_to_keep = c("C","D"),
                         separator = "",
                         complete_blank = FALSE,
                         outfolder = "") {

  # 0. Get data and prep data -------------------------------------------------

  calibration_data <- utils::read.csv(calibration_csv, header = TRUE, check.names = FALSE)

  # Location for saved outputs -------------------------------------------------

  # check if parent directory exists
  parent_folder <- dirname(outfolder)
  if(!dir.exists(parent_folder)){
    message("Error: Please specify a valid path for the location 'outfolder' where files should be saved.")
    return()
  }

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # Get types of measure ----------------------------------------------------

  well_idx <- which(names(calibration_data) == "well")
  row_idx <- which(names(calibration_data) == "row")
  measures <- names(calibration_data)[(well_idx+1):(row_idx-1)]
  ### columns of measurements are looked for BETWEEN "well" column and "row" column

  # include: if defined, only include 'measures' with this substring
  if(!is.null(include_only)){

    # to include/exclude from tibble
    to_include <- measures[which(grepl(include_only, measures))]
    to_exclude <- setdiff(measures, to_include)
    calibration_data <- calibration_data[ , !names(calibration_data) %in% to_exclude]

    # redefine 'measures' list
    measures <- to_include

  }
  # exlude: if defined, only include 'measures' without this substring
  if(!is.null(exclude)){

    # to include/exclude from tibble
    to_exclude <- measures[which(grepl(exclude, measures))]
    to_include <- setdiff(measures, to_exclude)

    # head(calibration_data)
    calibration_data <- calibration_data[ , !names(calibration_data) %in% to_exclude]

    # redefine 'measures' list
    measures <- to_include

  }

  # Fix order if requested ---------------------------------------------------

  if(gain_fix == TRUE){

    for (gain in measures){

      if (!grepl("\\d{3}", gain)){ # for gains with fewer than 2 digits (in a row) in their name
        matches <- regmatches(gain, regexpr("([[:digit:]]){2}", gain)) # find location (regexpr) and identity (regmatches) ofthe 2 digits
        replacement <- paste0("0", matches) # add 0 before it
        newgain <- sub(pattern="([[:digit:]]){2}", replacement=replacement, x=gain)

        # replace:

        # names(calibration_data)
        colnum <- which(names(calibration_data) == gain)
        names(calibration_data)[colnum] <- newgain
        # names(calibration_data)

        # measures
        measure_num <- which(measures == gain)
        measures[measure_num] <- newgain
        # measures
      }
    }

  }

  # Rename channels if requested --------------------------------------------------

  if(!is.null(rename_from) & !is.null(rename_to)){

    for (meas in measures){ # for each measure-column..

      if (grepl(rename_from, meas)){ # if there's a match
        matches <- regmatches(meas, regexpr(rename_from, meas)) # find location (regexpr) and identity (regmatches) of rename_from
        newmeas <- sub(pattern=rename_from, replacement=rename_to, x=meas)

        # replace:

        # names(calibration_data)
        colnum <- which(names(calibration_data) == meas)
        names(calibration_data)[colnum] <- newmeas
        # names(calibration_data)

        # measures
        measure_num <- which(measures == meas)
        measures[measure_num] <- newmeas
        # measures
      }
    }

  }

  # Remove rows if requested -------------------------------------------------

  # always remove rows without protein annotated (empty wells)
  calibration_data <- calibration_data %>%
    dplyr::filter(.data$protein != "")

  # remove rows if requested
  if(subset_rows == TRUE){
    calibration_data <- calibration_data %>%
      dplyr::filter(.data$row %in% rows_to_keep)
    # dplyr::filter(.data$row %in% c("C","D"))
  }

  # Complete blank/buffer rows if requested ---------------------------------------
  # required for data that has been carried over from protein quant assays:
  # as input df only contains concentration data on the FP itself, not the blanks

  if(complete_blank){

    # if protein is "none", edit that row, otherwise leave data as is
    calibration_data <- calibration_data %>%
      # dplyr::mutate(concentration_ngul = ifelse(is.na(concentration_ngul), 0, concentration_ngul))
      dplyr::mutate(concentration_ngul = ifelse(.data$protein == "none", 0, .data$concentration_ngul))
    calibration_data

  }

  # Molecules calculation ----------------------------------------------

  # work out molecules from concentration_ngul and MW and volume
  calibration_data <- calibration_data %>%
    dplyr::mutate(mass_ng = .data$concentration_ngul*.data$volume) %>% # ng = ng/ul*ul
    dplyr::mutate(moles = .data$mass_ng*10^-9/.data$mw_gmol1) %>% # mol = g/(g/mol)
    dplyr::mutate(molecules = .data$moles*6*10^23) # molecules = mol * molecules/mol
  calibration_data

  # Data transformation steps -------------------------------------------------------

  # 1. Reshape data
  # 2. Normalise data
  # 3. Remove saturated values
  # 4. Summarise (take means)
  # 5. Find cfs fits
  # 6. Sensitivity calcs

  # 1. Reshape data -------------------------------------------------------
  # into long format

  raw_values <- calibration_data %>%
    tidyr::pivot_longer(tidyselect::all_of(measures),
                        names_to = "measure",
                        values_to = "raw_value")
  raw_values

  # 2. Normalise data -------------------------------------------------------

  norm_values <- c()

  for(calib in unique(raw_values$calibrant)){

    for(meas in measures){

      # create temp df. required so we can add in the extra columns.
      temp_meas_calib_values <- raw_values %>%
        dplyr::filter(.data$calibrant == calib) %>% # filter by calibrant
        dplyr::filter(.data$measure == meas) # filter rows by measure/gain
      temp_meas_calib_values

      # find mean of blanks
      blanks <- temp_meas_calib_values %>%
        dplyr::filter(.data$protein == "none") # find blanks
      blanks
      mean_raw_blanks <- mean(blanks$raw_value, na.rm = TRUE) # calc mean
      mean_raw_blanks

      # add new columns
      temp_meas_calib_values <- temp_meas_calib_values %>%
        dplyr::mutate(raw_blanks = mean_raw_blanks) %>%
        dplyr::mutate(normalised_value = .data$raw_value - .data$raw_blanks)
      temp_meas_calib_values

      # rbind normalised values
      norm_values <- rbind(norm_values, temp_meas_calib_values)
      norm_values

    } # measures
  } # calibrants

  # 3. Remove saturated values  -------------------------------------------------------
  # saturated values refers to those not in the linear range
  # saturation check POST normalisation allows us to retain more data points

  all_values <- c()
  non_sat_values <- c()

  for(calib in unique(norm_values$calibrant)){

    temp_calib_values <- norm_values %>%
      dplyr::filter(.data$calibrant == calib)
    temp_calib_values

    concentrations <- sort(unique(temp_calib_values$molecules))
    concentrations

    fold_dilution <- concentrations[3] / concentrations[2]
    # function is flexible to fold dilution bc works this out from the metadata
    fold_dilution

    high_saturation_threshold <- fold_dilution * 0.75 # saturation defined as 0.75*dilution factor
    high_saturation_threshold # 1.5

    # adding columns to record - dilution factor, max conc, dilution id
    temp_calib_values$dilution_ratio <- 1 / fold_dilution
    temp_calib_values$max_concentration <- max(concentrations)
    temp_calib_values$dilution_idx <- - log(temp_calib_values$max_concentration / temp_calib_values$molecules) / log(temp_calib_values$dilution_ratio)
    # complex way to get to numbering dilutions as #0-#10.

    temp_calib_values

    for(meas in measures){

      temp_meas_calib_values <- temp_calib_values %>%
        dplyr::filter(.data$measure == meas) # filter by gain
      temp_meas_calib_values

      # Need to refresh these values:
      # previous mean_raw_blanks value was for raw blank mean value of gain120 from prev section, ie 1000s.
      # here, mean_norm_blanks should all be 0, as we've normalised, and sd_norm_blanks will be variable.

      # find mean of blanks
      blanks <- temp_meas_calib_values %>%
        dplyr::filter(.data$protein == "none")
      blanks
      mean_norm_blanks <- mean(blanks$normalised_value, na.rm = TRUE) # calc mean
      mean_norm_blanks # 0

      # find sd of blanks
      sd_norm_blanks <- stats::sd(blanks$normalised_value, na.rm = TRUE)
      sd_norm_blanks

      # add new columns
      temp_meas_calib_values <- temp_meas_calib_values %>%
        dplyr::mutate(norm_blanks = mean_norm_blanks) %>%
        dplyr::mutate(norm_blanks_sd = sd_norm_blanks)
      temp_meas_calib_values

      # Join all_values table before some values get changed to NA
      all_values <- rbind(all_values, temp_meas_calib_values)
      # NB1. all_values joins together a just-created temp_meas_calib_values with new columns but unchanged normalised_values column
      # NB2. non_sat_values joins together a temp_meas_calib_values with new columns and saturation-screened (NA inserted) normalised_values column

      for(rplct in unique(temp_meas_calib_values$replicate)){
        # ...for each replicate ie for each dilution series (1,2,3,4) of difft concs...

        # Prevent error in next section if replicate section has no rows
        if(nrow(dplyr::filter(temp_meas_calib_values, replicate == rplct)) == 0){ next }

        prev_value <- 0
        this_value <- NA # add starting value of this_value so that it doesn't trip up on is.na(this_value) when starting at conc=0
        # previously we relied on what happened in top section to provide some value

        for(i in 1:length(concentrations)){

          current_conc_precise <- concentrations[i]

          # Prevent error in next section
          if(nrow(dplyr::filter(temp_meas_calib_values, temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct)) == 0){ next }

          if(current_conc_precise != 0){
            # Added for experiments where multiple conc 0s exist.
            # (Can't add this around the following "next" statement as "next" only refers to current if/for loop.)

            this_value <- temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"]
            this_value
          }

          if(is.na(this_value)){ next }

          # change saturated values to NA
          if(current_conc_precise != 0){

            ## check high saturation (or for values less than zero!)
            # makes sense to focus on this if going from conc 0 upwards.
            # prev value starts at 0, so first q checks if first non-zero conc gives positive result at all.
            # next iterations check if value is at least 1.5* the previous value, ie gives a MINIMUM.
            # this makes sense if starting from low concs as what you expect to encounter is a "lag phase to log phase" sigmoid curve (at least in the raw values)
            # ..in which values tend to not increase much at first, so the first priority is to make sure that they start to roughly double from the prev value.

            # ideally this_value is MORE THAN (1.5x previous) ... # this if statement is the PROBLEM condition!
            if(this_value <= prev_value * high_saturation_threshold){
              temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"] <- NA
            }

            ## check low saturation
            if(this_value <= mean_norm_blanks + 2 * sd_norm_blanks){
              temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"] <- NA
            }
          }
          prev_value <- this_value
          prev_value

        } # for i in 1:length conc

        # what i currently have - a 18*34 df of current calib/gain with some NAs..
        temp_meas_calib_values
        # what i want to do - get rid of anomalies at the bottom...

        # 3B. Check sat values in decreasing order too -------------------------------------

        concentrations_dec <- sort(unique(temp_calib_values$molecules), decreasing = TRUE)
        concentrations_dec

        # Previously worked out params:
        # fold_dilution <- concentrations[3] / concentrations[2]
        fold_dilution # 2
        # temp_calib_values$dilution_ratio <- 1 / fold_dilution
        temp_calib_values$dilution_ratio[1] # 0.5
        # high_saturation_threshold <- fold_dilution * 0.75
        high_saturation_threshold # 1.5
        # temp_calib_values$max_concentration <- max(concentrations)
        temp_calib_values$max_concentration[1] # 1.25e+15
        # temp_calib_values$dilution_idx <- - log(temp_calib_values$max_concentration / temp_calib_values$concentration) / log(temp_calib_values$dilution_ratio)

        # Find first non-NA values to be first comparator (prev_value) value:
        for(i in 1:length(concentrations_dec)){

          # record conc
          current_conc_precise <- concentrations_dec[i]

          # Prevent error in next section
          if(nrow(dplyr::filter(temp_meas_calib_values, temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct)) == 0){ next }

          # check if normalised value has a value
          if(current_conc_precise != 0){
            # Added for experiments where multiple conc 0s exist.
            # (Can't add this around the following "next" statement as "next" only refers to current if/for loop.)

            this_value <- NA # make sure
            this_value <- temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"]
            this_value
          }

          # if it has a value, break for loop - you have your starting conc and starting value
          if(!is.na(this_value)){ break }

        } # for i in 1:length

        prev_value <- this_value
        prev_value

        # Narrow the range of concs to look at
        first_conc <- which(concentrations_dec == current_conc_precise) + 1
        first_conc
        last_conc <- which(concentrations_dec == dplyr::last(concentrations_dec))
        last_conc
        concentrations_dec_nonoverflow <- concentrations_dec[first_conc:last_conc]
        concentrations_dec_nonoverflow

        # Remove saturated values2 - top to bottom
        for(i in 1:length(concentrations_dec_nonoverflow)){
          # starts from highest conc

          current_conc_precise <- concentrations_dec_nonoverflow[i]

          # Prevent error in next section
          if(nrow(dplyr::filter(temp_meas_calib_values, temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct)) == 0){ next }

          if(current_conc_precise != 0){
            # Added for experiments where multiple conc 0s exist.
            # (Can't add this around the following "next" statement as "next" only refers to current if/for loop.)

            this_value <- temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"]
            this_value
          }

          # if value is NA
          # essentially - let's pretend this_value was the highest tolerable value, which we pass on to prev_value before proceeding
          if(is.na(this_value)){
            prev_value <- prev_value / high_saturation_threshold
            next
          }

          # change saturated values to NA
          if(current_conc_precise != 0){

            ## check high saturation (or for values less than zero!)
            # for first section, iterations check if value is at least 1.5* the previous value, ie gives a MINIMUM.
            # if(this_value <= prev_value * high_saturation_threshold){ ...there is a problem... }
            # this makes sense if starting from low concs as what you expect to encounter is a "lag phase to log phase" sigmoid curve (at least in the raw values)
            # ..in which values tend to not increase much at first, so the first priority is to make sure that they start to roughly double from the prev value.

            # From top down, you have the opposite issue. You want to test if the lower value is at least 1.5x less than higher value.
            # You are setting a MAXIMUM.

            # ideally this_value is LESS THAN (prev/1.5) ... # this if statement is the PROBLEM condition!
            if(this_value >= prev_value / high_saturation_threshold){
              temp_meas_calib_values[temp_meas_calib_values$molecules == current_conc_precise & temp_meas_calib_values$replicate == rplct,"normalised_value"] <- NA
            }
          } # if conc 0

          prev_value <- this_value
          prev_value

          # all conc loop does is iteratively EDIT a pre-existing df (temp_meas_calib_values) so nothing needs joining here

        } # for conc in concentrations_dec_nonoverflow

      } # for rplct in replicate

      # New temp_meas_calib_values created in each loop of measures, so needs joining here.
      non_sat_values <- rbind(non_sat_values, temp_meas_calib_values)
      # temp blocks are per calibrant, per gain. so you'll get calibrant1: gain40-120; calibrant2: gain40-120..
      # NB1. all_values joins together a just-created temp_meas_calib_values with new columns but unchanged normalised_values column
      # NB2. non_sat_values joins together a temp_meas_calib_values with new columns and saturation-screened (NA inserted) normalised_values column

      # as for loop runs (potentially) multiple times, once for each calibrant..
      # non_sat_values simply provides a way to rowbind the calibrants' datasets
      # back together (this time without saturated values)

    } # for meas in measure

  } # for each calibrant

  non_sat_values

  # 4. Summarise  -------------------------------------------------------

  # calculate mean of n replicates, for raw_values and normalised_values columns, for non_sat_values and all_values dfs
  # a. summ_values_nonsat uses non_sat_values and therefore does not include saturated values.
  summ_values_nonsat <- non_sat_values %>%

    # For each dilution and wavelength
    dplyr::group_by(.data$concentration_ngul, .data$molecules, .data$dilution_idx,
                    .data$measure,
                    .drop = FALSE) %>% # don't remove groups w no values

    # Take mean of raw and normalised values
    dplyr::mutate(raw_value = mean(.data$raw_value, na.rm = TRUE)) %>%
    dplyr::mutate(normalised_value = mean(.data$normalised_value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%

    # Tidy
    dplyr::select(-c(.data$replicate, .data$well, .data$row, .data$column
                     # .data$kfactor_1cm, .data$kfactor_well,
                     # .data$pathlength_each, .data$pathlength_blanks, .data$pathlength_volume, .data$pathlength,
                     # .data$raw_cm1_blanks
    )) %>%
    dplyr::distinct() # remove duplicate rows
  # dplyr::arrange(dplyr::desc(.data$dilution_idx)) # arrange by dilution, starting from 1
  summ_values_nonsat

  # b. summ_values_all uses all_values and therefore includes both sat/non_sat points
  summ_values_all <- all_values %>%

    # For each dilution and wavelength
    dplyr::group_by(.data$concentration_ngul, .data$molecules, .data$dilution_idx,
                    .data$measure,
                    .drop = FALSE) %>% # don't remove groups w no values

    # Take mean of raw and normalised values
    dplyr::mutate(raw_value = mean(.data$raw_value, na.rm = TRUE)) %>%
    dplyr::mutate(normalised_value = mean(.data$normalised_value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%

    # Tidy
    dplyr::select(-c(.data$replicate, .data$well, .data$row, .data$column
                     # .data$kfactor_1cm, .data$kfactor_well,
                     # .data$pathlength_each, .data$pathlength_blanks, .data$pathlength_volume, .data$pathlength,
                     # .data$raw_cm1_blanks
    )) %>%
    dplyr::distinct() # remove duplicate rows
  # dplyr::arrange(dplyr::desc(.data$dilution_idx)) # arrange by dilution, starting from 1
  summ_values_all

  # 5. Find CFS fits  -------------------------------------------------------
  # conversion factors (cfs) defined here as cf = rfu/molecules
  # can be found (i) by taking mean of conversion factors at each dilution, or (ii) by adding in pipetting error (beta)
  # and fitting a 'pipetting error model' to minimise error between the cf and rfu/molecules ratio
  # where actual molecule number in each dilution is adjusted (b_i) to take account of this pipetting error (beta)
  # as in 2020 Beal et al Comms Bio https://doi.org/10.1038/s42003-020-01127-5

  ### Think of these as temporary dfs, just for the fit. The final dfs will be summ_values_nonsat/fit_values. ###
  # Remove 0s for fit value fn below, but keep 0s in separate table for rbinding at the end...
  trimmed_values_blanks <- summ_values_nonsat %>% dplyr::filter(.data$protein == "none") # save first into new df
  trimmed_values <- summ_values_nonsat %>% dplyr::filter(.data$protein != "none") # overwrite long_values df
  trimmed_values
  trimmed_values <- stats::na.omit(trimmed_values) # ESSENTIAL to allow stats::optim to work. NAs and 0s need to be removed.
  trimmed_values

  ### Model Fit

  ### create table of fit values (conversion factors) by fitting data to some model
  fit_values <- c()

  for(calib in unique(trimmed_values$calibrant)){ # for each calibrant

    temp_calib_values <- trimmed_values %>% dplyr::filter(.data$calibrant == calib) # temporarily subset df for current calibrant

    for(meas in unique(temp_calib_values$measure)){ # for each measure/gain

      temp_meas_calib_values <- temp_calib_values %>%
        dplyr::filter(.data$measure == meas) # temp filter by gain

      # if number of points in any one measure/gain is less than three, don't bother trying to fit a model.
      if(nrow(temp_meas_calib_values) < 3){
        next
      }

      model <- 0

      error_func <- function(x){
        data <- temp_meas_calib_values

        # this function is here for the optim() function below.
        # optim stipulates that the function provided to optim must take the parameters of optim(par = ...) as its first argument
        # so x = c(1e-10,0)

        cf <- x[1] # initial value of cf is always 1e-10
        beta <- x[2] # initial value of beta is always 0
        error <- 0

        for(i in data$dilution_idx){ # for each dilution
          data_i <- data[data$dilution_idx == i,] # subset data

          # 2020 Beal et al Comms Bio: The model for systematic pipetting error modifies the intended dilution factor dilution_ratio
          # with the addition of an unknown bias beta, such that the expected biased population b_i for the ith dilution level is:
          b_i <- data_i$max_concentration * (1 - data_i$dilution_ratio - beta) *
            (data_i$dilution_ratio + beta) ^ (data_i$dilution_idx - 1) # see 2020 Beal et al Comms Bio

          # The squared error is:
          e_i <- abs(log10(cf * b_i / data_i$normalised_value))^2

          # The sum squared error over all valid dilution levels is:
          error <- error + e_i
        }

        return(error)
      }

      # 2020 Beal et al Comms Bio:
      # Simultaneously fit beta and the scaling factor cf to minimize the sum squared error over all valid dilution levels.

      ## Minimise the error using optim:
      # optim() finds the parameter values in its 1st arg that minimise the output of its 2nd, fn arg
      # here:
      # argument1: par = c(1e-10,0) # Initial values for the parameters to be optimized over.
      # argument2: fn = error_func
      # returns a list with following components:
      # res$par # The best set of parameters found.
      # res$value # the value of the function corresponding to par
      # res$convergence # an integer code. if 0, this indicates successful completion!
      res <- stats::optim(c(1e-10,0), error_func)
      # !! requires no 0 or NA values in normalised_value column.

      ## Residuals
      # 2020 Beal et al Comms Bio: The residuals for this fit are then the absolute ratio of fit-predicted to observed net mean .. for all valid levels
      # `(bi/Sp) / mean(Oi)-mean(B)` which here would be `(b_i/cf) / normalised_value`
      # but 2020 Fedorec et al used `cf * b_i/normalised_value` in error model, rather than Beal's `1/cf * b_i/normalised_value`
      # so residuals here should be `cf*b_i/normalised_value`
      # as the 'prediction' for the rfu are given by cf*molecules, where cf is defined as cf=rfu/molecules

      # Get residuals at all dilutions
      residuals_table <- c()
      for(i in temp_meas_calib_values$dilution_idx){

        # nb
        # res$par[1] is cf
        # res$par[2] is beta

        data_i <- temp_meas_calib_values[temp_meas_calib_values$dilution_idx == i,] # subset data
        b_i <- data_i$max_concentration * (1 - data_i$dilution_ratio - res$par[2]) *
          (data_i$dilution_ratio + res$par[2]) ^ (data_i$dilution_idx - 1) # work out b_i again using res$par values for beta
        # e_i <- abs(log10(cf * b_i / data_i$normalised_value))^2
        # error <- error + e_i

        residual_i <- res$par[1] * b_i/data_i$normalised_value

        table_i <- data.frame(
          dilution_idx = data_i$dilution_idx,
          residual = abs(residual_i) # take absolute of residual
        )
        residuals_table <- rbind(residuals_table, table_i)
      }
      residuals_table

      # 2020 Beal et al Comms Bio: Criteria for valid fluorescein dilution defined as
      # Systematic pipetting error has geometric mean absolute residual <1.1-fold.
      # Geometric mean function:
      gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
        if(any(x < 0, na.rm = TRUE)){
          return(NaN)
        }
        if(zero.propagate){
          if(any(x == 0, na.rm = TRUE)){
            return(0)
          }
          exp(mean(log(x), na.rm = na.rm))
        } else {
          exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
      }
      # Geometric mean of the absolutes of each residual
      residual <- gm_mean(residuals_table$residual)
      residual

      if(res$convergence == 0){
        ### save parameters of model - 1st parameter = conversion factor, 2nd = beta

        new_fit <- data.frame(instrument = temp_meas_calib_values$instrument[1],
                              plate = temp_meas_calib_values$plate[1],
                              seal = temp_meas_calib_values$seal[1],
                              channel_name = temp_meas_calib_values$channel_name[1],
                              channel_ex = temp_meas_calib_values$channel_ex[1],
                              channel_em = temp_meas_calib_values$channel_em[1],
                              media = temp_meas_calib_values$media[1],

                              calibrant = calib, # calib is the current element in calibrant for loop
                              measure = meas, # meas is the current element in measure (eg GFP040) for loop

                              cf = res$par[1], beta = res$par[2],
                              residuals = residual
        )
        fit_values <- rbind(fit_values, new_fit)
      } # res$convergence

    }  # for each measure/gain
  } # for each calibrant

  ### Prepare final dfs ----------------------------------------------

  if(!is.null(fit_values)){

    # add cf and beta columns to the right of long_values tibble:
    trimmed_values <- dplyr::full_join(trimmed_values, fit_values)
    # put the blanks back
    trimmed_values <- dplyr::bind_rows(trimmed_values_blanks, trimmed_values) # avoid rbind error due to trimmed_values_blanks missing fit_values' columns
    # Gain needed for standard plots as well as others. Add now:
    trimmed_values$gain <- as.numeric(gsub(paste0(trimmed_values$channel_name[1], separator), "", trimmed_values$measure))
    # trimmed_values$gain # 40  50  60  70  80  90 100
    # Move gain to right of measure
    trimmed_values <- trimmed_values %>%
      dplyr::relocate(gain, .after = .data$measure)

    # Add a check for channel_name and channel in measure being the same thing!!
    if(grepl(fit_values$channel_name[1], fit_values$measure[1])){

      # Gain needed for standard plots as well as others. Add now:
      fit_values$gain <- as.numeric(gsub(paste0(fit_values$channel_name[1], separator), "", fit_values$measure))
      # Move gain to right of measure
      fit_values <- fit_values %>%
        dplyr::relocate(gain, .after = .data$measure)
      # Sort by gain (as sometimes gains get mixed up..)
      fit_values <- fit_values[order(fit_values$gain),]
      fit_values

    } else {
      message("Entries provided for the 'channel name' column in the metadata csv are not contained within the 'measure' column in the data csv.
              This can happen when the 'channel name' column in the metadata is not adjusted between different scans,
              or if the naming schemes for the channel names provided in the plate reader methods change, from eg 'GFP' to 'green' or 'greenamber'.")
      fit_values$gain <- NA
      fit_values <- fit_values %>%
        dplyr::relocate(gain, .after = .data$measure)

    }

  } else {
    message("Warning: No conversion factors could be fitted.")
  }

  # end of 5. CFS fit ----

  # 6. Sensitivity calcs -----------------------------------------------

  # Add variables: min/max_normflu, min/max_mols -------------------------------------------------

  if(sensitivity_plots){

    trimmed_values <- trimmed_values %>%
      dplyr::mutate(min_normflu = 1) %>% # normalised min fluorescence value that can be distinguished from buffer (0) is 1.
      dplyr::mutate(min_mols = .data$min_normflu/.data$cf)
    # MAX
    trimmed_values <- trimmed_values %>%
      dplyr::group_by(.data$measure) %>%
      dplyr::mutate(max_normflu = 100000 - # raw max = 100,000 bc depends on instrument
                      .data$raw_blanks) %>% # normalised max fluorescence = raw max - blank (for that measure)
      dplyr::mutate(max_mols = .data$max_normflu/.data$cf)
    trimmed_values
  }

  # Amend fit_values to include these values -------------------------------------------------

  if(sensitivity_plots){

    fit_values2 <- c()

    for(meas in unique(trimmed_values$measure)){ # for each measure/gain

      # Data from trimmed_values
      temp_trimmed_values <- trimmed_values %>%
        dplyr::filter(.data$measure == meas) # temp filter by gain

      new_fit <- data.frame(instrument = temp_trimmed_values$instrument[1],
                            plate = temp_trimmed_values$plate[1],
                            seal = temp_trimmed_values$seal[1],
                            channel_name = temp_trimmed_values$channel_name[1],
                            channel_ex = temp_trimmed_values$channel_ex[1],
                            channel_em = temp_trimmed_values$channel_em[1],

                            media = temp_trimmed_values$media[1],
                            calibrant = temp_trimmed_values$calibrant[2], # [2] might help if a BSA/none gets in there
                            measure = meas, # meas is the current element in measure (eg GFP040) for loop
                            gain = temp_trimmed_values$gain[1],

                            # new columns:
                            min_normflu = temp_trimmed_values$min_normflu[2],
                            min_mols = dplyr::first(stats::na.omit(temp_trimmed_values$min_mols)), # first non NA
                            max_normflu = temp_trimmed_values$max_normflu[2],
                            max_mols = dplyr::first(stats::na.omit(temp_trimmed_values$max_mols)), # first non NA

                            # This doesn't work if BSA rows in fluor assays.
                            cf = dplyr::first(stats::na.omit(temp_trimmed_values$cf)), # first value that isn't NA
                            beta = dplyr::first(stats::na.omit(temp_trimmed_values$beta)), # first value that isn't NA
                            residuals = dplyr::first(stats::na.omit(temp_trimmed_values$residuals)) # first value that isn't NA
      )
      fit_values2 <- rbind(fit_values2, new_fit)

    } # for each measure/gain

    # Overwrite standard fit_values to allow (a) saving (b) return - of the expanded fit values table
    fit_values <- fit_values2
    fit_values

  } # incl_flu, sens_plots

  # 7. Save CSVs  -------------------------------------------------------

  if(more_csvs){

    # # 1. Transformed raw data
    # # redundant bc of norm_values
    # utils::write.csv(x = raw_values, file = file.path(outfolder, "1_raw_values.csv"), row.names = FALSE)

    # # 2. Normalised data
    # # norm_values = raw_values with 2 extra columns: raw_blanks and normalised_values
    # # redundant bc of all_values
    # utils::write.csv(x = norm_values, file = file.path(outfolder, "2_normalised_values.csv"), row.names = FALSE)

    # 3. Add more columns and make two dfs:
    # 3a. All values - for extra 5 columns and no data removal
    # all_values = norm_values with 5 extra columns: dilution_ratio, max_conc, dilution_idx, norm_blanks, norm_blanks_sd
    # no data removal
    utils::write.csv(x = all_values, file = file.path(outfolder, "3a_all_values.csv"), row.names = FALSE)
    # 3b. Non saturated values - for extra 5 columns and saturated data removal (replacement with NAs)
    # non_sat_values = norm_values with 5 extra columns: dilution_ration, max_conc, dilution_idx, norm_blanks, norm_blanks_sd
    # but also with saturated values replaced by NAs
    utils::write.csv(x = non_sat_values, file = file.path(outfolder, "3b_non_sat_values.csv"), row.names = FALSE)

    # 4. Summary values (means)
    # 4a. summ_values_all - norm_values summarised to give means per replicate but without taking out the saturated values
    utils::write.csv(x = summ_values_all, file = file.path(outfolder, "4a_summ_values_all.csv"), row.names = FALSE)
    # 4b. summ_values_nonsat - with saturated values (NAs), from non_sat_values
    utils::write.csv(x = summ_values_nonsat, file = file.path(outfolder, "4b_summ_values_nonsat.csv"), row.names = FALSE)

    # 5. Summarised values (with NAs missing) plus fits
    # trimmed_values = summ_values_nonsat, minus all rows with NA (saturated) normalised values, PLUS CF/BETA
    # NOT REDUNDANT: can't plot fits with other df else fit line extends further than its relevant points!
    # also reqd for plotting raw values that are not saturated bc na.omit removes entire rows where saturated norm values are.
    utils::write.csv(x = trimmed_values, file = file.path(outfolder, "5_trimmed_values.csv"), row.names = FALSE)

    # 6. Fits
    # fit_values = just the CF/BETA fits for each gain (with min extra cols).
    utils::write.csv(x = fit_values, file = file.path(outfolder, "6_conversion_factors.csv"), row.names = FALSE)

  } else {

    # save conversion factors only
    filename <- gsub(".csv", "_cfs.csv", basename(calibration_csv))
    utils::write.csv(fit_values, file.path(outfolder, filename), row.names = FALSE)

  }

  # 8. Make plots  -------------------------------------------------------

  # Standard calibration plots ------

  # top down diagonal calibration plot
  flu_plt <- ggplot2::ggplot(data = trimmed_values %>%
                               dplyr::filter(.data$molecules != 0)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution_idx,
                                     y = .data$normalised_value)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$dilution_idx,
                                    y = .data$cf * .data$max_concentration *
                                      (1 - .data$dilution_ratio - beta) *
                                      (.data$dilution_ratio + beta) ^ (.data$dilution_idx - 1))) +
    ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10",
                                breaks = c(1e+00, 1e+01, 1e+02, 1e+03, 1e+04, 1e+05)) +
    ggplot2::scale_x_continuous("Dilution # (1:2 dilution series)",
                                breaks = seq(0,22,2)) +
    ggplot2::facet_wrap(~measure) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio = 1,
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) ###
  flu_plt
  plotname <- gsub(".csv", "_cfs_flu1.pdf", basename(calibration_csv))
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = flu_plt,
                  width = 18, height = 18, units = "cm")

  # molecule number on x, linear x/y
  flu_plt2 <- ggplot2::ggplot(data = trimmed_values %>%
                                dplyr::filter(.data$molecules != 0)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$molecules,
                                     y = .data$normalised_value)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$molecules,
                                    y = .data$cf * .data$max_concentration *
                                      (1 - .data$dilution_ratio - beta) *
                                      (.data$dilution_ratio + beta) ^ (.data$dilution_idx - 1))) +
    ggplot2::scale_y_continuous("Normalised Fluorescence") +
    ggplot2::scale_x_continuous("Molecules") +
    ggplot2::facet_wrap(~measure) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio = 1,
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0))
  flu_plt2
  plotname <- gsub(".csv", "_cfs_flu2.pdf", basename(calibration_csv))
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = flu_plt2,
                  width = 18, height = 18, units = "cm")

  # molecule number on x, log-scale x/y
  flu_plt3 <- ggplot2::ggplot(data = trimmed_values %>%
                                dplyr::filter(.data$molecules != 0)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$molecules,
                                     y = .data$normalised_value)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$molecules,
                                    y = .data$cf * .data$max_concentration *
                                      (1 - .data$dilution_ratio - beta) *
                                      (.data$dilution_ratio + beta) ^ (.data$dilution_idx - 1))) +
    ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10",
                                breaks = c(1e+00, 1e+01, 1e+02, 1e+03, 1e+04, 1e+05)) +
    ggplot2::scale_x_continuous("Molecules", trans = "log10",
                                breaks = c(1e+08, 1e+09, 1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15, 1e+16)) + ###
    ggplot2::facet_wrap(~measure) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio = 1,
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0))
  flu_plt3
  plotname <- gsub(".csv", "_cfs_flu3.pdf", basename(calibration_csv))
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = flu_plt3,
                  width = 18, height = 18, units = "cm")

  # Gain vs cfs plot ----------------------------------------

  # Fit cf to Gain relation model
  model <- stats::lm(log10(cf) ~ poly(gain, 2), data = fit_values)
  # This can fail if fitvalues5 overwrites fit_values with NAs instead of cf values.

  # Plot points and model
  plt_gain <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = fit_values$gain,
                                    y = 10^stats::predict(model, fit_values))) +
    ggplot2::geom_point(ggplot2::aes(x = fit_values$gain,
                                     y = fit_values$cf)) +
    ggplot2::scale_x_continuous("Gain", limits = c(40,120), breaks = c(40,60,80,100,120)) +
    ggplot2::scale_y_continuous("Conversion factor\n(fluorescence/molecule)",
                                # cf = rfu/molecule. molec * cf = relative flu units (eg 10^15 molec * 10^-10 cf = 10^5 rfu)
                                # can call this ... rfu/molecule ... or fluorescence/molecule...
                                trans = "log10",
                                limits = c(1e-13, 1e-7),
                                breaks = c(1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07)) +
    ggplot2::theme_bw() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(aspect.ratio = 1,
                   panel.grid.minor = ggplot2::element_blank())
  plt_gain

  plotname <- gsub(".csv", "_cfs_flu4.pdf", basename(calibration_csv))
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plt_gain,
                  height = 10, width = 10, units = "cm")

  # Plots of all values -----------------------------------------

  if(more_plots){

    # Plot1. all points, raw.
    # log scale (can't see it well on linear scale)
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summ_values_nonsat %>%
                            # raw_values and norm_values equivalent for raw value columns.
                            # summ_values_nonsat loses a few points bc of taking of the mean (tidier)
                            # for trimmed_values unfortunately i had to remove all rows where normalised_value is NA so some raw values disappeared too
                            dplyr::filter(.data$molecules != 0),
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$raw_value),
                          colour = "red") +
      ggplot2::geom_hline(data = summ_values_nonsat %>%
                            dplyr::filter(.data$molecules == 0),
                          ggplot2::aes(yintercept = .data$raw_value),
                          colour = "red", alpha = 0.5) +
      ggplot2::scale_y_continuous("Raw Fluorescence", trans = "log10", breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::scale_x_continuous("Molecules", trans = "log10",
                                  breaks = c(1,1e+4,1e+8,1e+12,1e+16)) + ###
      ggplot2::labs(title = "Raw data with LOD") +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, ###
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
      ggplot2::coord_cartesian(xlim = c(1, 1e+16), ylim = c(1, 100000))
    flu_plt3
    plotname <- gsub(".csv", "_cfs_flu7a.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 18, height = 18, units = "cm")

    # Plot2. all points, normalised.
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summ_values_nonsat %>%
                            # raw_values and norm_values equivalent for raw value columns.
                            # summ_values_nonsat loses a few points bc of taking of the mean (tidier)
                            # for trimmed_values unfortunately i had to remove all rows where normalised_value is NA so some raw values disappeared too
                            dplyr::filter(.data$molecules != 0),
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$normalised_value),
                          colour = "red") +
      ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10", breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::scale_x_continuous("Molecules", trans = "log10",
                                  breaks = c(1,1e+4,1e+8,1e+12,1e+16)) + ###
      ggplot2::labs(title = "Normalised data") +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, ###
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
      ggplot2::coord_cartesian(xlim = c(1, 1e+16), ylim = c(1, 100000))
    flu_plt3
    plotname <- gsub(".csv", "_cfs_flu7b.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 18, height = 18, units = "cm")

    # Comparisons of non-sat/sat values --------------------------------------------

    # Plot 1 - all values vs non-sat values with this version
    # log scale (can't see it well on linear scale)
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summ_values_nonsat %>%
                            # raw_values and norm_values equivalent for raw value columns.
                            # summ_values_nonsat loses a few points bc of taking of the mean (tidier)
                            # for trimmed_values unfortunately i had to remove all rows where normalised_value is NA so some raw values disappeared too
                            dplyr::filter(.data$molecules != 0),
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$raw_value),
                          colour = "red") +
      ggplot2::geom_hline(data = summ_values_nonsat %>%
                            dplyr::filter(.data$molecules == 0),
                          ggplot2::aes(yintercept = .data$raw_value),
                          colour = "red", alpha = 0.5) +
      ggplot2::geom_point(data = trimmed_values %>%
                            # summ_values_nonsat best for all raw values as some raw values removed by the na.omit
                            # but trimmed_values best for non_sat raw_values, as raw values weren't checked for saturation at all!
                            # the only way to work out non-sat raw values is to remove rows that had NA in normalised_values
                            # ie what i did for trimmed_values.
                            dplyr::filter(.data$molecules != 0),
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$raw_value),
                          colour = "green") + # overlay the green!
      # removed fit line here as that's technically to normalised values not raw values
      ggplot2::scale_y_continuous("Raw Fluorescence", trans = "log10", breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::scale_x_continuous("Molecules", trans = "log10",
                                  breaks = c(1,1e+4,1e+8,1e+12,1e+16)) + ###
      ggplot2::labs(
        title = "All points vs non-saturated points",
        caption = "all = red; non-sat = green"
      ) +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, ###
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
      ggplot2::coord_cartesian(xlim = c(1, 1e+16), ylim = c(1, 100000))
    flu_plt3

    plotname <- gsub(".csv", "_cfs_flu9a.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 18, height = 18, units = "cm")

    # Plot2 - all points (red) vs non-sat values (green points) and green cfs line
    flu_plt3 <- ggplot2::ggplot() +
      # all points
      ggplot2::geom_point(data = summ_values_all, # all
                          # summ_values_all required for this. summarised version of norm_values that was not put through sat check
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$normalised_value),
                          colour = "red") +
      # non-saturated points
      ggplot2::geom_point(data = summ_values_nonsat, # sat points removed post normalisation
                          # can use either summ_values_nonsat or trimmed_values here
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$normalised_value),
                          colour = "green") +
      ggplot2::geom_line(data = trimmed_values,
                         ggplot2::aes(x = .data$molecules,
                                      y = .data$cf * .data$max_concentration *
                                        (1 - .data$dilution_ratio - beta) *
                                        (.data$dilution_ratio + beta) ^ (.data$dilution_idx - 1)),
                         colour = "green") +
      ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10", breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::scale_x_continuous("Molecules", trans = "log10", breaks = c(1,1e+4,1e+8,1e+12,1e+16)) + ###
      ggplot2::labs(
        title = "All points vs non-saturated points",
        caption = "all=red; non-sat = green"
      ) +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, ###
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
      ggplot2::coord_cartesian(xlim = c(1, 1e+16), ylim = c(1, 100000))
    flu_plt3

    plotname <- gsub(".csv", "_cfs_flu9b.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 18, height = 18, units = "cm")

  } # more plots

  # Sensitivity plots -------------------------------------------------------

  if (sensitivity_plots) {

    # raw blank values vs gain ----------------------------------------------

    # Plot raw blank values against Gain:
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = trimmed_values %>%
                            dplyr::filter(.data$molecules == 0),
                          ggplot2::aes(x = gain,
                                       y = .data$raw_value)) +
      ggplot2::geom_line(data = trimmed_values %>%
                           dplyr::filter(.data$molecules == 0),
                         ggplot2::aes(x = gain,
                                      y = .data$raw_value)) +
      ggplot2::scale_x_continuous("Gain", limits = c(40,120)) +
      ggplot2::scale_y_continuous("Raw fluorescence", trans = "log10", limits = c(1,100000),
                                  breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::ggtitle("Raw buffer values vs Gain") +
      ggplot2::labs(colour = "") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1,
                     panel.grid.minor = ggplot2::element_blank())
    flu_plt3
    plotname <- gsub(".csv", "cfs_flu11.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 10, height = 10, units = "cm")

    # sensitivity plots: detection range -------------------------------------------------------

    # Plot min/max values on normalised log plot
    flu_plt3 <- ggplot2::ggplot(data = trimmed_values) +
      ggplot2::geom_point(data = trimmed_values %>%
                            dplyr::filter(.data$normalised_value >=1),
                          ggplot2::aes(x = .data$molecules,
                                       y = .data$normalised_value)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$min_normflu),
                          colour = "red", alpha = 0.5) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$max_normflu),
                          colour = "red", alpha = 0.5) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$min_mols),
                          colour = "black", alpha = 0.5) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$max_mols),
                          colour = "black", alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(x = .data$molecules,
                                      y = .data$cf * .data$max_concentration *
                                        (1 - .data$dilution_ratio - beta) *
                                        (.data$dilution_ratio + beta) ^ (.data$dilution_idx - 1)),
                         colour = "black") +
      ggplot2::scale_y_continuous("Normalised Fluorescence", trans = "log10", breaks = c(1,10,100,1000,10000,100000)) +
      ggplot2::scale_x_continuous("Molecules", trans = "log10",
                                  breaks = c(1,1e+4,1e+8,1e+12,1e+16)) + ###
      ggplot2::labs(
        title = "Normalised data with min/max detection limits",
        caption = "black lines=detection limits"
      ) +
      ggplot2::facet_wrap(~measure) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, ###
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
      ggplot2::coord_cartesian(xlim = c(1, 1e+16), ylim = c(0.5, 100000))
    flu_plt3
    plotname <- gsub(".csv", "cfs_flu12.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 18, height = 18, units = "cm")

    # Plot min and max with gain
    flu_plt3 <- ggplot2::ggplot(data = fit_values2) +
      ggplot2::geom_line(ggplot2::aes(x = gain,
                                      y = .data$min_mols,
                                      colour = "min mols")) +
      ggplot2::geom_line(ggplot2::aes(x = gain,
                                      y = .data$max_mols,
                                      colour = "max mols")) +
      ggplot2::geom_ribbon(ggplot2::aes(x = gain, ymin = .data$min_mols, ymax = .data$max_mols, fill = "detection range"),
                           # fill = "grey50",
                           alpha = 0.1) +
      ggplot2::scale_x_continuous("Gain", limits = c(40,120)) +
      ggplot2::scale_y_continuous("Molecules", trans = "log10",
                                  limits = c(1,1e+18), breaks = c(1,1e+4,1e+8,1e+12,1e+16)) +
      ggplot2::labs(title = "Limits of detection") +

      # Legend just for fill
      ggplot2::guides(colour = FALSE) + # remove legend for colours
      ggplot2::scale_color_manual(breaks = c("min mols", "max mols"),
                                  values = c("grey80", "grey80")) + # specify colours for lines
      ggplot2::labs(fill = "") + # include legend for fill but remove its title
      ggplot2::scale_fill_manual(breaks = c("detection range"),
                                 values = c("grey50")) + # specify colours for fill

      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1,
                     legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank())
    flu_plt3
    plotname <- gsub(".csv", "cfs_flu13.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 12, height = 12, units = "cm")

    # Plot sensitivity of each gain (using min_mols)
    # Work out min_mols of gain40
    min_mols_g40 <- as.numeric(fit_values2 %>%
                                 dplyr::filter(gain == 40) %>%
                                 dplyr::select(.data$min_mols))
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = fit_values2,
                          ggplot2::aes(x = gain,
                                       y = min_mols_g40/.data$min_mols),
                          colour = "black") +
      ggplot2::geom_line(data = fit_values2,
                         ggplot2::aes(x = gain,
                                      y = min_mols_g40/.data$min_mols),
                         colour = "black") +
      ggplot2::scale_x_continuous("Gain", limits = c(40,120)) +
      ggplot2::scale_y_continuous("molecules (Gain40) / molecules (Gain)", trans = "log10") +
      ggplot2::ggtitle("Relative sensitivity of each gain") +
      ggplot2::labs(colour = "") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1,
                     panel.grid.minor = ggplot2::element_blank())
    flu_plt3
    plotname <- gsub(".csv", "cfs_flu14.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 10, height = 10, units = "cm")

    # Plot dynamic range of each gain
    flu_plt3 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = fit_values2,
                          ggplot2::aes(x = gain,
                                       y = .data$max_mols/.data$min_mols),
                          colour = "black") +
      ggplot2::geom_line(data = fit_values2,
                         ggplot2::aes(x = gain,
                                      y = .data$max_mols/.data$min_mols),
                         colour = "black") +
      ggplot2::scale_x_continuous("Gain", limits = c(40,120)) +
      ggplot2::scale_y_continuous("max/min detectable molecules", limits = c(8.9e+4,1.01e+5),
                                  breaks = c(9.0e+4, 9.2e+4, 9.4e+4, 9.6e+4, 9.8e+4, 10e+4)) +
      ggplot2::ggtitle("Dynamic range by gain") +
      ggplot2::labs(colour = "") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1,
                     panel.grid.minor = ggplot2::element_blank())
    flu_plt3
    plotname <- gsub(".csv", "cfs_flu16.pdf", basename(calibration_csv))
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = flu_plt3,
                    width = 10, height = 10, units = "cm")

  } # sensitivity plots

  # Return ----------------------------------------

  return(fit_values)

}
