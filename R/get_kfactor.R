#' Get k-factor of a buffer
#'
#' Calculates k-factor of a given buffer/temperature combination according to
#' data from Thermo Scientific for different buffers and temperatures. The
#' k-factor is defined as the A975-A900 for a 1cm pathlength. The calculations
#' are made under the assumption that the fold changes for buffers from water,
#' and for temperatures from 25oC may be combined to get an approximate estimate
#' of the k-factor for a given condition.
#'
#' @param buffer_used string corresponding to buffer. Must match a buffer entry
#'   in dataset `kfactors_buffers_data`. Default is "water".
#' @param concentration_used numeric value of concentration of buffer used.
#'   Default is 0 (pure water), so it needs changing if buffer isn't water. Take
#'   care to ensure the units are correct (buffers specified in M require
#'   concentrations in M not mM). This value is ignored if buffer = "water".
#' @param temperature_used numeric value of temperature in oC. Default is 25.
#'
#' @export
#'
#' @examples get_kfactor(buffer_used = "TBS", concentration_used = 0.005, temperature_used = 30)
get_kfactor <- function(buffer_used = "water", concentration_used = 0, temperature_used = 25
){

  # Summary --------------------------------------------------

  message(paste0("\nCalculating k-factor for ", buffer_used, " at concentration ", concentration_used, " at temperature ", temperature_used, " oC."))

  # Get k-factor data -------------------------------------------------

  # Get buffers data
  kfactors_buffers_data <- fpcountr::kfactors_buffers_data

  # Get temperature data
  kfactors_temperature_data <- fpcountr::kfactors_temperature_data

  # Get k-factor - reference -------------------------------------------------

  # Set reference value: water at 25oC according to kfactors_temperature_data = 0.172
  reference_kfactor <- kfactors_temperature_data %>%
    dplyr::filter(.data$temperature == 25) %>%
    dplyr::select("kfactor") %>%
    as.numeric
  reference_kfactor # should be 0.172

  message(paste0("\nReference k-factor ", reference_kfactor, "."))

  # Get k-factor - buffer fold change -------------------------------------------------

  # Subset for buffer and find fold change

  if(nrow(dplyr::filter(kfactors_buffers_data, .data$buffer == buffer_used))>0){

    # Subset by buffer
    subset_data <- kfactors_buffers_data %>%
      dplyr::filter(.data$buffer == buffer_used)
    subset_data

    cat("\nK-factors available for given buffer: \n")
    print(subset_data)

    if(buffer_used == "water" | nrow(dplyr::filter(subset_data, .data$concentration == concentration_used))>0){

      # Subset by concentration
      # For water, concentration value is ignored
      if(buffer_used == "water"){
        buffer_fc_to_use <- subset_data %>%
          dplyr::select("fold_change") %>%
          as.numeric()
        buffer_fc_to_use
      } else {

        # For others:
        buffer_fc_to_use <- subset_data %>%
          dplyr::filter(.data$concentration == concentration_used) %>%
          dplyr::select("fold_change") %>%
          as.numeric()
        buffer_fc_to_use
      }

    } else {

      # Subset by concentration, but conc not found

      # 1. Interpolate with concs we have

      # get data for water as well as buffer
      subset_data_water <- kfactors_buffers_data %>%
        dplyr::filter(.data$buffer == "water")
      subset_data_water
      subset_data_water$concentration <- 0
      subset_data_water
      subset_data # from above
      subset_data <- rbind(subset_data_water, subset_data)
      subset_data

      cat("\nValues used for model (kfactor ~ concentration): \n")
      print(subset_data)

      # build model
      # Fit curve as (Y ~ X):
      model1 <- stats::lm(fold_change ~ concentration, data = subset_data)
      model1

      # print(model1)

      plot1 <- ggplot2::ggplot() +
        ggplot2::geom_point(data = subset_data,
                            ggplot2::aes(x = .data$concentration, y = .data$fold_change)) +
        ggplot2::geom_line(data = subset_data,
                           ggplot2::aes(x = .data$concentration,
                                        y = stats::predict(model1, subset_data))) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          aspect.ratio = 1,
          panel.grid.minor = ggplot2::element_blank()
        )
      plot1
      # plotname <- paste0("kfactor_foldchange_conc_model_", buffer_used, ".pdf")
      # ggplot2::ggsave(file.path(outfolder, plotname),
      #                 plot = plot1,
      #                 width = 8, height = 8, units = "cm")

      # Choose k-factor
      predicted_kfactor_foldchange <- as.numeric(stats::predict(model1, data.frame(concentration = concentration_used)))
      predicted_kfactor_foldchange
      buffer_fc_to_use <- predicted_kfactor_foldchange

    } # concentration_used

  } else {

    # Subset by buffer - but rows don't exist
    message("\nBuffer not in data table:\n")
    print(kfactors_buffers_data)

    message("\nUsing kfactor of water instead...")
    buffer_fc_to_use <- 1
    buffer_fc_to_use

  } # buffer used

  message(paste0("\nChange in k-factor required for given buffer: ", round(buffer_fc_to_use, 3), "."))

  # Get k-factor - temperature fold change -------------------------------------------------

  # Subset for temperature and find fold change

  if(nrow(dplyr::filter(kfactors_temperature_data, .data$temperature == temperature_used))>0){

    # Subset by temperature
    subset_data2 <- kfactors_temperature_data %>%
      dplyr::filter(.data$temperature == temperature_used)
    subset_data2

    cat("\nRows found for temperature: \n")
    print(subset_data2)

    # Find fold change
    temperature_fc_to_use <- subset_data2 %>%
      dplyr::select("fold_change") %>%
      as.numeric()
    temperature_fc_to_use

  } else {

    # Subset by temperature - but rows don't exist
    # message(paste0("\n", temperature_used, " oC not in data."))

    # 1. Interpolate with concs we have

    cat("\nValues used for model (fold_change ~ temperature): \n")
    print(kfactors_temperature_data)

    # build model
    # Fit curve as (Y ~ X):
    model2 <- stats::lm(fold_change ~ temperature, data = kfactors_temperature_data)
    model2

    # print(model2)

    plot1 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = kfactors_temperature_data,
                          ggplot2::aes(x = .data$temperature, y = .data$fold_change)) +
      ggplot2::geom_line(data = kfactors_temperature_data,
                         ggplot2::aes(x = .data$temperature,
                                      y = stats::predict(model2, kfactors_temperature_data))) +
      ggplot2::scale_x_continuous("temperature (oC)") +
      ggplot2::scale_y_continuous("fold change over 25oC", limits = c(0.8,1.2)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        panel.grid.minor = ggplot2::element_blank()
      )
    plot1
    # plotname <- paste0("kfactor_foldchange_conc_model_", buffer_used, ".pdf")
    # ggplot2::ggsave(file.path(outfolder, plotname),
    #                 plot = plot1,
    #                 width = 8, height = 8, units = "cm")

    # Choose k-factor
    predicted_kfactor_foldchange <- as.numeric(stats::predict(model2, data.frame(temperature = temperature_used)))
    predicted_kfactor_foldchange
    temperature_fc_to_use <- predicted_kfactor_foldchange

  } # temperature used

  message(paste0("\nChange in k-factor required for given temperature: ", round(temperature_fc_to_use, 3), "."))

  # Get k-factor - k-factor -------------------------------------------------

  kfactor_to_use <- reference_kfactor * buffer_fc_to_use * temperature_fc_to_use
  kfactor_to_use

  message(paste0("\nOverall k-factor: ", round(kfactor_to_use, 3), "."))

  # Return -------------------------------------------------

  return(kfactor_to_use)
}
