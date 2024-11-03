#' Get path length in a typical 96-well plate
#'
#' Calculates the path length of an input `test_volume` by comparing it to a
#' linear model using internal data of pathlengths vs volume. ... To work out
#' the pathlength in a microplate assay, there are two options. Option 1:
#' Measure the pathlength each time you do an assay by taking A975-A900 for each
#' well as the `kfactor(well)`, and calculating `pathlength` = `kfactor(well)` /
#' `kfactor(1cm pathlength)`. The 1cm k-factors can be derived using
#' `get_kfactor()` function. Option 2 (used here): Create a standardised dataset
#' of pathlength vs volume by measuring the A975-A900 of defined volumes of
#' buffer, and calculating the pathlengths of each (the `pathlength_water_data`
#' dataset) and use this to create a linear model of pathlength vs volume.
#' Experimental data suggests that aqueous buffers give very similar pathlength
#' values for a given volume, so here we use data from water to approximate all
#' aqueous buffers.
#'
#' @param test_volume numeric value of volume whose pathlength is required, in
#'   microlitres (ul)
#' @param plot logical. Should the function plot the model and prediction?
#'   Defaults to FALSE.
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @export
#'
#' @examples get_pathlength(200)
get_pathlength <- function(test_volume, plot = FALSE, outfolder = "."){

  # Location for saved outputs -------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # Get volume data -------------------------------------------------

  pathlength_volume_data <- fpcountr::pathlength_water_data
  pathlength_volume_data

  data.to.plot <- pathlength_volume_data
  plot1 <- ggplot2::ggplot(data = data.to.plot) +
    ggplot2::geom_point(ggplot2::aes(x = volume, y = .data$pathlength)) +
    ggplot2::scale_x_continuous("volume (ul)", limits = c(0,300)) +
    ggplot2::scale_y_continuous("path length (cm)", limits = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(angle= 90),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot1

  # Get pathlength -------------------------------------------------

  # cat("\nValues used for model (pathlength ~ volume): \n")
  # print(pathlength_volume_data)

  # Make pathlength ~ volume model
  model2 <- stats::lm(pathlength ~ volume, data = pathlength_volume_data)
  model2

  # Predict path length based on volume that was used
  predicted_pathlength <- as.numeric(stats::predict(model2, data.frame(volume = test_volume)))
  predicted_pathlength

  # Plot -------------------------------------------------

  if(plot){

    small_df <- data.frame(volume = test_volume,
                           pathlength = predicted_pathlength)

    # Make prediction table from model to plot a wider range to volume = 0
    df2 <- data.frame(volume <- seq(0,300,25)) # create list of potential volumes
    names(df2) <- c("volume")
    df2
    # predict the path length
    df2 <- df2 %>%
      dplyr::mutate(predicted_pathlength = stats::predict(model2, data.frame(volume = df2$volume)))
    df2

    # plot3 path length
    plot1 <- ggplot2::ggplot() +

      # black: pathlength vs volume data
      ggplot2::geom_point(data = pathlength_volume_data,
                          ggplot2::aes(x = volume, y = .data$pathlength)) +
      # black: pathlength vs volume model using volumes that go to 0
      ggplot2::geom_line(data = df2,
                         ggplot2::aes(x = volume,
                                      y = stats::predict(model2, df2))) +
      # red: predicted pathlength for given volume
      ggplot2::geom_point(data = small_df,
                          ggplot2::aes(x = volume, y = .data$pathlength),
                          colour = "red", shape = 1, size = 2) +
      # ggplot2::geom_vline(xintercept = small_df$volume[1], colour = "red") +
      # ggplot2::geom_hline(yintercept = small_df$pathlength[1], colour = "red") +

      ggplot2::scale_x_continuous("volume (ul)", limits = c(0,300)) +
      ggplot2::scale_y_continuous("path length (cm)", limits = c(0, 1)) +

      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        axis.text.x = ggplot2::element_text(angle= 90),
        panel.grid.minor = ggplot2::element_blank()
      )
    plot1
    plotname <- paste0("pathlength_volume_model_prediction.pdf")
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = plot1,
                    width = 8, height = 8, units = "cm")
  }

  # Return -------------------------------------------------

  return(predicted_pathlength)
}
