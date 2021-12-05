#' Normalise fluorescence against negative well autofluorescence
#'
#' Used by `process_plate` function for fluorescence normalisation. Remains
#' virtually unchanged from `flopr::flu_norm`, except that `af_model` can be set
#' to `NULL` to normalise to blank wells instead, and files are saved to a
#' specified `outfolder`.
#'
#' @param pr_data a long data.frame containing you plate reader data with OD
#' normalised
#' @param neg_well the well coordinates of a non-fluorescent control
#' @param blank_well the well coordinates of a media blank
#' @param flu_name the column name of the fluorescence chanel to normalise
#' @param af_model model used to fit negative control autofluorescence.
#' For now these include "polynomial", "invers_poly", "exponential", "spline" or "loess".
#' If set to NULL, no model is used, and fluorescence is normalised akin to OD: by subtracting the value for the blanks.
#' @param data_csv path to the original data. Used for saving normalisation curve plots.
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @export
#'
#' @return
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
flu_norm <- function(pr_data, neg_well, blank_well, flu_name, af_model, data_csv,
                     outfolder) {

  ### function 'normalises' fluor raw values to autofluorescence model
  # af_model used here
  # model predicted to explain relationship between normalised_OD and GFP in neg wells (autofluorescence)

  pr_data$v1 <- pr_data[, flu_name] # copies fluor data into new column
  # head(pr_data[, "GFP"]) # raw GFP
  # head(pr_data[, "v1"]) # copied over GFP

  # option for af_model to be NULL --------------------------------------

  if(is.null(af_model)){

    # when af_model is NULL:
    # normalise as OD is normalised (to the fluorescence of the blank media wells)

    # remove background optical density signal -------------------------------------

    pr_data <- pr_data %>%
      dplyr::group_by(.data$time) %>%
      dplyr::mutate(v1 = .data$v1 - mean(.data$v1[.data$well %in% blank_well]))
    # for each timepoint, normalises GFP to mean of blank wells
    # normalised_GFP = (copied over GFP)- mean(GFP in blank wells)
    # head(pr_data[, "GFP"], 12) # raw GFP data
    # head(pr_data[, "v1"], 12) # normalised GFP data

    names(pr_data)[ncol(pr_data)] <- paste0("normalised_", flu_name)
    # head(pr_data[, paste("normalised_", flu_name, sep = "")], 12) # normalised GFP data

    return(as.data.frame(pr_data))

  }

  # fit autofluorescence model to negative control --------------------------

  negative_data <- pr_data %>% dplyr::filter(.data$well %in% neg_well)

  if (af_model == "polynomial") {
    model <- stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ 2 + c),
                        start = c(a = 1, b = 1, c = 1), data = negative_data)

  } else if (af_model == "inverse_poly") {
    model <- stats::nls(normalised_OD ~ (a * v1 + b * v1 ^ 2 + c),
                        start = c(a = 1, b = 1, c = 1), data = negative_data)

  } else if (af_model == "exponential") {
    ## ae^(bx) + c
    ## intial parameter estimation
    model_0 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]),
                  b = stats::coef(model_0)[2],
                  c = -1)

    model <- stats::nls(v1 ~ (a * exp(b * normalised_OD) + c),
                        start = start, data = negative_data)

  } else if (af_model == "bi_exponential") {
    ## a exp(bx) + c exp(dx) + e

    model_0 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]*0.2),
                  b = stats::coef(model_0)[2]*0.2,
                  c = exp(stats::coef(model_0)[1])*0.8,
                  d = stats::coef(model_0)[2]*0.8,
                  e = 1)

    model <- stats::nls(v1 ~ (a * exp(b * normalised_OD) +
                                c * exp(d * normalised_OD) + e),
                        start = start,
                        data = negative_data)

  } else if (af_model == "linear_exponential") {
    ## ax + be^cx + d

    model_01 <- stats::lm(v1 ~ normalised_OD, data = negative_data)
    model_02 <- stats::lm(log(v1) ~ normalised_OD, data = negative_data)
    start <- list(a = stats::coef(model_01)[2],
                  b = exp(stats::coef(model_02)[1]),
                  c = stats::coef(model_02)[2],
                  d = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD +
                                b * exp(c * normalised_OD) + d),
                        start = start,
                        data = negative_data)

  } else if (af_model == "power") {
    ## ax^b + c
    model_0 <- stats::lm(log(v1) ~ log(normalised_OD), data = negative_data)
    start <- list(a = exp(stats::coef(model_0)[1]),
                  b = stats::coef(model_0)[2],
                  c = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD ^ b + c),
                        start = start, data = negative_data)

  } else if (af_model == "linear_power") {
    ## ax + bx^c + d
    model_01 <- stats::lm(v1 ~ normalised_OD, data = negative_data)
    model_02 <- stats::lm(log(v1) ~ log(normalised_OD), data = negative_data)
    start <- list(a = stats::coef(model_01)[2],
                  b = exp(stats::coef(model_02)[1]),
                  c = stats::coef(model_02)[2],
                  d = 1)

    model <- stats::nls(v1 ~ (a * normalised_OD + b * normalised_OD ^ c + d),
                        start = start,
                        data = negative_data)
  } else if (af_model == "loess") {
    model <- stats::loess(v1 ~ normalised_OD,
                          data = negative_data,
                          span = 0.5)
  } else if (af_model == "spline") {
    model <- mgcv::gam(v1 ~ s(normalised_OD), data = negative_data)
  }

  # plot model fit curves ---------------------------------------------------

  if (af_model == "polynomial" | af_model == "power" |
      af_model == "exponential" | af_model == "bi_exponential" |
      af_model == "linear_exponential" | af_model == "linear_power" |
      af_model == "loess") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = stats::predict(model,
                                                         negative_data))) +
      ggplot2::geom_point(ggplot2::aes(x = negative_data$normalised_OD,
                                       y = negative_data$v1)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  } else if (af_model == "spline") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = mgcv::predict.gam(model, negative_data))) +
      ggplot2::geom_point(ggplot2::aes(x = negative_data$normalised_OD,
                                       y = negative_data$v1)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  }
  else if (af_model == "inverse_poly") {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = negative_data$normalised_OD,
                                      y = ((- (stats::coef(model)[1]) +
                                              sqrt((stats::coef(model)[1]) ^ 2 -
                                                     4 *
                                                     (stats::coef(model)[2]) *
                                                     (stats::coef(model)[3]) +
                                                     4 *
                                                     (stats::coef(model)[2]) *
                                                     negative_data$normalised_OD)) /
                                             (2 * (stats::coef(model)[2]))))) +
      ggplot2::geom_point(ggplot2::aes(y = negative_data$v1,
                                       x = negative_data$normalised_OD)) +
      ggplot2::scale_x_continuous("normalised_OD") +
      ggplot2::scale_y_continuous(flu_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  }

  plotname <- gsub(".csv", paste("_norm-curve_", flu_name, ".pdf", sep = ""), basename(data_csv)) # changed: plot location
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plt,
                  height = 16, width = 24, units = "cm")

  # normalise fluorescence data ---------------------------------------------

  if (af_model == "polynomial" | af_model == "power" |
      af_model == "exponential" | af_model == "bi_exponential" |
      af_model == "linear_exponential" | af_model == "linear_power" |
      af_model == "loess") {

    pr_data$v1 <- pr_data$v1 - stats::predict(model, pr_data)
    # normalised fluorescence using the fitted autofluorescence model values, which act as the expected baseline fluorescence here.

  } else if (af_model == "spline") {
    pr_data$v1 <- pr_data$v1 - mgcv::predict.gam(model, pr_data)
    # normalised fluorescence using the fitted autofluorescence model values, which act as the expected baseline fluorescence here.
  }
  else if (af_model == "inverse_poly") {
    pr_data$v1 <- pr_data$v1 - ((- (stats::coef(model)[1]) +
                                   sqrt((stats::coef(model)[1]) ^ 2 - 4 *
                                          (stats::coef(model)[2]) *
                                          (stats::coef(model)[3]) +
                                          4 * (stats::coef(model)[2]) *
                                          pr_data$normalised_OD)) /
                                  (2 * (stats::coef(model)[2])))
    # normalised fluorescence using the fitted autofluorescence model values, which act as the expected baseline fluorescence here.
  }

  # rename normalised fluorescence column
  names(pr_data)[ncol(pr_data)] <- paste("normalised_", flu_name, sep = "")

  return(pr_data)
}
