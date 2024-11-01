#' Get FP concentrations using BCA method
#'
#' Get protein's concentration from a *dilution series* measured with the BCA
#' assay. Takes two input data sets, the BCA assay data and the A562 baseline
#' data. The A562 baseline data is necessary for proteins that might naturally
#' absorb in this range. Script normalises BCA assay data for A562 baseline,
#' then to blank values of BCA. Data from the BSA standards (identified by
#' protein column containing the word "BSA"), are used to construct a standard
#' curve of concentration in ng/ul vs normalised A562 values, which is used to
#' predict FP concentrations at each dilution. The FP's concentration vs
#' dilution values are then used to predict the FP concentration at each
#' dilution in one of two ways. Where option is set to `fit`, a linear model is
#' fitted between FP concentration and dilution, and the fitted values are
#' exported. Where option is set to `highest`, FP concentrations are taken only
#' from the highest concentration/dilution specified and exported.
#'
#' @param microbca_data_csv path of the CSV file of your microBCA data.
#'   required.
#' @param a562_baseline_csv path of the CSV file of your A562 baseline data.
#'   Optional. If data is missing, use NULL. If NULL is specified, the value of
#'   0 is assigned as baseline for all wells. Default is NULL.
#' @param calibr string specifying the value of the 'calibrant' column to assess
#'   with this function. Function subsets the data by the value specified here.
#'   This works by taking all the rows with specified string in the calibrant
#'   column and discarding all other rows (which means that the blanks relevant
#'   to the specified calibrant need to be specified as calibrant = `calibr`,
#'   protein = "none", otherwise they will be removed).
#' @param buffer string specifying the value of the 'media' column to assess
#'   with this function. Function subsets the data by the value specified here.
#' @param protein_seq character string of protein sequence using 1-letter code.
#'   Required for MW calculation.
#' @param option string specifying how to choose the predicted concentration to
#'   use. Default is "highest", in which the mean predicted conc of the highest
#'   dilution (i.e. neat) is used, and is multiplied by the dilution to determine
#'   the concentration of the other dilutions. The alternative, "fit", fits a
#'   `y=mx` linear model and uses that for the concentration determination.
#' @param outfolder path to the folder where output files should be saved.
#'   Defaults to the current working directory.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#'
#' @examples
#' \dontrun{
#'   bca_concs <- get_conc_bca(
#'     microbca_data_csv = "bca_data_parsed.csv", a562_baseline_csv = "a562_data_parsed.csv",
#'     calibr = "mCherry", buffer = "T5N15_pi", protein_seq = protein_seq, option = "highest",
#'     outfolder = "protquant_microbca/mCherry_T5N15pi"
#'   )
#' }
get_conc_bca <- function(microbca_data_csv, a562_baseline_csv = NULL,
                         calibr, buffer, protein_seq,
                         option = "highest", # "highest" or "fit"
                         outfolder = "."){

  # Get data ---------------------------------------------------------------------------------------------------

  bca_data <- utils::read.csv(microbca_data_csv, header = TRUE, check.names = FALSE)
  # Rename data column to something helpful
  col_idx <- which(names(bca_data) == "Raw data")
  names(bca_data)[col_idx] <- "raw_bca"
  # Rename volume column to allow joining
  bca_data <- bca_data %>%
    dplyr::rename(volume_bca = .data$volume)

  if(is.null(a562_baseline_csv)){

    bca_data$raw_baseline <- 0
    bca_data$volume_a562 <- bca_data$volume_bca
    joined_data <- bca_data

  } else {

    a562_baseline_data <- utils::read.csv(a562_baseline_csv, header = TRUE, check.names = FALSE)
    # Rename data column to something helpful
    col_idx <- which(names(a562_baseline_data) == "Raw data")
    names(a562_baseline_data)[col_idx] <- "raw_baseline"
    # Rename volume column to allow joining
    a562_baseline_data <- a562_baseline_data %>%
      dplyr::rename(volume_a562 = .data$volume)
    # # Easiest to join data as early as possible
    # colnames(bca_data)
    # colnames(a562_baseline_data)

    # Low risk join: will def only have one of each well
    joined_data <- dplyr::left_join(bca_data, a562_baseline_data) # keeps all rows of x, but not y
    joined_data
  }

  # Location for saved outputs ---------------------------------------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  # Data analysis ordering ----------------------------------------------------------------

  # One option is to first subset data by protein, then do normalisations to the blanks.
  # The trouble with this is that blanks are typically placed at the end of each row - A12, B12, C12 etc,
  # so each BSA+FP dataset will have a different combination of blanks left in subset_data which affects the blank mean used for the BSA fit!
  # A better option then is to do the normalisations first, then subset later (happens naturally with data_standards/data_tests).

  # Subset data for your buffer of choice ----------------------------------------------------------------

  # Take only the rows with the correct buffer
  # To avoid the large errors that occur when BSA and FP results in different buffers are compared
  joined_data <- joined_data %>%
    dplyr::filter(.data$media == buffer)

  # Normalise --------------------------------------------------------------------------------------------------

  # Step1. For each well, subtract raw_baseline (baseline A562 absorbance) from raw_bca (BCA value)

  norm_data <- joined_data %>%
    dplyr::mutate(base_norm_bca = .data$raw_bca - (.data$volume_bca/.data$volume_a562)*.data$raw_baseline)
  ## nb. volume differences matter because they effect the path length. the path length varies with volume proportionally (see pathlength_water_data).
  norm_data

  # Step2. Take mean of baseline-normalised blanks/buffers, subtract this from baseline-normalised values of proteins

  # Identify negatives/blanks
  blanks_bca <- subset(norm_data, .data$protein == "none")
  blanks_bca
  mean_blanks_bca <- mean(blanks_bca$base_norm_bca, na.rm = TRUE)
  mean_blanks_bca

  # Process data by normalising to blanks
  processed_bca <- norm_data %>%
    dplyr::mutate(mean_blanks_bca = mean_blanks_bca) %>%
    # make new column called 'mean_blanks_bca' which contains the value of object mean_blanks_bca from above
    dplyr::mutate(normalised_bca = .data$base_norm_bca - mean_blanks_bca)
  processed_bca

  # Save processed BCA data - BSA and FP ---------------------------------------------------------------------

  filename <- gsub(".csv", "_processed.csv", basename(microbca_data_csv))
  utils::write.csv(processed_bca, file.path(outfolder, filename), row.names = FALSE)

  # Basic Plots ---------------------------------------------------------------------

  # Remove blanks for plotting
  data.to.plot <- subset(processed_bca, .data$protein != "none" & .data$protein != "")
  data.to.plot

  # Plot BSA and FP(s) - dilution vs raw/normalised values
  plot3 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$raw_bca,
                                     colour = "raw")) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$base_norm_bca,
                                     colour = "baseline normalised")) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$normalised_bca,
                                     colour = "baseline and blank normalised")) +

    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ylab("A562") +
    ggplot2::scale_color_discrete("data") +
    ggplot2::facet_wrap(protein ~ .) +

    ggplot2::labs(title = "microBCA overview") +

    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      # legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      axis.text.x = ggplot2::element_text(angle = 90),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot3
  plotname <- "plot1_BCA_norm.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot3,
                  height = 16, width = 16, units = "cm")

  # Plot BSA vs conc
  plot3 <- ggplot2::ggplot(subset(data.to.plot, .data$protein == "BSA")) +
    ggplot2::geom_point(ggplot2::aes(x = .data$concentration_ngul,
                                     y = .data$normalised_bca)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_y_continuous(limits = c(0, 2.5)) +
    ggplot2::xlab("concentration (ng/ul)") +
    ggplot2::ylab("A562") +
    ggplot2::labs(title = "BSA standards") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot3
  plotname <- "plot2_BSA_conc.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot3,
                  height = 12, width = 12, units = "cm")

  # # Plot BSA vs molecules
  # plot3 <- ggplot2::ggplot(subset(data.to.plot, protein == "BSA")) +
  #   ggplot2::geom_point(ggplot2::aes(x = molecules_microBCA,
  #                                    y = normalised_bca)) +
  #   ggplot2::geom_hline(yintercept = 0) +
  #   ggplot2::scale_y_continuous(limits = c(0, 2.5)) +
  #   ggplot2::xlab("molecules") +
  #   ggplot2::ylab("A562") +
  #   ggplot2::labs(title = "BSA standards") +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(
  #     aspect.ratio = 1,
  #     strip.background = ggplot2::element_blank(),
  #     strip.text = ggplot2::element_text(face = "bold", hjust = 0),
  #     panel.grid.minor = ggplot2::element_blank()
  #   )
  # plot3
  # plotname <- "plot2b_BSA_mols.pdf"
  # ggplot2::ggsave(file.path(outfolder, plotname),
  #                 plot = plot3,
  #                 height = 12, width = 12, units = "cm")

  # Visualising A562 baselines (side question) -----------------------------------------------------------

  # Plot raw A562 baselines...
  plot3 <- ggplot2::ggplot(data.to.plot) +
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$raw_baseline,
                                     colour = .data$protein)) +
    ggplot2::scale_y_continuous(limits = c(0,0.15)) +
    ggplot2::xlab("dilution") +
    ggplot2::ylab("A562") +
    ggplot2::labs(title = "A562 baselines") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot3
  plotname <- "plot3_A562_baselines.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot3,
                  height = 12, width = 12, units = "cm")

  # Normalise...

  # Identify negatives/blanks
  mean_blanks_a562 <- mean(blanks_bca$raw_baseline, na.rm = TRUE)
  mean_blanks_a562
  sd_blanks_a562 <- stats::sd(blanks_bca$raw_baseline, na.rm = TRUE)
  sd_blanks_a562

  # Process data by normalising to blanks
  processed_a562 <- processed_bca %>%
    dplyr::mutate(mean_blanks_baseline = mean_blanks_a562) %>%
    dplyr::mutate(sd_blanks_baseline = sd_blanks_a562) %>%
    dplyr::mutate(normalised_baseline = .data$raw_baseline - .data$mean_blanks_baseline)
  processed_a562

  # Remove blanks for plotting
  data.to.plot3 <- subset(processed_a562, .data$protein != "none" & .data$protein != "")
  data.to.plot3

  # Plot A562 baselines, visualised normalised to its own blanks...
  plot3 <- ggplot2::ggplot(data.to.plot3) +

    # blanks: mean line and sd rect (mean +- 2*sd)
    ggplot2::annotate(geom = "rect",
                      xmin = 0, xmax = 1, ymin = 0-2*sd_blanks_a562, ymax = 0+2*sd_blanks_a562,
                      fill = "grey", alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0) +

    # data
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$normalised_baseline,
                                     colour = .data$protein)) +

    ggplot2::scale_y_continuous(limits = c(-0.02,0.02)) +
    ggplot2::xlab("dilution") +
    ggplot2::ylab("Normalised A562") +
    ggplot2::labs(title = "A562 baselines, normalised to blanks",
                  caption = "shading represents blank mean +/- 2sd") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank()
    )
  plot3
  plotname <- "plot3b_A562_baselines_norm.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot3,
                  height = 12, width = 12, units = "cm")

  # Fit curve to BSA standards --------------------------------------------------------------------------------------------------------

  # Fit polynomial curve to BSA vs concentration relationship, to obtain standard curve -----------------------------------------------

  # Separate out standards
  data_standards <- subset(processed_bca, .data$protein == "BSA")
  data_standards

  # Fit curve with A562 on X and conc on Y:

  # using poly(..., raw = TRUE) gets you coefficients you can use to rebuild the fit as y=ax^2+bx+c
  model <- stats::lm(concentration_ngul ~ poly(normalised_bca, degree = 2, raw = TRUE), data = data_standards)
  model
  # example:
  # Coefficients:
  # (Intercept)  poly(normalised_bca, degree = 2, raw = TRUE)1  poly(normalised_bca, degree = 2, raw = TRUE)2
  # -2.536                                         59.534                                         12.476
  # which means y = 12x^2 + 60x - 3

  # Save coefficients of BSA standard curve model
  bsa_fit_coefs <- data.frame(
    poly2 = as.numeric(model[[1]][3]),
    poly1 = as.numeric(model[[1]][2]),
    intercept = as.numeric(model[[1]][1])
  )
  bsa_fit_coefs
  filename <- gsub(".csv", "_BSA_fit_coeffs.csv", basename(microbca_data_csv))
  utils::write.csv(bsa_fit_coefs, file.path(outfolder, filename), row.names = FALSE)

  # Plot BSA standard curve with polynomial fit
  plot3 <- ggplot2::ggplot(data_standards) +
    ggplot2::geom_point(ggplot2::aes(x = .data$normalised_bca,
                                     y = .data$concentration_ngul)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$normalised_bca,
                                    y = stats::predict(model, data_standards))) +

    ggplot2::xlab("A562") +
    ggplot2::ylab("concentration (ng/ul)") +
    ggplot2::labs(title = "BSA standard curve with polynomial fit") +

    ggplot2::coord_flip() + # reverse x and y!

    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  plot3
  plotname <- "plot6_BSA_curvefit.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plot3,
                  height = 12, width = 12, units = "cm")

  # Get concentration estimates for chosen FP ----------------------------------------------------------------------------------------

  # Get FP data
  data_tests <- subset(processed_bca, .data$protein == calibr)
  data_tests

  # Add column
  # Given model that takes x=normalised_BCA, y=conc, and given normalised_BCA, what is predicted conc?
  data_tests$predicted_conc <- stats::predict(model, data.frame(normalised_bca = data_tests$normalised_bca))
  data_tests

  # Plot concentration estimates on BSA curve
  plt2 <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = data_standards$normalised_bca,
                                     y = data_standards$concentration_ngul)) + # data
    ggplot2::geom_line(ggplot2::aes(x = data_standards$normalised_bca,
                                    y = stats::predict(model, data_standards))) + # model
    ggplot2::geom_point(ggplot2::aes(x = data_tests$normalised_bca,
                                     y = data_tests$predicted_conc),
                        colour = "red") + # data
    ggplot2::xlab("A562") +
    ggplot2::ylab("concentration (ng/ul)") +
    ggplot2::labs(title = "Predicted concentration vs BSA curve") +
    # ggplot2::scale_x_continuous(limits = c(0,0.5)) +
    # ggplot2::scale_y_continuous(limits = c(0,2.5)) +

    ggplot2::coord_flip() + # reverse x and y!

    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    )
  plt2
  plt3 <- plt2 + ggplot2::coord_cartesian(xlim = c(0,0.5), ylim = c(0, 50))
  plt3
  plotname <- "plot7_estimates_vs_curve.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plt3,
                  height = 12, width = 12, units = "cm")

  # Plot concentration estimates vs dilutions
  plt4 <- ggplot2::ggplot(data_tests) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") + # to guide the eyes
    ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                     y = .data$predicted_conc)) + # data
    ggplot2::ylab("predicted concentration") +
    ggplot2::labs(title = "Predicted concentration vs dilution") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid.minor = ggplot2::element_blank(),
    )
  plt4
  plotname <- "plot8_estimates_vs_dilns.pdf"
  ggplot2::ggsave(file.path(outfolder, plotname),
                  plot = plt4,
                  height = 12, width = 12, units = "cm")

  # Complete table ---------------------------------------------------------------------------------------------------

  if(option == "highest"){

    ### Protein concentration prediction decision here is to take the mean of the two highest-dilutions,
    # as higher concentrations will be better estimated by the BCA than lower concentrations.

    # Chosen prediction for highest dilution
    predictedconc_neat <- mean(data_tests$predicted_conc[data_tests$dilution == 1], na.rm = TRUE)
    predictedconc_neat

    ### Complete table

    # Using data_tests as a starting point for a concentration table is no good
    # because this misses out all anomalies (where protein is set to "") and any dilutions not assessed by microBCA (eg col9-11).
    # It should be neater to start from the original BCA data with its data layout, subset for (1) the buffer in question, as for joined_data
    # and (2) for the calibrant in question, similar to what we do with data_tests...
    # ..data_tests subsets for protein == calibr, indicating all wells to be used in microbca calcs
    # ..here we subset for calibrant == calibr, allowing us access to dilutions too low for the microbca assay

    ## Prep base df
    # NB Need to use a base df that contains columns like well, row, column - use input df
    names(bca_data)
    bca_data_tidy <- bca_data %>%
      dplyr::select(c( .data$media:.data$well )) %>% # take all cols between A:B - https://suzan.rbind.io/2018/01/dplyr-tutorial-1/
      dplyr::select(c( -.data$volume_bca )) # remove volume as downstream assays may have used different volumes
    names(bca_data_tidy)
    bca_data_tidy

    # subset
    conc_data <- bca_data_tidy %>%
      dplyr::filter(.data$calibrant == calibr) %>%
      dplyr::filter(.data$media == buffer)
    conc_data

    # complete MW column
    protein_seq
    conc_data$mw_gmol1 <- fpcountr::get_mw(protein = protein_seq)

    # add concs
    conc_data <- conc_data %>%
      # Concentration = predicted concentration for the highest dilution, multiplied down for each dilution
      dplyr::mutate(concentration_ngul = predictedconc_neat * .data$dilution)
    conc_data

    # fill in protein column where dilution exists (for ease of future joining)
    conc_data <- conc_data %>%
      dplyr::mutate(protein = ifelse(.data$protein == "" & !is.na(.data$dilution), calibr, .data$protein))
    # if protein column is empty AND dilution is not NA, add the calibr into the protein column, else, leave protein entry as it is
    conc_data

  }

  if(option == "fit"){

    ## Take mean of each set of replicates
    # Doesn't make a difference whether fit is done between raw or mean values,
    # but mean is used to decide whether a dilution should be trimmed or not
    # as BCA assay is reported to be linear only >2ng/ul.
    data_tests_fit <- data_tests %>%
      dplyr::group_by(.data$dilution) %>%
      dplyr::mutate(mean_predicted_conc = mean(.data$predicted_conc, na.rm = TRUE))
    data_tests_fit

    plt4 <- ggplot2::ggplot(data_tests_fit) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey") + # to guide the eyes
      ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                       y = .data$mean_predicted_conc)) + # data
      ggplot2::scale_x_continuous(limits = c(0,1)) +
      ggplot2::ylab("predicted concentration") +
      ggplot2::labs(title = "Predicted concentration (mean) vs dilution") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        panel.grid.minor = ggplot2::element_blank(),
      )
    plt4

    # Remove anomalous values
    data_tests_fit_trimmed <- data_tests_fit %>%
      dplyr::filter(.data$mean_predicted_conc > 1)
    data_tests_fit_trimmed

    plt4 <- ggplot2::ggplot(data_tests_fit_trimmed) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey") + # to guide the eyes
      ggplot2::geom_point(ggplot2::aes(x = .data$dilution,
                                       y = .data$mean_predicted_conc)) + # data
      ggplot2::scale_x_continuous(limits = c(0,1)) +
      ggplot2::ylab("predicted concentration") +
      ggplot2::labs(title = "Predicted concentration (mean, trimmed) vs dilution") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        panel.grid.minor = ggplot2::element_blank(),
      )
    plt4

    message("There are ", length(unique(data_tests_fit_trimmed$dilution)), " dilutions for this calibrant whose predicted concentration is over 1ng/ul.")
    message("Using ", length(unique(data_tests_fit_trimmed$dilution)), " data points for the fit.")

    # Now do second linear fit. y=mx -----------------------------------------------------------

    model3 <- stats::lm(predicted_conc ~ dilution + 0, data = data_tests_fit_trimmed)
    model3
    # Coefficients:
    #   dilution
    # 22.11

    # Save fit
    conc_fit_coefs <- data.frame(
      m = as.numeric(model3[[1]][1])
    )
    conc_fit_coefs

    # Report fit
    message("Fit coefficients: y = (", signif(conc_fit_coefs$m, 3), ")*x")

    # Plot fit
    plt4 <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, colour = "grey") + # to guide the eyes
      ggplot2::geom_point(ggplot2::aes(x = data_tests_fit$dilution,
                                       y = data_tests_fit$predicted_conc),
                          colour = "black") + # all data
      ggplot2::geom_point(ggplot2::aes(x = data_tests_fit_trimmed$dilution,
                                       y = data_tests_fit_trimmed$predicted_conc),
                          colour = "red") + # trimmed data used for the fit
      ggplot2::geom_line(ggplot2::aes(x = data_tests_fit$dilution, # extrapolated to all values of dilution tested
                                      y = stats::predict(model3, data_tests_fit)),
                         colour = "red") +

      ggplot2::scale_x_continuous("dilution", limits = c(0,1)) +
      ggplot2::ylab("predicted concentration") +
      ggplot2::labs(title = "Predicted concentration vs dilution, with fit",
                    caption = "red: mean values used for the fit; black: all points") +

      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        panel.grid.minor = ggplot2::element_blank(),
      )
    plt4
    plotname <- "plot8b_estimates_vs_dilns_withfit.pdf"
    ggplot2::ggsave(file.path(outfolder, plotname),
                    plot = plt4,
                    height = 12, width = 12, units = "cm")

    ### Protein concentration prediction:

    # take the gradient of the linear fit as the 'conversion' coefficient between the dilution and the predicted concentration
    conc_fit_coefs$m

    ### Complete table

    # Using data_tests as a starting point for a concentration table is no good
    # because this misses out all anomalies (where protein is set to "") and any dilutions not assessed by microBCA (eg col9-11).
    # It should be neater to start from the original BCA data with its data layout, subset for (1) the buffer in question, as for joined_data
    # and (2) for the calibrant in question, similar to what we do with data_tests...
    # ..data_tests subsets for protein == calibr, indicating all wells to be used in microbca calcs
    # ..here we subset for calibrant == calibr, allowing us access to dilutions too low for the microbca assay

    ## Prep base df
    # NB Need to use a base df that contains columns like well, row, column - use input df
    names(bca_data)
    bca_data_tidy <- bca_data %>%
      dplyr::select(c( .data$media:.data$well )) %>% # take all cols between A:B - https://suzan.rbind.io/2018/01/dplyr-tutorial-1/
      dplyr::select(c( -.data$volume_bca )) # remove volume as downstream assays may have used different volumes
    names(bca_data_tidy)
    bca_data_tidy

    # subset
    conc_data <- bca_data_tidy %>%
      dplyr::filter(.data$calibrant == calibr) %>%
      dplyr::filter(.data$media == buffer)
    conc_data

    # complete MW column
    protein_seq
    conc_data$mw_gmol1 <- fpcountr::get_mw(protein = protein_seq)

    # add concs
    conc_data <- conc_data %>%
      # Concentration = predicted concentration for the highest dilution, multiplied down for each dilution
      dplyr::mutate(concentration_ngul = conc_fit_coefs$m * .data$dilution)
    conc_data

    # fill in protein column where dilution exists (for ease of future joining)
    conc_data <- conc_data %>%
      dplyr::mutate(protein = ifelse(.data$protein == "" & !is.na(.data$dilution), calibr, .data$protein))
    # if protein column is empty AND dilution is not NA, add the calibr into the protein column, else, leave protein entry as it is
    conc_data

  }

  # Save and Return -----------

  # conc_data - predictions for calibrations
  # FP only, all dilutions, no BCA data columns
  filename <- "protein_concs_bca.csv"
  utils::write.csv(conc_data, file.path(outfolder, filename), row.names = FALSE)

  # data_tests - predictions for each well individually
  # only microbca dilutions, includes all BCA calculation columns and individual conc predictions in last column
  filename <- gsub(".csv", "_predicted_concs.csv", basename(microbca_data_csv))
  utils::write.csv(data_tests, file.path(outfolder, filename), row.names = FALSE)

  # -----------------------------------------------------------------------------

  return(conc_data)

}
