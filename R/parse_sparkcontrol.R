#' Parse Tecan Spark/SparkControl data from endpoint and kinetic assays
#'
#' Parses data exported from a Tecan Spark plate reader using SparkControl
#' software. Handles standard (endpoint) or timecourse (kinetic) data containing
#' absorbance and/or fluorescence readings, but cannot handle spectra, such as
#' absorbance scans. Parsing consists of data extraction, data tidying, and data
#' joining to relevant metadata. Originally based on `flopr::spark_parse()`
#' function.
#'
#' @param data_csv path to CSV file from Tecan Spark plate reader
#' @param metadata_csv path to CSV file containing metadata
#' @param timeseries logical. Is the data a timeseries? Defaults to FALSE.
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   parsed_calib_plate <- parse_sparkcontrol(
#'     data_csv = "calibrations/20210104_calibration_data.csv",
#'     metadata_csv = "calibrations/20210104_calibration_metadata.csv",
#'     timeseries = FALSE
#'   )
#'
#'   parsed_data <- parse_sparkcontrol(
#'     data_csv = "data/20210104_data.csv",
#'     metadata_csv = "data/20210104_metadata.csv",
#'     timeseries = TRUE
#'   )
#' }
parse_sparkcontrol <- function(data_csv, metadata_csv, timeseries = FALSE) {

  data <- utils::read.table(data_csv, sep = ",", blank.lines.skip = TRUE,
                            header = FALSE, stringsAsFactors = FALSE)

  metadata <- utils::read.csv(metadata_csv)

  if(timeseries == TRUE){

    ## find start of data section
    start_time_idx <- which(data[, 1] == "Start Time")
    if (length(start_time_idx) > 1) { # if there are multiple Start Time wells
      start_time_idx <- start_time_idx[length(start_time_idx)] # pick the last one
    }

    next_block_start_idx <- start_time_idx + 2

    end_of_file <- FALSE

    all_data <- c()

    while (!end_of_file) {

      ## find what is being measured
      block_name <- data[next_block_start_idx, 1]

      ## check if we are at the end of the file
      if (block_name == "End Time") {
        end_of_file <- TRUE
        break
      }

      ## find where the end of the current measurement block is
      block_end_idx <- next_block_start_idx
      while (data[block_end_idx, 1] != "") {
        block_end_idx <- block_end_idx + 1
      }

      ## grab the data only for that measurement
      new_block <- data[(next_block_start_idx + 2):(block_end_idx - 1), ]
      colnames(new_block) <- data[(next_block_start_idx + 1),]
      # new_block[1:6,1:10] # check

      ## trim unnecessary variables
      new_block <- new_block[,-c(1,3)]
      # new_block[1:6,c(1:10)] # check

      ## fix time column
      # rename time column
      colnames(new_block)[1] <- "time" # new - rename time column
      # convert time to minutes, round to nearest minute # new
      new_block <- new_block %>%
        dplyr::mutate(time = as.numeric(.data$time)) %>% # change to numeric else division fails
        dplyr::mutate(time = round(.data$time/60, 1))
      # new_block[1:6,c(1:10)] # check

      ## create wells column
      # pivot table to create 'wells' column, which is needed to do the join w the metadata
      new_block <- new_block %>%
        tidyr::pivot_longer(cols = c(2:ncol(new_block)),
                            names_to = "well",
                            values_to = "value") %>%
        dplyr::mutate(value = as.numeric(.data$value)) # ensure all data columns are numeric
      # new_block[1:6,] # check

      ## add metadata for each well
      joined_block <- dplyr::full_join(metadata, new_block)
      # utils::head(joined_block)

      ## append name of reading/measure
      joined_block$measure <- block_name
      # joined_block[1:6,c(1:3,ncol(joined_block)-2,ncol(joined_block)-1,ncol(joined_block))] # check

      all_data <- rbind(all_data, joined_block)
      # all_data[1:6,c(1:3,ncol(all_data)-2,ncol(all_data)-1,ncol(all_data))] # check

      next_block_start_idx <- block_end_idx + 1
    }

    # rearrange data ----------------------------------------------------------

    ## create separate column for each of the measures/readings
    wide_data <- all_data %>%
      tidyr::pivot_wider(names_from = .data$measure, values_from = .data$value) %>% # pivot to create new columns for each reading/measure
      dplyr::mutate(row = substr(x = .data$well, start = 1, stop = 1)) %>% # extract first element of well. this is row.
      dplyr::mutate(column = as.numeric(substr(x = .data$well, start = 2, stop = nchar(.data$well)))) %>% # extract column number.
      dplyr::arrange(dplyr::across(c(.data$row, .data$column))) # sort rows by row > column

    # write parsed data to csv ------------------------------------------------
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = wide_data, file = out_name, row.names = FALSE)

    return(wide_data)
  } # timeseries true

  if (timeseries == FALSE){

    message("Standard/Endpoint data from SparkControl is assumed to always be exported in 8x12 matrix format.")

    ## find location of each reading
    block_start_idx <- which(data[, 1] == "<>") + 1
    block_start_idx

    ## find location of each reading/measurement name
    names_idx <- which(data[, 1] == "Name")
    names_idx <- names_idx[-1] # skip first
    names_idx

    all_data <- c()

    ## for each reading...
    for(i in 1:length(block_start_idx)){

      ## record reading name
      block_name <- data[names_idx[i], 2]
      block_name

      ## grab the data
      new_block <- data[block_start_idx[i]:(block_start_idx[i]+7), 2:13] # get data - 8x12
      new_block <- new_block %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), as.numeric)) # make sure data is all numeric
      new_block

      ## convert matrix format to tidy format, assuming we have a 8x12 matrix # sparkcontrol
      # assign row and column names
      rownames(new_block) <- c(LETTERS[1:8])
      colnames(new_block) <- c(1:12)
      new_block
      # convert matrix to dataframe with rownames preserved
      new_block <- new_block %>%
        tibble::rownames_to_column(var = "row") %>%
        tidyr::pivot_longer(cols = -row, names_to = "column", values_to = "value") %>%
        dplyr::mutate(well = paste0(.data$row, .data$column)) %>%
        dplyr::select(.data$well, .data$value)
      new_block

      ## join data and metadata
      joined_block <- dplyr::full_join(metadata, new_block)
      joined_block$measure <- block_name # add reading name
      # utils::head(joined_block)

      ## combine with all data
      all_data <- rbind(all_data, joined_block)

    } # for each reading

    # rearrange data ----------------------------------------------------------

    # create separate column for each of the measures/readings
    wide_data <- all_data %>%
      tidyr::pivot_wider(names_from = .data$measure, values_from = .data$value) %>% # pivot to create new columns for each reading/measure
      dplyr::mutate(row = substr(x = .data$well, start = 1, stop = 1)) %>% # extract first element of well. this is row.
      dplyr::mutate(column = as.numeric(substr(x = .data$well, start = 2, stop = nchar(.data$well)))) %>% # extract column number.
      dplyr::arrange(dplyr::across(c(.data$row, .data$column))) # sort rows by row > column

    # write parsed data to csv ------------------------------------------------
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = wide_data, file = out_name, row.names = FALSE)

    return(wide_data)
  } # timeseries false

}
