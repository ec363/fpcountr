#' Parse Tecan Spark/SparkControl data from absorbance spectrum data
#'
#' Parses absorbance spectrum data exported from a Tecan Spark plate reader
#' using SparkControl software. Parsing consists of data extraction, data
#' tidying, and data joining to relevant metadata. Based on
#' `parse_magellan_spectrum()`.
#'
#' @param data_csv path to CSV file from Tecan Spark plate reader
#' @param metadata_csv path to CSV file containing metadata
#' @param wellstart character string representing first well recorded in data
#'   file. Defaults to "A1".
#' @param wellend character string representing last well recorded in data file.
#'   Defaults to "H12".
#' @param save_file logical. Would you like to save the parsed file as a CSV?
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   parsed_data <- parse_sparkcontrol_spectrum(
#'     data_csv = "data/20210104_data.csv",
#'     metadata_csv = "data/20210104_metadata.csv",
#'     wellstart = "A1", wellend = "H12",
#'     save_file = TRUE
#'   )
#' }
parse_sparkcontrol_spectrum <- function(data_csv, metadata_csv, wellstart = "A1", wellend = "H12", save_file = FALSE){

  data <- utils::read.table(data_csv, sep = ",", blank.lines.skip = TRUE,
                            header = FALSE, stringsAsFactors = FALSE)

  metadata <- utils::read.csv(metadata_csv)

  # Spectrum data ----------------------------------------------------------

  ## Work out the number of rows of data
  # find first row
  data_start <- which(data[, 1] == "Wavel.") + 1
  # find last row
  data_end <- which(data[, 1] == "End Time")
  data_end <- data_end[length(data_end)] - 2 # pick the last one, - 2
  total_rows <- data_end - data_start + 1 # total number of rows of data is last row-first row, +1 to be inclusive of first row
  message(paste0(total_rows, " wells identified."))

  ## Grab the list of wells
  first_well <- which(data[(data_start-1), ] == wellstart)
  last_well <- which(data[(data_start-1), ] == wellend)
  data_wells <- as.character(data[(data_start-1),(first_well:last_well)])
  # data_wells

  ## Grab the data
  data_subset <- data[data_start:data_end,2:ncol(data)] # from first data row to last data row, from col2 to last column
  # str(data_subset)
  # Convert columns to numeric (they are grabbed as character originally)
  data_numeric <- data_subset %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.numeric))
  # str(data_numeric)
  # Flip matrix (from rows = wavelengths and cols = wells, to rows = wells and cols = wavelengths)
  data_flipped <- as.data.frame(t(as.matrix(data_numeric)))

  # Add wavelength column names from col1 of original data
  names(data_flipped) <- data[(data_start:data_end),1]

  # Join well list and data
  data_joined <- cbind(data_wells, data_flipped)
  names(data_joined)[1] <- "well"
  # data_joined[1:6,1:10]

  ## Join metadata to data
  joined_data <- dplyr::left_join(x = metadata, y = data_joined, by = "well")
  # joined_data[1:6,1:20]

  # rearrange data ----------------------------------------------------------

  # add row and data columns at the end
  joined_data <- joined_data %>%
    dplyr::mutate(row = substr(x = .data$well, start = 1, stop = 1)) %>% # extract first element of well. this is row.
    dplyr::mutate(column = as.numeric(substr(x = .data$well, start = 2, stop = nchar(.data$well)))) %>% # extract column number.
    dplyr::arrange(dplyr::across(c(.data$row, .data$column))) # sort rows by row > column

  # write parsed data to csv ------------------------------------------------
  if(save_file){
    out_name <- gsub(".csv", "_parsed.csv", data_csv)
    utils::write.csv(x = joined_data, file = out_name, row.names = FALSE)
  }

  return(joined_data)

}
