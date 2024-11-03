#' Parse Tecan Spark/Magellan data from absorbance spectrum data
#'
#' Parses absorbance spectrum data exported from a Tecan Spark plate reader
#' using Magellan software. Parsing consists of data extraction, data tidying,
#' and data joining to relevant 'plate layout' metadata.
#'
#' @param data_csv path to CSV file from Tecan Spark plate reader
#' @param layout_csv path to CSV file containing plate layout information
#' @param wellstart character string representing first well recorded in data
#'   file. Defaults to "A1".
#' @param wellend character string representing last well recorded in data file.
#'   Defaults to "H12".
#'
#' @return a data.frame containing the parsed plate reader data
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   parsed_data <- parse_magellan_spectrum(
#'     data_csv = "data/20210104_data.csv",
#'     layout_csv = "data/20210104_data_layout.csv",
#'     wellstart = "A1", wellend = "H12"
#'   )
#' }
parse_magellan_spectrum <- function(data_csv, layout_csv,
                                    wellstart = "A1", wellend = "H12"
) {

  data <- utils::read.table(data_csv, sep = ",", blank.lines.skip = TRUE,
                            header = FALSE, stringsAsFactors = FALSE)

  plate_layout <- utils::read.csv(layout_csv)

  # Spectrum data ----------------------------------------------------------

  ## Work out the number of rows of data
  # find first row
  data_start <- which(data[, 1] == wellstart)
  # find last row
  data_end <- which(data[, 1] == wellend)
  total_rows <- data_end - data_start +1 # total number of rows of data is last row-first row, +1 to be inclusive of first row
  message(paste0(total_rows, " wells identified."))

  ## Get data
  # Grab the list of wells
  data_wells <- data[data_start:data_end,1] # from first data row to last data row, col1
  # Grab the data
  data_subset <- data[data_start:data_end,2:ncol(data)] # from first data row to last data row, from col2 to last column
  # str(data_subset)
  # Convert columns to numeric (which from spark magellan are character when they contain "overflow" entries)
  data_numeric <- data_subset %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.numeric))
  # str(data_numeric)
  # Join well list and data
  data_joined <- cbind(data_wells, data_numeric)
  # Add column names from row2 of original data
  names(data_joined) <- data[2,1:ncol(data)]
  names(data_joined)[1] <- "well"

  ## Join plate_layout to data
  joined_data <- dplyr::left_join(x = plate_layout, y = data_joined, by = "well")

  # rearrange data ----------------------------------------------------------
  joined_data$row <- substr(x = joined_data$well, start = 1, stop = 1)
  joined_data$column <- as.numeric(substr(x = joined_data$well, start = 2,
                                          stop = nchar(joined_data$well)))
  joined_data <- dplyr::arrange_at(joined_data, dplyr::vars(.data$row,
                                                            .data$column))
  joined_data

  # write parsed data to csv ------------------------------------------------
  out_name <- gsub(".csv", "_parsed.csv", data_csv)
  utils::write.csv(x = joined_data, file = out_name, row.names = FALSE)

  return(joined_data)

}
