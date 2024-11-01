# plate format functions

#' Find wells of a given plate format
#'
#' Works out all wells in that plate type and returns as a list.
#'
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#' @return a list containing the names of all the wells of a specified plate
#'
#' @export
#' @examples
#'   find_wells(96)
#'   find_wells(24)
#'
find_wells <- function(plate_type = 96){

  # find number of rows and columns
  plate_format <- find_plate_format(plate_type = plate_type)

  # work out row names and column numbers from plate format
  rows <- find_rows(plate_type = plate_type)
  columns <- find_columns(plate_type = plate_type)

  # find combinations of rows and columns, i.e. well names
  df <- expand.grid(row = rows, column = columns)
  df <- df %>%
    dplyr::arrange(row) %>%
    dplyr::mutate(well = paste0(row, .data$column))
  df

  all_wells <- df$well # save list of wells

  return(all_wells)
}

#' Find row and column configuration of a given plate type
#'
#' Internal function.
#'
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#' @keywords internal
find_plate_format <- function(plate_type) {

  ## Print warning if plate type is unusual
  usual_types <- c(6,12,24,48,96,384)
  if(!plate_type %in% usual_types){
    message("Warning: plate type specified is not a standard size. Check outputs carefully.")
  }

  ## Rows and columns are specified by a pair of numbers whose product is n, with the closest in value to one another
  # find first pair starting at sqrt(n) and going down, stop at first
  for (i in floor(sqrt(plate_type)):1) { # loop through all numbers the the first number below sqrt(n)
    if (plate_type %% i == 0) {  # if 'i' (rows) is a divisor of 'n'
      j <- plate_type / i # find j (columns)
      rows_columns <- c(i, j) # save pair
      break() # finish as soon as you find first one - exit for loop
    }
  }

  return(rows_columns)

}

#' Find rows of a given plate type
#'
#' Internal function.
#'
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#' @keywords internal
#' @export
find_rows <- function(plate_type = 96){

  # find number of rows
  plate_format <- find_plate_format(plate_type = plate_type)
  row_number <- plate_format[1]

  # find row names
  rows <- LETTERS[seq(1,row_number)]
  return(rows)

}

#' Find columns of a given plate type
#'
#' Internal function.
#'
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#' @keywords internal
#' @export
find_columns <- function(plate_type = 96){

  # find number of rows
  plate_format <- find_plate_format(plate_type = plate_type)
  column_number <- plate_format[2]

  # find column numbers
  columns <- seq(1,column_number)
  return(columns)

}

#' List of all wells in a specific row of a given plate type
#'
#' @param rows character string, or list of character strings, representing row(s) such as "A" or "H"
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#'
#' @return a list containing the names of all the wells in the specified row(s) of a specified plate
#' @export
#' @examples
#'  find_wells_in_row("A")
#'  find_wells_in_row(c("A","D"), plate_type = 24)
find_wells_in_row <- function(rows,
                              plate_type = 96
                              # put plate_type last so you can call the function with positional arguments only when using the default 96-well plate
){

  # find all wells in plate
  all_wells <- find_wells(plate_type)

  # make pattern
  rows <- toupper(rows) # make them capitals
  row_pattern <- paste(rows, collapse = "|")  # combine letters into an "A|B" style pattern
  pattern <- paste0("^[", row_pattern, "]")

  # subset
  row_wells <- all_wells[grep(pattern, all_wells)] # subset
  return(row_wells)
}

#' List of all wells in a specific column of a given plate type
#'
#' @param columns number or list of numbers, representing column(s) such as `1` or `12`. character strings also accepted.
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#'
#' @return a list containing the names of all the wells in the specified column(s) of a specified plate
#' @export
#' @examples
#'  find_wells_in_column("1")
#'  find_wells_in_column("12")
#'  find_wells_in_column(1)
#'  find_wells_in_column(c(1,3), plate_type = 6)
find_wells_in_column <- function(columns,
                                 plate_type = 96
                                 # put plate_type last so you can call the function with positional arguments only when using the default 96-well plate
){

  # find all wells in plate
  all_wells <- find_wells(plate_type)

  # find row ids # need row ids to make pattern here
  plate_format <- find_plate_format(plate_type)
  row_number <- plate_format[1]
  # first row is always A
  last_row <- LETTERS[row_number] # find last row
  # make row pattern
  row_pattern <- paste0("A-",last_row)

  # make column pattern
  column_pattern <- paste(columns, collapse = "|") # combine numbers into a "1|2" style pattern
  pattern <- paste0("^[", row_pattern, "](", column_pattern, ")$")

  # subset
  column_wells <- all_wells[grep(pattern, all_wells)]
  return(column_wells)
}

#' List of all wells in specific rows and columns of a given plate type
#'
#' @param rows character string, or list of character strings, representing row(s) such as "A" or "H"
#' @param columns number or list of numbers, representing column(s) such as `1` or `12`. character strings also accepted.
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#'
#' @return a list containing the names of all the wells in specified row(s) and column(s) of a specified plate
#' @export
#' @examples
#'  find_wells_in_rowcolumn(c("A","C"), c(2,3))
#'  find_wells_in_rowcolumn(c("A","C"), c(2,3), plate_type = 24)
find_wells_in_rowcolumn <- function(rows, columns,
                                    plate_type = 96
                                    # put plate_type last so you can call the function with positional arguments only when using the default 96-well plate
){

  # find all wells in plate
  all_wells <- find_wells(plate_type)

  # make pattern
  rows <- toupper(rows) # make them capitals
  row_pattern <- paste(rows, collapse = "|")  # combine letters into an "A|B" style pattern
  column_pattern <- paste(columns, collapse = "|")  # combine numbers into a "1|2" style pattern
  pattern <- paste0("^[", row_pattern, "](", column_pattern, ")$")

  # subset
  wells <- all_wells[grep(pattern, all_wells)]

  return(wells)
}
