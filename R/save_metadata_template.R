#' Save a template metadata file
#'
#' Save a CSV file containing the minimum columns required by the functions that
#' parse and process the given `data_type` in `fpcountr`. One of three data
#' types can be chosen: absorbance spectrum data for calibrations
#' ("absspectrum"), fluorescence data for calibrations ("fluordata") and
#' experimental data ("exptdata") The `well` column is populated with all the
#' wells in a given plate type specified by `plate_type`. The file is saved to
#' the location specified by `outfolder`.
#'
#' @param data_type type of data, from absorbance spectrum data for calibrations
#'   ("absspectrum"), fluorescence data for calibrations ("fluordata") and
#'   experimental data ("exptdata").
#' @param plate_type type of plate. numeric, i.e. `96` for 96-well plate.
#' @param outfolder path to folder where output files should be saved.
#'
#' @importFrom dplyr %>%
#'
#' @export
#' @return No return value, called for side effects.
#' @examples
#' \dontrun{
#'   save_metadata_template(data_type = "exptdata", plate_type = 96, outfolder = ".")
#' }
save_metadata_template <- function(
    data_type, # absspectrum, fluordata, exptdata
    plate_type = 96,
    outfolder = ""
){

  # input checks
  if(!(data_type %in% c("absspectrum", "fluordata", "exptdata"))){
    message('Error: Please specify a valid data_type: "absspectrum", "fluordata" or "exptdata".')
    return()
  }
  if(!dir.exists(outfolder)){
    message("Error: Please specify a valid path for the location 'outfolder' where the file should be saved.")
    return()
  }

  message("Loading minimal template for chosen data type...")
  message("Plate type set to: ", plate_type, "-well plate...")

  # pick data ----

  # v1 load data from vignette folder
  # data_file <- file.path("vignettes", "data", paste0("vignetteparsing_meta_", data_type, ".csv")) # failed to find
  # data_file <- file.path("inst", "extdata", paste0("vignetteparsing_meta_", data_type, ".csv")) # failed to find
  # df <- read.csv(data_file)

  # v2 load data from package's internal data (in R/sysdata.rda, created with usethis::use_data(..internal = TRUE))
  # NB. INTERNAL DATA FILES LIVE IN R/sysdata.rda !
  if(data_type == "absspectrum"){df <- meta_absspectrum}
  if(data_type == "fluordata"){df <- meta_fluordata}
  if(data_type == "exptdata"){df <- meta_exptdata}

  # find filename ----
  filename <- paste0("template_metadata_for_", data_type, "_data.csv")
  filename

  # build template df ----

  # grab column names from metadata requirement file
  column_names <- df %>%
    dplyr::select(.data$metadata_variables_expected) %>%
    dplyr::filter(.data$metadata_variables_expected != "") %>% # remove empties
    dplyr::filter(.data$metadata_variables_expected != "row" & .data$metadata_variables_expected != "column") %>% # remove "row" and "column"
    dplyr::distinct() %>%
    as.vector()
  column_names # list
  column_names <- unlist(column_names) # required to make df the right no of columns
  column_names # char array

  # build list of wells
  wells_list <- find_wells(plate_type = plate_type)
  wells_list

  # make df
  template_df <- data.frame(matrix(ncol = length(column_names), nrow = length(wells_list)))
  template_df
  colnames(template_df) <- column_names # add colnames
  template_df <- dplyr::relocate(template_df, .data$well, .after = tidyselect::last_col()) # move 'well' to RHS
  template_df$well <- wells_list # add wells

  # save df
  # if(outfolder != "."){
  #   # make folder if it doesn't exist already
  #   ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)
  #   message("Saving minimal template in folder: ", outfolder, " ...")
  # } else {
  #   message("Saving minimal template in current folder...")
  # }
  utils::write.csv(template_df, file.path(outfolder, filename), row.names = FALSE)
}
