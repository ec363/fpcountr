#' Get properties of an FP from FPbase
#'
#' Gets properties of an FP from FPbase, including the maximal excitation
#' wavelength and extinction coefficient, and saves these as csv.
#'
#' @param slug Name of FP, which is used to find the FP-relevant data lines in
#'   the FPbase datasets. This argument is called `slug` because for FPbase
#'   retrievals, this needs to be an exact match to the slug specified by
#'   FPbase. (To find this, navigate to the FPbase entry for your FP and copy
#'   the part after `https://www.fpbase.org/protein/` (without the trailing
#'   `/`).
#' @param verbose logical. Should the function print messages to allow the user
#'   to check its progress? Defaults to TRUE.
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#' @param filename How to name the output csv file. Requires ".csv" at the end.
#'   Defaults to "fp_properties.csv".
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#'
#' @examples get_properties("mcherry")
get_properties <- function(slug, verbose = TRUE, outfolder = ".", filename = "fp_properties.csv"){

  # Location for saved outputs ---------------------------------------------------------------------------------------------------

  # make folder if it doesn't exist already
  ifelse(test = !dir.exists(file.path(outfolder)), yes = dir.create(file.path(outfolder)), no = FALSE)

  ## Get FP properties -----------------------------------------------------------------------------------------------------------

  # In case of FPbase url error, wrap FPbase call in tryCatch
  tryCatch(
    expr = {

      all_fp_df <- "empty"

      # 1) fpbase version ---
      # All FP properties
      all_fp_df <- jsonlite::fromJSON("https://www.fpbase.org/api/proteins/?format=json")

      # Get FP entry
      fp_entry <- all_fp_df %>%
        dplyr::filter(slug == {{slug}})
      fp_entry

      # FP properties
      fp_properties <- fp_entry[[11]][[1]]
      fp_properties

      if(verbose){
        message("FP properties retrieved from FPbase:")
        print(fp_properties)
      }

      # Warning mesg for proteins with multiple states
      if(nrow(fp_properties)>1){
        message("FP has multiple states, taking first:")
        fp_properties <- fp_properties[1,]
        print(fp_properties)
      }
      if(nrow(fp_properties)<1){
        return(message("Error: no entries matching slug in FPbase."))
      }
      # fpbase version ---


    },
    error = function(e){
      # Do this if an error is caught...
      message('Error in call to FPbase:')
      print(e)
    },
    warning = function(w){
      # Do this if an warning is caught...
      message('Warning in call to FPbase:')
      print(w)
    },
    finally = {
      # do this at the end before quitting the tryCatch structure...
      # empty
    }
  )

  # If FPbase download fails, use package data instead:
  # 2) pre-downloaded version ---
  if(all_fp_df == "empty"){

    # All FP properties - from pre-downloaded data
    all_fp_df <- fpcountr::fpbase_data
    all_fp_df

    # Get FP entry
    fp_entry <- all_fp_df %>%
      dplyr::filter(slug == {{slug}})
    fp_entry

    # FP properties
    # fp_properties <- fp_entry[[11]][[1]] # don't need this here
    fp_properties <- fp_entry
    fp_properties

    if(verbose){
      message("FP properties retrieved from FPbase data stored in package:")
      print(fp_properties)
    }

    # Warning mesg for proteins with multiple states
    if(nrow(fp_properties)>1){
      message("FP has multiple states, taking first:")
      fp_properties <- fp_properties[1,]
      print(fp_properties)
    }
    if(nrow(fp_properties)<1){
      return(message("Error: no entries matching slug in FPbase."))
    }

  }
  # pre-downloaded version ---

  # Save and return ------------------------------

  write.csv(fp_properties, file.path(outfolder, filename), row.names = FALSE)

  return(fp_properties)
}
