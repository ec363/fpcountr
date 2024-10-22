#' @import ghql
#' @import jsonlite

# global GraphQL client
client <- ghql::GraphqlClient$new(url = "https://www.fpbase.org/graphql/")
# global protein ID mapping cache
protein_id_map_cache <- NULL

#' Get detailed info for a specific protein at FPbase
#'
#' @param name_or_slug_or_id A string representing the case-insensitive name,
#'   slug, or FPbase accession ID of the protein to retrieve.
#' @return A nested list containing details of the protein.
#'
#' @examples
#' \dontrun{
#' get_protein_detail("mcherry")
#' }
#'
#' @export
get_protein_detail <- function(name_or_slug_or_id) {
  # Get the protein ID from the mapping
  protein_id <- protein_id_map()[[name_or_slug_or_id]]
  if (is.null(protein_id)) {
    stop("Not recognized as a protein name, slug, or id: ", name_or_slug_or_id)
  }
  return(get_protein_detail_by_id(protein_id))
}

#' Get detailed info for a specific dye at FPbase
#'
#' @param name_or_slug_or_id A string representing the case-insensitive name
#'   of the dye to retrieve.
#' @return A list containing details of the dye.
#'
#' @examples
#' \dontrun{
#' get_dye_detail("Alexa Fluor 488")
#' }
#'
#' @export
get_dye_detail <- function(name_or_slug) {
  qry <- ghql::Query$new()
  qry$query("dyeDetail", "query getDyeDetail($dyeName: String!)
    {
      dye(name:$dyeName) {
        id
        name
        slug
        exMax
        emMax
        extCoeff
        qy
        brightness
        pka
        lifetime
        manufacturer
        url
        spectra {
          subtype
          data
        }
      }
    }
  ")
  # Execute the GraphQL query
  tryCatch(
    {
      response <- client$exec(
        qry$queries$dyeDetail, list(dyeName = name_or_slug)
      )
    },
    error = function(e) {
      stop("Error in call to FPbase API: ", e$message)
    }
  )
  # Parse the result as JSON
  d <- jsonlite::fromJSON(response)
  # if the data has a keys "errors", print the first error message
  if ("errors" %in% names(d)) {
    stop("Error in call to FPbase API: ", d$errors[[1]])
  }
  return(d$data$dye)
}


#' Get detailed info for a specific protein at FPbase
#'
#' @param protein_id A string representing the ID of the protein to retrieve.
#' @return A list containing details of the protein.
#'
#' @examples
#' \dontrun{
#' get_protein_detail_by_id("ZERB6")
#' }
get_protein_detail_by_id <- function(protein_id) {
  # Initialize the GraphQL client

  # Define the GraphQL query
  # Use https://www.fpbase.org/graphql/ to play with queries and paste
  # your desired query below.  The query here shows the structure of the
  # data that will be returned.
  qry <- ghql::Query$new()
  qry$query("proteinDetail", "query getProteinDetail($protId: String!)
    {
      protein(id: $protId) {
        id
        name
        slug
        aliases
        seq
        seqValidated
        genbank
        agg
        switchType
        primaryReference{
          doi
          journal
          date
          title
        }
        states {
          slug
          exMax
          emMax
          extCoeff
          qy
          brightness
          pka
          lifetime
          maturation
          spectra {
            subtype
            data
          }
        }
      }
    }
  ")

  # Execute the GraphQL query
  tryCatch(
    {
      response <- client$exec(
        qry$queries$proteinDetail, list(protId = protein_id)
      )
    },
    error = function(e) {
      stop("Error in call to FPbase API: ", e$message)
    }
  )

  # Parse the result as JSON
  d <- jsonlite::fromJSON(response)
  # if the data has a keys "errors", print the first error message
  if ("errors" %in% names(d)) {
    stop("Error in call to FPbase API: ", d$errors[[1]])
  }

  return(d$data$protein)
}

#' Memoized function to fetch the protein ID mapping
#'
#' Facilitates easy lookup of protein IDs from protein names and slugs.
#' Names and slugs are case-insensitive.
#'
#' @return A list containing the mapping from protein slug/name to ID.
#'
#' @export
protein_id_map <- function() {
  # If the cache is empty, fetch the protein ID mapping
  if (is.null(protein_id_map_cache)) {
    protein_id_map_cache <<- fetch_protein_id_map()
  }
  return(protein_id_map_cache)
}

#' Fetch Protein ID Mapping from FPbase
#'
#' Facilitates easy lookup of protein IDs from protein names and slugs. There is
#' a memoization version of this function, `protein_id_map` which should be used
#' most of the time.
#'
#' @return A list containing the mapping from protein slug/name to ID.
#'
#' @examples
#' \dontrun{
#' map <- fetch_protein_id_map()
#' map[["mClover"]]
#' }
#'
fetch_protein_id_map <- function() {
  qry <- ghql::Query$new()
  qry$query("proteinList", "{ proteins { id name slug } }")

  # Execute the query and parse the JSON response
  response <- client$exec(qry$queries$proteinList)
  all_proteins <- jsonlite::fromJSON(response)$data$proteins

  mapping <- list()

  # Iterate over each row of the proteins data frame
  for (i in seq_len(nrow(all_proteins))) {
    # Map the protein name and slug to the protein ID
    name <- all_proteins[i, "name"]
    id <- all_proteins[i, "id"]
    slug <- all_proteins[i, "slug"]
    mapping[[slug]] <- id
    mapping[[name]] <- id
    # also map lower case names to the id
    mapping[[tolower(name)]] <- id
  }

  return(mapping)
}
