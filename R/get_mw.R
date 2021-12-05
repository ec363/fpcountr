#' Get molecular weight of a protein
#'
#' Get molecular weight of protein in Daltons from its sequence alone.
#'
#' @param protein character string of protein sequence using the 1-letter code.
#'
#' @export
#'
#' @importFrom dplyr %>%
get_mw <- function(protein){

  # Table of AA masses ----------------------------------

  ## how to access data that is within a package - simply name it
  # eg # fpcountR::aa_mass_data
  # eg # aa_mass_data
  # source: https://web.expasy.org/findmod/findmod_masses.html#AA

  # Get AAs ----------------------------------
  amino_acids <- strsplit(protein, split="")
  # class(amino_acids) # list
  amino_acids <- unlist(amino_acids)
  # class(amino_acids) # character array

  # Get MW ----------------------------------

  mw <- 0

  # For each AA...
  for(aa in aa_mass_data$X1.letter.code){

    new_aa <- aa

    num_aa <- length(which(amino_acids == new_aa)) # count AA
    mass_aa <- aa_mass_data %>% # get mass of AA
      dplyr::filter(X1.letter.code == new_aa) %>%
      dplyr::select(average) %>%
      as.numeric()
    added_mass <- num_aa*mass_aa # find total mass of all of those AAs

    mw <- mw + added_mass # add

  }

  # Add one molecule of H2O ----------------------------------

  # As per Expasy protocol for protein MW:
  # https://web.expasy.org/compute_pi/pi_tool-doc.html
  # "Protein Mw is calculated by the addition of average isotopic masses of amino
  # acids in the protein and the average isotopic mass of one water molecule.
  # Molecular weight values are given in Dalton (Da)."

  mass_water <- 18.01524
  mw <- mw + mass_water

  return(mw)
}
