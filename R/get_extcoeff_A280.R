#' Get a protein's molar extinction coefficient at A280
#'
#' The purpose of this function is to work out the theoretical molar extinction
#' coefficient of a protein at 280nm (EC280, M-1cm-1) using only the protein's
#' primary sequence, using the ProtParam method (Pace values). A full
#' explanation of the method can be found at
#' https://web.expasy.org/protparam/protparam-doc.html
#'
#' @param protein character string of protein sequence using 1-letter code
#' @param disulphides logical. Does protein have disulphides?
#' @param showWarnings logical. Should function show warnings?
#' @param showMessages logical. Should function show messages?
#' @param protein_name character string of protein name. Optional.
#' @param buffer character string of buffer. Optional.
#' @param mol_weight numerical value for molecular weight (g/mol). Optional. If
#'   specified the function gives extinction coefficients for 1% (10mg/ml) and
#'   0.1% (1mg/ml) solutions too.
#' @param save logical. Should function save csv file of output?
#' @param outfolder path to folder where output files should be saved. Defaults
#'   to current working directory.
#'
#' @export
#'
#' @importFrom dplyr %>%
get_extcoeff_A280 <- function(protein, # "AAAWYCAAA"
                              disulphides = FALSE, showWarnings = TRUE, showMessages = TRUE,
                              protein_name = "-", buffer = "-", mol_weight = NULL,
                              save = TRUE, outfolder = "."
                              ){

  # Get AAs ----------------------------------
  amino_acids <- strsplit(protein, split="")
  amino_acids <- unlist(amino_acids)

  # Count AAs ----------------------------------
  num_trp <- length(which(amino_acids =="W"))
  num_tyr <- length(which(amino_acids =="Y"))
  num_cys <- length(which(amino_acids =="C"))
  if(disulphides == FALSE){num_cys <- 0}

  # Get molar ext coeff ----------------------------------

  # ProtParam method, see function description above for details.
  extcoeff <- 5500*num_trp + 1490*num_tyr + 125*0.5*num_cys
  # NB. 125 is absorbance of one cystine which is the product of TWO cysteines

  # Get other ext coeff ----------------------------------

  if(!is.null(mol_weight)){
    # Get MW
    mw <- mol_weight

    # # E(1%) = E(10mg/ml)
    # extcoeff2 <- extcoeff/mw*10

    # E(0.1%) = E(1mg/ml)
    extcoeff3 <- extcoeff/mw
  } else {
    mw <- "-"
    # extcoeff2 <- "-"
    extcoeff3 <- "-"
  }

  # Populate EC table ----------------------------------

  EC280_table <- data.frame(
    protein_name = protein_name,
    buffer = buffer,
    disulphides = disulphides,

    molecular_weight_gmol1 = mw,

    E_molar = extcoeff,
    E_molar_units = "M-1cm-1",
    # E_1pc = extcoeff2,
    # E_1pc_units = "10(mgml)-1(cm)-1",
    E_0.1pc = extcoeff3,
    E_0.1pc_units = "(mgml)-1(cm)-1",

    protein_seq = protein
    )
  EC280_table

  # Save table ----------------------------------
  if(save){
    csvname <- paste0("extcoeff_A280_", protein_name, ".csv")
    utils::write.csv(x = EC280_table, file = file.path(outfolder, csvname), row.names = FALSE)
  }

  # Warnings ----------------------------------
  if(showWarnings){
    message("get_extcoeff_A280. Warnings:")
    message("Proteins with no Trp residues may have >10% error.")
    message("Not suitable for proteins with chromophores that absorb at 280nm.")
    message("")
  }
  if(showMessages){
    message("get_extcoeff_A280. Messages:")
    # message("The extinction coefficient here is a molar extinction coefficient with units M-1cm-1.")
    # message("Beer's law: A = E*C*L, where E = molar extinction coefficient (M-1cm-1), C = concentration (M), L = path length (cm).")
    # message("")
    message("E_molar is the molar extinction coefficient for a 1M solution, with units of (M)-1(cm)-1.")
    # message("E_1pc is the mass extinction coefficient for a 10mg/ml solution, with units of 10(mgml)-1(cm)-1.")
    message("E_0.1pc is the mass extinction coefficient for a 1mg/ml solution, with units of (mgml)-1(cm)-1.")
  }

  return(EC280_table)
}

