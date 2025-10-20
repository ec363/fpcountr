#' Amino acid masses
#'
#' A dataset containing the masses of all amino acids.
#'
#' @format A data frame with 22 rows and 5 variables: \describe{
#'   \item{amino.acid}{amino acid name}
#'   \item{X1.letter.code}{amino acid in 1-letter code format}
#'   \item{X3.letter.code}{amino acid in 3-letter code format}
#'   \item{monoisotopic}{monoisotopic mass, Da}
#'   \item{average}{average mass, Da}
#'   }
#' @source \url{https://web.expasy.org/findmod/findmod_masses.html#AA}
"aa_mass_data"

#' k-factors of different buffers
#'
#' A dataset containing the k-factors (k-factor = A975-A900 at 1cm pathlength)
#' of a range of aqueous buffers with salts and solvents, obtained from figure 5
#' of a Thermo Scientific application note on 'Microplate Based Pathlength
#' Correction Method for Photometric DNA Quantification Assay' (Lampinen,
#' Raitio, Perälä, Oranen and Harinen, 2012, ANMR_MUGO_0412). The dataset was
#' extracted using WebPlotDigitizer (`https://automeris.io/WebPlotDigitizer/`),
#' tidied to allow easier handling for programmatic purposes, and a column was
#' added to calculate fold change in k-factors for all buffers compared to
#' water. Source: `https://static.thermoscientific.com/images/D20827~.pdf`
#'
#' @format A data frame with 10 rows and 6 variables: \describe{
#'   \item{buffer}{name of the buffer} \item{concentration}{concentration of the
#'   buffer, units are variable} \item{units}{units of the concentration column,
#'   M is molar, pc is percentage i.e. w/v (%)} \item{description}{string
#'   containing buffer-concentration-units information all in one}
#'   \item{kfactor}{k-factor value (A975-A900 at 1cm pathlength)}
#'   \item{fold_change}{fold change in k-factor from water} }
"kfactors_buffers_data"

#' k-factors of water at different temperatures
#'
#' A dataset containing the k-factors (k-factor = A975-A900 at 1cm pathlength)
#' of water at different temperatures, obtained from figure 6 of a Thermo
#' Scientific application note on 'Microplate Based Pathlength Correction Method
#' for Photometric DNA Quantification Assay' (Lampinen, Raitio, Perälä, Oranen
#' and Harinen, 2012, ANMR_MUGO_0412). The dataset was extracted using
#' WebPlotDigitizer (`https://automeris.io/WebPlotDigitizer/`), and a column was
#' added to calculate fold change in k-factors for temperatures compared to
#' 25oC. Source: `https://static.thermoscientific.com/images/D20827~.pdf`
#'
#' @format A data frame with 10 rows and 6 variables: \describe{
#'   \item{temperature}{temperature in oC} \item{kfactor}{k-factor value
#'   (A975-A900 at 1cm pathlength)} \item{fold_change}{fold change in k-factor
#'   from 25oC} }
"kfactors_temperature_data"

#' pathlength vs volume data
#'
#' A dataset containing experimentally derived pathlengths for microplate wells
#' containing water at different volumes. Pathlengths were derived by using the
#' equation `pathlength (cm) = kfactor_well / kfactor_1cm`. `Kfactor_well` is
#' defined as A975-A900 value for each well (there were 4 replicates per volume
#' per experiment). `Kfactor_1cm` values were derived using the `get_kfactor()`
#' function for water at 26oC.
#'
#' @format A data frame with 10 rows and 6 variables: \describe{
#'   \item{volume}{volume in ul} \item{pathlength}{pathlength in cm}
#'   \item{expt}{experiment number denoting experimental repeat} }
"pathlength_water_data"

#' fluorescent protein properties from fpbase
#'
#' A dataset containing some of the properties of all the 'basic' (single
#' fluorescent state) fluorescent proteins in the fpbase database (downloaded
#' 2022.10.04).
#'
#' @format A data frame with 593 rows and 16 variables: \describe{
#'   \item{url}{FPbase url} \item{name}{name} \item{stokes}{Stokes shift (nm)}
#'   \item{slug}{slug - lower case name used in url} \item{ipg_id}{Identical
#'   Protein Group ID at Pubmed} \item{agg}{oligomerisation state}
#'   \item{ex_max}{wavelength of excitation peak (nm)} \item{em_max}{wavelength
#'   of emission peak (nm)} \item{ext_coeff}{excitation coefficient (EC,
#'   M-1cm-1)} \item{qy}{quantum yield} \item{pka}{pKa}
#'   \item{brightness}{brightness} \item{bleach}{rate of
#'   photobleaching/photostability (s)} \item{maturation}{maturation half time
#'   (min)} \item{lifetime}{fluorescence lifetime (ns)} \item{cofactor}{cofactor
#'   required for fluorescence (if any)}}
#' @source \url{https://www.fpbase.org/api/proteins/basic/?format=json}
"fpbase_data"

#' protein sequence data
#'
#' A table containing the protein sequences used in the FPCountR paper.
#'
#' @format A data frame with 4 rows and 3 variables: \describe{
#'   \item{name}{protein name} \item{slug}{protein slug - lower case name used
#'   in url of fpbase entry} \item{sequence}{primary sequence in 1-letter code}
#'   }
"protein_seq_data"

#' cell quench data
#'
#' A table containing experimental data for informing how cells quench fluorescence of E. coli cells.
#'
#' @format A data frame with 34 rows and 12 variables: \describe{
#'  \item{channel_name}{name of filter set used}
#'  \item{channel_ex}{excitation filter}
#'  \item{channel_em}{emission filter}
#'  \item{media}{media}
#'  \item{calibrant}{calibrant}
#'  \item{protein}{protein}
#'  \item{cells_od600}{accurate OD600 cm-1 values}
#'  \item{cells_od700}{accurate OD700 cm-1 values}
#'  \item{volume}{volume in ul}
#'  \item{measure}{filter set plus gain}
#'  \item{gain}{gain}
#'  \item{autof_norm_value}{fluorescence values normalised to remove the contribution of cellular autofluorescence}
#'  \item{norm_fc_value}{fold change between the fluorescence sample and its relevant negative control without cells}
#'   }
"cell_quench_data"

