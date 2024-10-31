#' View k-factors
#'
#' Displays the data tables used for k-factor calculations in `get_kfactors()`, along with notes.
#'
#' @export
#' @examples
#' view_kfactors()
view_kfactors <- function(){
  message("\nData from https://static.thermoscientific.com/images/D20827~.pdf :")

  # Get buffers data
  kfactors_buffers_data <- fpcountr::kfactors_buffers_data
  print(kfactors_buffers_data)

  message("")

  # Get temperature data
  kfactors_temperature_data <- fpcountr::kfactors_temperature_data
  print(kfactors_temperature_data)

  message("\nThe k-factor of water (A975-A900 at 1cm) at ~25oC is 0.172.")
  message("Increasing the temperature has a small effect: 30oC = 3% increase to 0.176.")
  message("The volume of water is minimally effected by temperature, so path lengths calculated at 25oC should apply to other temperatures.")

  message("\nMost aqueous buffers with low salt will be ~ 0.170 and can therefore be approximated to the k-factor of pure water.")
  message("eg.")
  message("153mM (0.9%) NaCl = 0.168, TE = 0.169, 50mM TBS = 0.166")
  message("1M Tris-HCl = 0.157, 2M H2SO4 = 0.154")

  message("\nOther solvents, like DMSO, EtOH and MeOH have a bigger effect.")
  message("eg.")
  message("10% DMSO = 0.148")
  message("20% EtOH = 0.092, 20% MeOH = 0.100")
  message("100% EtOH and MeOH have negative k-factors, and k-factors with these solvents are no longer reliable methods for calculating path lengths.")
}
