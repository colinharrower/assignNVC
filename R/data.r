#' Psuedo-quadrat data for GB NVC Comunities
#'
#' @description  A dataset containing pseudo-quadrat data for NVC
#' communities generated by random sampling from NVC species
#' list. The psuedo-quadrat data comes from a study by Tipping \emph{et al.} 2013.
#'
#' @details When generating psudeo-quadrants the  probabilities
#' of species being in a sample were set to match the
#' constancy scores of the species in that given NVC
#' community/sub-community. The probabilities associated
#' with each constancy score were 0.1, 0.3, 0.5, 0.7 and 0.9 for
#' constancy scores 1, 2, 3, 4 and 5 respectively. The number
#' of psuedo-quadrats generated for each NVC community or
#' sub-community was set to match the number of samples originally
#' used to define the NVC and produce the contingency tables.
#'
#' The citation for the paper is:
#' Tipping, E., Henrys, P. A., Maskell, L. C., & Smart, S. M. (2013). Nitrogen deposition effects on plant species diversity; Threshold loads from field data. \emph{Environmental Pollution}, 179, 218-223.
#'
#'
#'
#' @format A data frame with 709695 rows and 6 variables:
#' \describe{
#'   \item{name}{species scientific name}
#'   \item{spp}{species code}
#'   \item{Pid2}{psuedo-quadrat ID}
#'   \item{Pid3}{ID created concatenating NVC and psuedo-quadrat ID `Pid2`}
#'   \item{Pid3}{psued-quadrat number}
#' }
#'
#' @seealso [nvc_pquads()] for updated version
#'
#' @source \url{http://dx.doi.org/10.1016/j.envpol.2013.04.008}
"ps_quad"

#' Psuedo-quadrat data for GB NVC Comunities (updated)
#'
#' @description  A dataset containing pseudo-quadrat data for NVC
#' communities generated by random sampling from NVC species
#' list. The psuedo-quadrat data is based on that from a study by Tipping \emph{et al.} 2013.
#' but has been updated to include additional NVCs
#'
#' @details When generating psudeo-quadrants the  probabilities
#' of species being in a sample were set to match the
#' constancy scores of the species in that given NVC
#' community/sub-community. The probabilities associated
#' with each constancy score were 0.1, 0.3, 0.5, 0.7 and 0.9 for
#' constancy scores 1, 2, 3, 4 and 5 respectively. The number
#' of psuedo-quadrats generated for each NVC community or
#' sub-community was set to match the number of samples originally
#' used to define the NVC and produce the contingency tables.
#'
#' The citation for the Tipping \emph{et al.} 2013 paper is:
#' Tipping, E., Henrys, P. A., Maskell, L. C., & Smart, S. M. (2013). Nitrogen deposition effects on plant species diversity; Threshold loads from field data. \emph{Environmental Pollution}, 179, 218-223.
#'
#'
#'
#' @format A data frame with 974601 rows and 6 variables:
#' \describe{
#'   \item{name}{species scientific name}
#'   \item{spp}{species code}
#'   \item{Pid2}{psuedo-quadrat ID}
#'   \item{Pid3}{ID created concatenating NVC and psuedo-quadrat ID `Pid2`}
#'   \item{Pid3}{psued-quadrat number}
#' }
#' @source \url{http://dx.doi.org/10.1016/j.envpol.2013.04.008}
"nvc_pquads"


#' NVC community floristic tables
#'
#' @description  A dataset of the NVC community floristic table information.
#' For each NVC community/sub community this lists the species present and their
#' respective constancy and frequency values.
#'
#' @details The data in this dataset can be used to generate new psuedo-quadrats
#' for each NVC community/sub community.
#'
#'
#' @format A data frame with 48194 rows and 5 variables:
#' \describe{
#'   \item{NVC}{The NVC community or sub-community code}
#'   \item{BRC}{species ID code as used historically by the BRC}
#'   \item{Species}{The species/taxa scientific name}
#'   \item{Constancy}{The constancy value for the species specifying how commonly the species are likely to be in
#'   samples from that NVC. Ranges from 1 to 5 where 1 is rare and 5 is constant}
#'   \item{freq}{The probability/frequency of occurrence, these are related to the constancy where 1 = 0.1, 2 = 0.3, 3 = 0.5, 4 = 0.7, 5 = 0.9}
#' }
#'
"NVC_communities"

#' NVC Sample sizes
#'
#' @description  A dataset listing the sizes of the samples used to define that NVC.
#' The psuedo-quadrat included with the package has been generated so that the number of
#' psuedo-quadrats for each NVC matches the original data sample sizes.
#'
#' @details Using the data in dataset a new set(s) of psuedo-quadrats that match the original sample
#' sizes can be generated
#'
#'
#' @format A data frame with 879 rows and 8 variables:
#' \describe{
#'   \item{NVC}{The NVC community or sub-community code}
#'   \item{LEVEL}{The level of the NVC given, 1 for community, 2 for sub-community}
#'   \item{NVC_COMMUNITY}{The NVC community that the NVC belongs to}
#'   \item{VEG_TYPE}{The vegetation type letter code for the NVC, i.e. the first letter(s)}
#'   \item{NVC_NUMBER}{The number of that NVC within the current `VEG_TYPE`, i.e. the digit(s) after the VEG_TYPE in the NVC code}
#'   \item{NVC_SUB_CODE}{The community code for that NVC, if applicable}
#'   \item{VERSION}{The version of that NVC as some NVC have been updated, with new sub-communities being added meaning that for some communities there is an original version and then an updated v2 version that include additional species/samples from the added subcategories}
#'   \item{N_SAMPLES}{The number of samples original used in the data used to define the NVC community/sub-community}
#' }
#'
"nvc_samp_sizes"
