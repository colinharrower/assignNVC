#' Generate psuedo-quadrat data from a data.frame of species with occurrence frequencies
#'
#' @description Produce data for psuedo-quadrat(s) from a species via random sampling using occurrence frequencies/probabilities of occurrence for each species to produce random samples of species that could have been obtained
#'
#' @param nq The number of psuedo-quadrats to create
#' @param sp_freq A 2 column \code{data.frame} containing the species names, or codes, along with the frequency probability of occurrence
#' @param spp_col A character string giving the name of the column in sp_freq that contains the species names/codes. Default is "species"
#' @param freq_col A character string giving the name of the column in sp_freq that contains the probility of occurrence for the give species. Default is "freq"
#'
#' @return A 2 column \code{data.frame} with the first column giving the pseudo-quadrat number and 2nd the species names/codes of species present in the psuedo-quadrat
#' @export
#'
#' @examples
#' # Setup a data.frame with a species list and occurrency frequencies
#'
#'   sp_freq <- data.frame(species = paste("Species",1:10),freq = c(0.9,0.8,0.4,0.4,0.4,0.3,0.1,0.1,0.1,0.1))
#'
#' # Now produce 5 pseudo-quadrat samples based on this composition data.frame
#'
#'   test_ps <- psuedo_quad(5,sp_freq)
#'
#'
psuedo_quad = function(nq, sp_freq, spp_col = "species", freq_col = "freq"){
  # Determine number of species in sp_freq
  n_spp = nrow(sp_freq)
  # Build a data.frame to hold data
  temp_ps = data.frame(ps_no = rep(1:nq,each = n_spp),sp_freq[rep(1:n_spp, nq),c(spp_col,freq_col)], occ = NA, row.names = NULL)
  # Now draw from bionomial for each row using the value in freq as the probability
  temp_ps[,"occ"] = stats::rbinom(nrow(temp_ps),1,prob = temp_ps[,freq_col])
  # Now subset to only data where species is included in the psuedo-quadrat (occ == 1)
  temp_ps = subset(temp_ps, occ == 1)
  # Reset row.names
  row.names(temp_ps) = NULL
  # Now return data (only returning ps_no and species columns)
  return(temp_ps[,c("ps_no",spp_col)])
}
