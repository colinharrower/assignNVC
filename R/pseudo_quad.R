#' Generate pseudo-quadrat data from a data.frame of species with occurrence frequencies
#'
#' @description Produce data for pseudo-quadrat(s) from a species via random sampling using occurrence frequencies/probabilities of occurrence for each species to produce random samples of species that could have been obtained
#'
#' @param sp_freq A 2 column \code{data.frame} containing the species names, or codes, along with the frequency probability of occurrence
#' @param nq An integer specifying the number of pseudo-quadrats to create. Where `sp_freq` contains data for multiple NVCs/groupings, indicated by values in `nvc_col` a vector of integer values, can be supplied if you want different numbers of pseudo-quadrats for the different NVC/groupings.  In this case the values in the `nq` vector must match the order in which the NVCs/groupings in `nvc_col` would be sorted, i.e. the 1st nq value is the number of quadrats for the 1st NVC/group in `nvc_col`, the 2nd nq value the number of quadrats for the 2nd NVC, etc.
#' @param spp_col A character string giving the name of the column in sp_freq that contains the species names/codes. Default is "species"
#' @param freq_col A character string giving the name of the column in sp_freq that contains the probability of occurrence for the give species. Default is "freq"
#' @param nvc_col A character string giving the name of the column in `sp_freq` that contains the NVC codes (or other groupings). When specified
#' it assumes that `sp_freq` contains separate species lists and associated frequencies for a number of different NVC communities or other groupings with that column specifying
#' the NVC/group to which each row belongs . Where this is not null a separate set of `nq` number of pseudo-quadrats will be generated
#' for each NVC/group in `nvc_col` of `sp_freq`
#'
#'
#' @return A 2 column \code{data.frame} with the first column giving the pseudo-quadrat number and 2nd the species names/codes of species present in the pseudo-quadrat
#' @export
#'
#' @examples
#' # Setup a data.frame with a species list and occurrency frequencies
#'
#'   sp_freq <- data.frame(species = paste("Species",1:10),freq = c(0.9,0.8,0.4,0.4,0.4,0.3,0.1,0.1,0.1,0.1))
#'
#' # Now produce 5 pseudo-quadrat samples based on this composition data.frame
#'
#'   test_ps <- pseudo_quad(sp_freq,5)
#'
#'
#' # Now an example where the data.frame contains species lists & frequencies for multiple NVCs each of which you want to generate pseudo-quadrats for
#'  # Load NVC_communites dataset
#'    data(NVC_communities)
#'  # Subset to the first 5 NVCs
#'    sp_freq = subset(NVC_communities,NVC %in% unique(NVC_communities$NVC)[1:5])
#'  # Create pseudo-quarats for these NVCs (note spp_col in NVC communities is Species not species as expected by default)
#'    test_ps = pseudo_quad(sp_freq,nq = c(10,5,5,5,5),nvc_col = "NVC",spp_col = "Species")
#'
pseudo_quad = function(sp_freq, nq, spp_col = "species", freq_col = "freq", nvc_col = NULL){
  # If nvc_col is not null then sp_freq contains data for multiple different sample types and nq quadrats from the species for that type should be generated
  if(!is.null(nvc_col)){
    n_nvc = length(unique(sp_freq[,nvc_col]))
    # If more than one grouping but nq only one value then repeat nq to create vector of same length as number of nvc
    if(n_nvc > 1 & length(nq) == 1){
      nq = rep(nq,n_nvc)
    }
  } else {
    n_nvc = 1
  }
  # Check that length of nq is the same as number of NVC/groupings
  stopifnot(length(nq) == n_nvc)
  # If more than one sample type found then use by to run pseudo_quad function independently for each sample type
  if(n_nvc > 1){
      # Setup temp object to hold data
        #ps_lt = by(sp_freq,INDICES = sp_freq[,nvc_col], FUN = pseudo_quad, nq = nq, spp_col = spp_col, freq_col = freq_col, nvc_col = nvc_col, simplify = FALSE)
        temp_sf = by(sp_freq, factor(sp_freq[,nvc_col],levels = unique(sp_freq[,nvc_col])), list)
        ps_lt = mapply(FUN = pseudo_quad, sp_freq = temp_sf,nq = nq, MoreArgs = list(spp_col = spp_col, freq_col = freq_col, nvc_col = nvc_col), SIMPLIFY = FALSE)
      # Now collapse the by list
        temp_ps = do.call("rbind",ps_lt)
        row.names(temp_ps) = NULL
  } else {
    # Determine number of species in sp_freq
      n_spp = nrow(sp_freq)
    # Build a data.frame to hold data
    if(!is.null(nvc_col)){
      temp_ps = data.frame(samp_id = sp_freq[1,nvc_col],ps_no = rep(1:nq,each = n_spp),sp_freq[rep(1:n_spp, nq),c(spp_col,freq_col)], occ = NA, row.names = NULL)
      # Now update samp_id column to name given in nvc_col
      names(temp_ps)[1] = nvc_col
    } else {
      temp_ps = data.frame(ps_no = rep(1:nq,each = n_spp),sp_freq[rep(1:n_spp, nq),c(spp_col,freq_col)], occ = NA, row.names = NULL)
    }
    # Now draw from bionomial for each row using the value in freq as the probability
      temp_ps[,"occ"] = stats::rbinom(nrow(temp_ps),1,prob = temp_ps[,freq_col])
    # Now subset to only data where species is included in the pseudo-quadrat (occ == 1)
      temp_ps = subset(temp_ps, occ == 1)
    # Reset row.names
      row.names(temp_ps) = NULL
  }
  # Now return data (only returning ps_no and species columns)
  ret_cols = c(nvc_col,"ps_no",spp_col) # if nvc_col is null this will produce a vector of column names without it
    return(temp_ps[,ret_cols])
}
