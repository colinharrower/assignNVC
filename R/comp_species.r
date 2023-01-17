#' Compare species names between the sample dataset and the comparison
#' pusedo-quadrat dataset
#'
#' @description Compare the species names used in the sample dataset against
#'   those used in the pseudo-quadrat dataset and highlight any species names
#'   from the sample not present in the pseudo-quadrat data it is to be compared
#'   against. Mismatched names may be the result of species in the focal data
#'   that are not included in or relevant to the NVC psueod-quadrat, in which
#'   case they can be safely ignored. However it is also possible that they may
#'   be due to differences in nomenclature used in the two datasets and that the
#'   species could be in both but under different names. Ideally the datasets
#'   should both use a  standardised nomenclature, however this may not often be
#'   possible.
#' @param samp_df A `data.frame` containing the sample data that you want to
#'   assign to NVC via comparison with psuedo-quadrat data
#' @param comp_df  A `data.frame` containing the comparison data samples you
#'   want to match your focal samples in `samp_df` to. The default value of
#'   `NULL` means that the comparison will use the NVC pseudo-quadrat data
#'   included with the package
#' @param spp_col A character vector specifying the name of the column in
#'   `samp_df` and `comp_df`that contain the species names or codes that you
#'   want to match on. The default is for the column to be named `species`.
#' @param verbose A logical variable determining whether the function should
#'   output text to the display
#'
#' @return A `character` vector containing the names of species from the sample
#'   data that were not matched to names in the psuedo-quadrant data
#' @export
#'
#' @details Currently the NVC pseudo-quadrat data uses relatively old
#'   nomenclature so there are likely to be issues with mismatched names when
#'   compared against datasets using relatively recent nomenclature. In future
#'   we plant to update the nomenclature used in the nvc psudo-quadrats to a
#'   more modern nomenclature and would like to ideally include or link to
#'   packages that can help with the standardisation of taxnomomies between the
#'   datasets but for now you should try to ensure that where possible your
#'   nomenclature matches the NVC nomenclature as fully as possible.
#'
#' @examples
#'
#' #' # Create a test dataset
#' test_samp = data.frame(ID = 1, species = c("Alchemilla alpina","Anthoxanthum odoratum","Blechnum spicant","Deschampsia flexuosa","Galium saxatile","Luzula sylvatica","Melampyrum pratense","Nardus stricta","Potentilla erecta","Rumex acetosella","Vaccinium myrtillus","Dicranum majus","Hylocomium splendens","Isothecium myosuroides","Plagiothecium undulatum","Pleurozium schreberi","Polytrichum commune","Racomitrium lanuginosum","Rhytidiadelphus loreus","Calypogeia fissa","Unidentified plant","Cladonia bellidiflora","Cladonia coccifera","Cladonia uncialis"))
#'
#' # Now run names comparison function on this test data
#' unmatched_spp = compare_spp_names(test_samp)
#'
compare_spp_names = function(samp_df, comp_df = NULL,spp_col = "species",verbose = TRUE){
  # if comp_df null then matching against pseudo-quadrats in built-in psuedo-quadrat dataset
  if(is.null(comp_df)){
    comp_df = nvc_pquads
  }
  # Get all unique species names from samp_df
    samp_names = unique(samp_df[,spp_col])
  # Get all unique species names from comp_df
    ps_names = unique(comp_df[,spp_col])
  # Find all names not in psuedo-quadrat/comparison df
    miss_names = samp_names[which(!samp_names %in% ps_names)]
    n_miss = length(miss_names)
    if(verbose){
        if(n_miss > 0){
          cat("warning: ",n_miss," name(s) from samp_df not found in psuedo-quadrat/comparison data\n",sep="")
        } else{
          cat("All names in samp_df present in psuedo-quadrat/comparison data\n")
        }
    }
    if(length(miss_names) == 0){
        miss_names = NULL
    }
  # Return unmatched names
    invisible(miss_names)
}
