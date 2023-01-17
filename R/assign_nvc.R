#' Determine most similar NVCs for each sample by comparing samples against NVC psuedo-quadrat data
#'
#' @param samp_df A `data.frame` containing the sample data that you want to assign to NVC via comparison with psuedo-quadrat data
#' @param comp_df A `data.frame` containing the comparison data samples you want to match your focal samples in `samp_df` to. The default value of `NULL` means that the comparison will use the NVC pseudo-quadrat data included with the package
#' @param spp_col A character vector specifying the name of the column in `samp_df` and `comp_df`that contain the species names or codes that you want to match on. The default is for the column to be named `species`.
#' @param samp_id A charcter vector specifying the name of the column in `samp_df` that contains the IDs identifying the data for each sample in the data you want assigned to NVC. The default is for the column to be named `ID`.
#' @param comp_id A charcter vector specifying the name of the column in `comp_df` that contains the IDs identifying the data for each sample in the comparison data. The default is for the column to be named `ID`.
#' @param top_n An integer specifying the number most similar comparison samples to be returned, or `NULL` to return all matches. The default `top_n = 20` will return the 20 most similar samples from `comp_df` for each sample in `samp_df`.
#'
#' @return a `data.frame` detailing the most similar NVC psuedo-quadrats for each sample in the focal data `samp_df`
#' @export
#'
#' @examples
#'
#' # Create a test dataset
#' test_samp = data.frame(ID = 1, species = c("Alchemilla alpina","Anthoxanthum odoratum","Blechnum spicant","Deschampsia flexuosa","Galium saxatile","Luzula sylvatica","Melampyrum pratense","Nardus stricta","Potentilla erecta","Rumex acetosella","Vaccinium myrtillus","Dicranum majus","Hylocomium splendens","Isothecium myosuroides","Plagiothecium undulatum","Pleurozium schreberi","Polytrichum commune","Racomitrium lanuginosum","Rhytidiadelphus loreus","Calypogeia fissa","Cladonia arbuscula","Cladonia bellidiflora","Cladonia coccifera","Cladonia uncialis"))
#'
#' # Now run assign_nvc on this data to get the top 10 most similar
#' test_jac = assign_nvc(test_samp, top_n=10)
#'
assign_nvc = function(samp_df, comp_df = NULL, spp_col = "species", samp_id = "ID", comp_id = "Pid3", top_n = 20){
  # if comp_df null then matching against pseudo-quadrats in built-in psuedo-quadrat dataset
  if(is.null(comp_df)){
    comp_df = nvc_pquads
  }
  # Determine number of samples in samp_df
    in_samps = unique(samp_df[,samp_id])
    n_samp = length(in_samps)
  # Produce list of unique names from comp_df (unique for each comparison sample not unique across all samples)
    nvc_pq_spp = tapply(comp_df[,spp_col], comp_df[,comp_id],unique)
  # Create object to hold output data
      out_jac = vector("list",n_samp)
  # Loop through samples and determine match with all comparions samples in comp_df
  for(i_s in 1:n_samp){
    # Get data from current sample
      cur_samp = subset(samp_df, samp_df[,samp_id] == in_samps[i_s])
    # Get vector of unique species in sample
      cur_spp = unique(cur_samp[,spp_col])
    # Determine similarity
      out_jac[[i_s]] = jaccard_sim(cur_spp, nvc_pq_spp, focal_id = in_samps[i_s], top_n = top_n)
  }

  # Collpase list to data.frame
    out_obj = do.call("rbind",out_jac)
  # Return object
    return(out_obj)
}



#' Determine the average similarities between a set of focal sample(s) and a set of NVC psuedo-quadrat sample(s)
#'
#' @param samp_df A `data.frame` containing the sample data that you want to assign to NVC via comparison with psuedo-quadrat data
#' @param comp_df A `data.frame` containing the comparison data samples you want to match your focal samples in `samp_df` to. The default value of `NULL` means that the comparison will use the NVC pseudo-quadrat data included with the package/
#' @param spp_col A character vector specifying the name of the column in `samp_df` and `comp_df`that contain the species names or codes that you want to match on. The default is for the column to be named `species`.
#' @param samp_id A charcter vector specifying the name of the column in `samp_df` that contains the IDs identifying the data for each sample in the data you want assigned to NVC. The default is for the column to be named `ID`.
#' @param comp_id A charcter vector specifying the name of the column in `comp_df` that contains the IDs identifying the data for each sample in the comparison data. The default is for the column to be named `ID`.
#' @param exc_zero_match A logical value determining whether NVC that have an average match of zero should be excluded from the returned results
#'
#' @return A `data.frame` returning the mean and standard deviation of jaccard similarities between each focal sample and the pseudo-quadrat samples for each NVC
#' @export
#'
#' @examples
#'
#' #' # Create a test dataset
#' test_samp = data.frame(ID = 1, species = c("Alchemilla alpina","Anthoxanthum odoratum","Blechnum spicant","Deschampsia flexuosa","Galium saxatile","Luzula sylvatica","Melampyrum pratense","Nardus stricta","Potentilla erecta","Rumex acetosella","Vaccinium myrtillus","Dicranum majus","Hylocomium splendens","Isothecium myosuroides","Plagiothecium undulatum","Pleurozium schreberi","Polytrichum commune","Racomitrium lanuginosum","Rhytidiadelphus loreus","Calypogeia fissa","Cladonia arbuscula","Cladonia bellidiflora","Cladonia coccifera","Cladonia uncialis"))
#'
#' # Now run assign_nvc on this data to get the top 10 most similar
#' test_avg = nvc_average_sim(test_samp, exc_zero_match = TRUE)
#'
nvc_average_sim = function(samp_df, comp_df = NULL, spp_col = "species", samp_id = "ID", comp_id = "Pid3",exc_zero_match = TRUE){
  # if comp_df null then matching against pseudo-quadrats in built-in psuedo-quadrat dataset
  if(is.null(comp_df)){
    comp_df = nvc_pquads
  }
  # Determine number of samples in samp_df
  in_samps = unique(samp_df[,samp_id])
  n_samp = length(in_samps)
  # Produce list of unique names from comp_df (unique for each comparison sample not unique across all samples)
  nvc_pq_spp = tapply(comp_df$name, comp_df[,comp_id],unique)
  # Create object to hold output data
  out_avg = vector("list",n_samp)
  # Loop through samples and determine match with all comparions samples in comp_df
  for(i_s in 1:n_samp){
    # Get data from current sample
      cur_samp = subset(samp_df, samp_df[,samp_id] == in_samps[i_s])
    # Get vector of unique species in sample
      cur_spp = unique(cur_samp[,spp_col])
    # Determine similarity
      cur_jac = jaccard_sim(cur_spp, nvc_pq_spp, focal_id = in_samps[i_s], top_n = NULL)
    # Determine average and SD of JAC_SIM for each NVC
      cur_avg = tapply(cur_jac$JAC_SIM, cur_jac$NVC, mean,na.rm = TRUE)
      cur_sd = tapply(cur_jac$JAC_SIM, cur_jac$NVC, sd,na.rm = TRUE)
      # Check the names are the same
        stopifnot(names(cur_avg) == names(cur_sd))
      # Build results data.frame
        cur_res = data.frame(FOCAL_ID = in_samps[i_s], NVC = names(cur_avg), MEAN_SIM = cur_avg, SD = cur_sd)
        # If exc_zero_match is true then remove any that have a mean sim of zero
        if(exc_zero_match){
          cur_res = subset(cur_res, MEAN_SIM != 0)
        }
        # Reorder so that most similar are at the top
          cur_res = cur_res[order(cur_res$MEAN_SIM, cur_res$SD, decreasing = TRUE),]
        # Reset row.names
          row.names(cur_res) = NULL
      out_avg[[i_s]] = cur_res
  }

  # Collpase list to data.frame
  out_obj = do.call("rbind",out_avg)
  # Return object
  return(out_obj)
}
