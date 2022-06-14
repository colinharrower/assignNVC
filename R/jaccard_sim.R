#' Calculate Jaccard similarity between a focal sample and one or more comparison sample(s)
#'
#' @description Calculate the Jaccard similarity values between the species recorded for a focal site or sample against those recorded for a series of comparison sites or samples (e.g. against the NVC pseudo-quadrats) and then return the \code{top_n} matches in \code{comp_spp_list} based on their similarity scores with the species in focal sample
#'
#' @param focal_spp A character vector containing species names for the focal sample that to be compared against against the species in \code{comp_spp_list}
#' @param comp_spp_list A \code{list} containing the species recorded in the comparison samples, wheere each element is a seperate comparison sample and contains a character vector with the species recorded for that comparison sample. The names of the list elements will be the names/identifiers for the comparison samples
#' @param focal_id An name or identifier for the focal sample. Default is \code{NULL}
#' @param top_n An integer determining how many of the top scoring comparison samples and their similarity scores will be returned in the final data.frame. The default is \code{50}
#'
#' @return A \code{data.frame} containing the comparison samples, similarity scores for the \code{top_n} matches from \code{comp_spp_list}
#' @export
#'
#' @examples
jaccard_sim = function(focal_spp, comp_spp_list, focal_id = NULL,top_n = 50){
  # Determine number of species in focal_spp
    n_foc = length(focal_spp)
  # Determine the number of species in each elemetn of comp_spp_list
    n_comp = sapply(comp_spp_list,length)
  # Determine the number in common between focal list and each pseudo-quadrat
    n_both = sapply(lapply(comp_spp_list,"%in%",focal_spp),sum)
  # Calculate the Jaccard similarity measure
    jac = n_both / (n_foc + n_comp - n_both)
  # Reorder the similarity measures in decreasing order
    jac = sort(jac, decreasing = TRUE)
  # Select the top nth (as specified by top_n arugment)
    sel_jac = head(jac,top_n)
  # Build output data.frame containing the Pid3, the Jaccard Similarity Measure and NVC for the top nth number of similarities
  if(is.null(focal_id)){
    ret_obj = data.frame(Pid3 = names(sel_jac),JAC_SIM = sel_jac, NVC = gsub("(.*)(P[[:digit:]]{1,})$","\\1",names(sel_jac)), row.names = NULL, stringsAsFactors = FALSE)
  } else {
    ret_obj = data.frame(FOCAL_ID = focal_id, Pid3 = names(sel_jac),JAC_SIM = sel_jac, NVC = gsub("(.*)(P[[:digit:]]{1,})$","\\1",names(sel_jac)), row.names = NULL, stringsAsFactors = FALSE)
  }
  # Return ret_obj
  return(ret_obj)
}
