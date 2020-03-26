#' Get variance for observations of an attribute
#'
#' @param attr_name character; the name of the attribute to compute variance for
#' @param attr_vec vector of characters; the vector of values for that attribute corresponding to the rows in embedding
#' @param embedding matrix; the embeddings, where columns correspond to PCs and rows correspond to cells (entries)
#' @param pcnum int; the PC to compute the variance for
#'
#' @return double; the scaled variance for that subset of observations with attr_name
#' 
#'
#' @examples
get_var_for_attr <- function(attr_name, attr_vec, embedding, pcnum){
  return(var(embedding[which(attr_vec == attr_name),pcnum]))
}

#' Compute scaled variance by attribute for a single embedding
#'
#' @param embedding matrix; the embeddings, where columns correspond to PCs and rows correspond to cells (entries)
#' @param attr_vec vector of characters; the vector of values for that attribute corresponding to the rows in embedding
#' @param pcnum int; the PC to compute the scaled variance for
#'
#' @return vector; scaled variances for each unique attribute type in attr_vec for this embedding
#' @export
#'
#' @examples
scaled_var_vec_by_attr <- function(embedding, attr_vec, pcnum = 1){
  overall_var <- var(embedding[,pcnum])
  uniq_attrs <- unique(attr_vec)
  var_by_attr <- sapply(uniq_attrs, get_var_for_attr, attr_vec, embedding, pcnum)
  return(var_by_attr/overall_var)
}

#' Compute scaled variance by attribute for a single embedding, multiple PCs
#'
#' @param embedding matrix; the embeddings, where columns correspond to PCs and rows correspond to cells (entries)
#' @param attr_vec vector of characters; the vector of values for that attribute corresponding to the rows in embedding
#' @param pcnums vector; the PCs to compute the variances for 
#'
#' @return data.frame; columns correspond to PCs and the rows correspond to each unique attribute type in attr_vec for this embedding, on that PC
#' @export
#'
#' @examples
scaled_var_df_by_attr <- function(embedding, attr_vec, pcnums = c(1:3)){
  overall_var <- apply(embedding[,pcnums], MARGIN = 2, FUN = var)
  uniq_attrs <- unique(attr_vec)
  raw_var_df <- data.frame(sapply(uniq_attrs, get_var_for_attr, attr_vec, embedding, pcnums[1]))
  for(i in 2:length(pcnums)){
    raw_var_df <- cbind(raw_var_df, sapply(uniq_attrs, get_var_for_attr, attr_vec, embedding, pcnums[i]))
  }
  scaled_var_df <- sweep(raw_var_df, MARGIN = 2,STATS = overall_var,FUN = '/')
  colnames(scaled_var_df) <- pcnums
  return(scaled_var_df)
}

#' Compute scaled variances by attribute to compare multiple embeddings on a single PC
#'
#' @param embeds_list list of matrices; each matrix is an embedding, where columns correspond to PCs and rows correspond to cells (entries). Ideally the embeddings are named appropriately, as names for the output will be assigned based on names of the list elements.
#' @param attr_vec attr_vec vector of characters; the vector of values for that attribute corresponding to the rows in embedding
#' @param pcnum int; the PC to compute the scaled variance for
#'
#' @return data.frame; columns correspond to different embeddings from the list and the rows correspond to each unique attribute type in attr_vec for that embedding
#' @export
#'
#' @examples
#' embed_list <- list('corralm' = reducedDim(sce,'corralm'), 'corral' = reducedDim(sce,'corral'))
#' scaled_var_comparison(embed_list, colData(sce)$Method, 2)
scaled_var_comparison <- function(embeds_list, attr_vec, pcnum = 1){
  embednames <- names(embeds_list)
  comp_df <- data.frame(scaled_var_vec_by_attr(embeds_list[[1]], attr_vec, pcnum))
  for(i in 2:length(embeds_list)){
    comp_df <- cbind(comp_df, scaled_var_vec_by_attr(embeds_list[[i]], attr_vec, pcnum))
  }
  colnames(comp_df) <- embednames
  return(comp_df)
}




