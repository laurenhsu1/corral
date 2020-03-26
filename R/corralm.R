#' Multi-table correspondence analysis
#'
#' @param mat_list list of input matrices; input matrices should be counts (raw or log). Matrices should be aligned row-wise by common features (either by sample or by gene)
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components (if using irlba); Default is 10
#'
#' @return list with the correspondence analysis matrix decomposition result
#' @export
#'
#' @examples
#' corralm(mat_list)
corralm_matlist <- function(mat_list, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mats <- lapply(mat_list, corral_preproc, rtype = 'indexed')
  concatted <- list2mat(matlist = preproc_mats, direction = 'c')
  result <- compsvd(concatted, method = method, ncomp = ncomp)
  # TODO add in the correspondence coordinates
  return(result)
}

#' Multi-table correspondence analysis
#'
#' @param sce SingleCellExperiment; containing the data to be integrated. Default is to use the counts, and to include all of the data in the integration. These can be changed by passing additional arguments. See sce2matlist function documentation for list of available parameters.
#' @param splitby character; name of the attribute from colData that should be used to separate the SCE
#' @param method character; the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric; number of components (if using irlba); Default is 10
#' @param ... for additional parameters, see the sce2matlist function
#'
#' @return list with the correspondence analysis matrix decomposition result
#' @export
#'
#' @examples
#' result <- corralm_sce(my_sce, 'platform')
#' 
#' #' # example on how to add UMAP/tsne based on corralm above, with 'scater' package
#' library(scater)
#' result <- runUMAP(result, dimred = 'corralm', name = 'corralm_UMAP')
#' result <- runTSNE(result, dimred = 'corralm', name = 'corralm_TSNE')
#' 
corralm_sce <- function(sce, splitby, method = c('irl','svd')[1], ncomp = 10, ...){
  mat_list <- sce2matlist(sce, splitby = splitby)
  svd_output <- corralm_matlist(mat_list, method = method, ncomp = ncomp)
  SingleCellExperiment::reducedDim(sce, 'corralm') <- svd_output$v
  return(sce)
}

