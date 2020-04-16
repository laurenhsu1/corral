#' Multi-table correspondence analysis (list of matrices)
#'
#' @param matlist list of input matrices; input matrices should be counts (raw or log). Matrices should be aligned row-wise by common features (either by sample or by gene)
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components (if using irlba); Default is 10
#'
#' @return list with the correspondence analysis matrix decomposition result
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' corralm(matlist)
corralm_matlist <- function(matlist, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mats <- lapply(matlist, corral_preproc, rtype = 'indexed')
  concatted <- list2mat(matlist = preproc_mats, direction = 'c')
  result <- compsvd(concatted, method = method, ncomp = ncomp)
  # TODO add in the correspondence coordinates
  return(result)
}

#' Multi-table correspondence analysis (SingleCellExperiment)
#'
#' @param sce SingleCellExperiment; containing the data to be integrated. Default is to use the counts, and to include all of the data in the integration. These can be changed by passing additional arguments. See \code{\link{sce2matlist}} function documentation for list of available parameters.
#' @param splitby character; name of the attribute from \code{colData} that should be used to separate the SCE
#' @param method character; the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric; number of components (if using irlba); Default is 10
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#' @param ... for additional parameters, see \code{\link{sce2matlist}}
#'
#' @return list with the correspondence analysis matrix decomposition result
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' result <- corralm_sce(my_sce, 'platform')
#' 
#' #' # example on how to add UMAP/tsne based on corralm above, with 'scater' package
#' library(scater)
#' result <- runUMAP(result, dimred = 'corralm', name = 'corralm_UMAP')
#' result <- runTSNE(result, dimred = 'corralm', name = 'corralm_TSNE')
#' 
corralm_sce <- function(sce, splitby, method = c('irl','svd')[1], ncomp = 10, whichmat = 'counts',...){
  if(missing(splitby)) {stop('If performing multi-table analysis with a single SCE, the splitby variable must be specified. \nUse corral to analyze as a single table.')}
  mat_list <- sce2matlist(sce, splitby = splitby, whichmat = whichmat)
  svd_output <- corralm_matlist(mat_list, method = method, ncomp = ncomp)
  SingleCellExperiment::reducedDim(sce, 'corralm') <- svd_output$v
  return(sce)
}

#' Multi-table correspondence analysis
#' 
#' This is a wrapper for \code{\link{corralm_matlist}} and \code{\link{corral_sce}}. See those functions for the list of possible parameters.
#'
#' @param inp list of matrices (any type), \code{SingleCellExperiment}, or list of \code{SingleCellExperiment}s. If using \code{SingleCellExperiment}, then include the \code{whichmat} argument to specify which slot to use (defaults to \code{counts}).
#' @param ... 
#'
#' @return For a list of matrices input, returns list with the correspondence analysis matrix decomposition result
#' For SingleCellExperiment input, returns the SCE with embeddings in the reducedDim slot 'corral'
#' For a list of SingleCellExperiments, returns a list of the SCEs with the embeddings in the respective reducedDim slot 'corral'
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
corralm <- function(inp,...){
  if(is(inp,'SingleCellExperiment')){
    corralm_sce(sce = inp, ...)
  }
  else if(is(inp,'list')){
    if(is(inp[[1]],'SingleCellExperiment')){
      matlist <- lapply(inp, assay, whichmat)
      res <- corralm_matlist(matlist = matlist, ...)
      add_embeddings2scelist(scelist = inp, embeddings = res$v)
    }
    else{
      corralm_matlist(matlist = inp, ...) 
    }
  }
}