#' Multi-table correspondence analysis (list of matrices)
#'
#' @param matlist (for \code{corralm_matlist}) list of input matrices; input matrices should be counts (raw or log). Matrices should be aligned row-wise by common features (either by sample or by gene)
#' @inheritParams compsvd
#' @return When run on a list of matrices, a list with the correspondence analysis matrix decomposition result, with indices corresponding to the concatenated matrices (in order of the list):
#' \describe{
#'     \item \code{d}: a vector of the diagonal singular values of the input \code{mat}
#'     \item \code{u}: a matrix of with the left singular vectors of \code{mat} in the columns
#'     \item \code{v}: a matrix of with the right singular vectors of \code{mat} in the columns. When cells are in the columns, these are the cell embeddings.
#' }
#' @rdname corralm
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' listofmats <- list(matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20),matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20))
#' result <- corralm_matlist(listofmats)
corralm_matlist <- function(matlist, method = c('irl','svd')[1], ncomp = 10, ...){
  .check_dims(matlist)
  preproc_mats <- lapply(matlist, corral_preproc, rtype = 'indexed')
  concatted <- list2mat(matlist = preproc_mats, direction = 'c')
  result <- compsvd(concatted, method = method, ncomp = ncomp)
  return(result)
}

#' Multi-table correspondence analysis (SingleCellExperiment)
#'
#' @param sce (for \code{corralm_sce}) SingleCellExperiment; containing the data to be integrated. Default is to use the counts, and to include all of the data in the integration. These can be changed by passing additional arguments. See \code{\link{sce2matlist}} function documentation for list of available parameters.
#' @param splitby character; name of the attribute from \code{colData} that should be used to separate the SCE
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#' @inheritParams compsvd
#' @param ... (additional arguments for methods)
#'
#' @return For SingleCellExperiment input, returns the SCE with embeddings in the reducedDim slot 'corralm' 
#' @rdname corralm
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' library(DuoClustering2018)
#' sce <- sce_full_Zhengmix4eq()[1:100,sample(1:3500,100,replace = FALSE)]
#' # for illustrative purposes only; would not actually use splitby = "phenoid"
#' # would instead use a batch or platform attribute
#' result <- corralm_sce(sce, splitby = 'phenoid')
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

#' Multi-table adaptation of correspondence analysis
#' 
#' This multi-table adaptation of correpondence analysis applies the same scaling technique and enables data alignment by finding a set of embeddings for each dataset within shared latent space.
#' 
#' \code{corralm} is a wrapper for \code{\link{corralm_matlist}} and \code{\link{corralm_sce}}, and can be called on any of the acceptable input types (see \code{inp} below).
#'
#' @param inp list of matrices (any type), \code{SingleCellExperiment}, or list of \code{SingleCellExperiment}s. If using \code{SingleCellExperiment}, then include the \code{whichmat} argument to specify which slot to use (defaults to \code{counts}).
#' @param ... (additional arguments for methods)
#'
#' @return For a list of \code{\link{SingleCellExperiment}}s, returns a list of the SCEs with the embeddings in the respective \code{reducedDim} slot 'corralm'
#' @rdname corralm
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' listofmats <- list(matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20),matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20))
#' corralm(listofmats)
#' 
#' library(DuoClustering2018)
#' sce <- sce_full_Zhengmix4eq()[1:100,sample(1:3500,100,replace = FALSE)]
#' # for illustrative purposes only; would not actually use splitby = "phenoid"
#' # would instead use a batch or platform attribute
#' result <- corralm(sce, splitby = 'phenoid')
#' 
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