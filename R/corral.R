#' Compute Singular Value Decomposition (SVD)
#'
#' @param preproc_mat matrix, pre-processed input; can be sparse or full [pre-processing can be performed using corral_preproc or pca_preproc from this package]
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components (if using irlba); Default is 10
#' @param ... other parameters for the svd and irlba functions can also be utilized. See those documentation for options.
#'
#' @return SVD result. List of matrices, dvu, d is  diagonal matrix with the eigenvalues
#' @export
#' 
#' @importFrom irlba irlba
#'
#' @examples
#' compsvd(my_matrix)
#' compsvd(my_matrix, method = 'svd')
#' compsvd(my_matrix, method = 'irl', ncomp = 20)
compsvd <- function(preproc_mat, method = c('irl','svd')[1], ncomp = 10, ...){
  if(method == 'irl'){
    return(irlba::irlba(preproc_mat, nv = ncomp, ...))
  }
  else if(method == 'svd'){
    return(svd(preproc_mat, nv = ncomp, ...))
  }
}


#' Correspondence analysis preprocessing (not intended to be called directly)
#'
#' This function performs the row and column scaling pre-processing operations, prior to SVD, for the corral methods. See \code{\link{corral}} for single matrix correspondence analysis and \code{\link{corralm}} for multi-matrix correspondence analysis.
#'
#' @param inp_mat matrix, counts or logcounts; can be sparse or full
#' @param rtype character indicating what type of residual should be computed; options are "indexed" and "standardized"; defaults to 'standardized'
#'
#' @return sparse matrix, preprocessed, upon which compsvd can be run to complete CA routine
#' @export
#' 
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' mat_corral <- corral_preproc(my_matrix)
#' corral_output <- compsvd(mat_corral)
corral_preproc <- function(inp_mat, rtype = c('standardized','indexed')[1]){
  if(!is(inp_mat, "dgCMatrix")) {sp_mat <- Matrix::Matrix(inp_mat, sparse = TRUE)}  # convert to a sparse matrix
  else {sp_mat <- inp_mat}
  N <- sum(sp_mat)
  sp_mat <- sp_mat/N
  row.w <- Matrix::rowSums(sp_mat)
  col.w <- Matrix::colSums(sp_mat)
  sp_mat <- sp_mat/row.w
  sp_mat <- sweep(sp_mat, 2, col.w, "/") - 1
  if (any(is.na(sp_mat))) {
    sp_mat <- apply(sp_mat, c(1, 2), fix_na)
  }
  if (rtype == 'indexed'){
    return(sp_mat)
  }
  else if (rtype == 'standardized'){
    sp_mat <- sp_mat * sqrt(row.w)
    sp_mat <- sweep(sp_mat, 2, sqrt(col.w), "*")
    return(sp_mat)
  }
}

#' Correspondence analysis on a single matrix
#'
#' @param inp_mat matrix, raw or lognormed counts (no negative values)
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components (if using irlba); Default is 10
#'
#' @return list with the correspondence analysis matrix decomposition result (u,v,d are the raw svd output; SCu and SCv are the standard coordinates; PCu and PCv are the principal coordinates)
#' @export
#'
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' result <- corral_mat(my_matrix)
#' result <- corral_mat(my_matrix, method = 'svd')
#' result <- corral_mat(my_matrix, method = 'irl', ncomp = 30)
corral_mat <- function(inp_mat, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mat <- corral_preproc(inp_mat,...)
  result <- compsvd(preproc_mat, method, ncomp)
  w <- get_w(inp_mat)
  result[['SCu']] <- sweep(result$u,1,sqrt(w$row.w),'/') # Standard coordinates
  result[['SCv']] <- sweep(result$v,1,sqrt(w$col.w),'/')
  result[['PCu']] <- sweep(result[['SCu']],2,result$d,'*') # Principal coordinates
  result[['PCv']] <- sweep(result[['SCv']],2,result$d,'*')
  return(result)
}

#' Correspondence analysis on a single matrix (SingleCellExperiment)
#'
#' @param inp_sce SingleCellExperiment; raw or lognormed counts (no negative values)
#' @param method character; the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric; number of components (if using irlba); Default is 10
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#' @param fullout boolean; whether the function will return the full \code{corral} output as a list, or a SingleCellExperiment; defaults to SingleCellExperiment (\code{FALSE}). To get back the \code{\link{corral_mat}}-style output, set this to \code{TRUE}.
#'
#' @return (default) SCE with the embeddings in the \code{reducedDim} slot \code{corral}. Also can return the same output as \code{\link{corral_mat}}.
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' result <- corral_sce(my_sce)
#' result <- corral_sce(my_sce, method = 'svd')
#' result <- corral_sce(my_sce, method = 'irl', ncomp = 30, whichmat = 'logcounts')
#' 
#' # example on how to add UMAP based on corral above, with 'scater' package
#' library(scater)
#' result <- runUMAP(result, dimred = 'corral', name = 'corral_UMAP')
#' result <- runTSNE(result, dimred = 'corral', name = 'corral_TSNE')
#' 
corral_sce <- function(inp_sce, method = c('irl','svd')[1], ncomp = 10, whichmat = 'counts', fullout = FALSE, ...){
  inp_mat <- assay(inp_sce, whichmat)
  svd_output <- corral_mat(inp_mat, ...)
  if(fullout){
    return(svd_output)
  }
  else{
    SingleCellExperiment::reducedDim(inp_sce,'corral') <- svd_output$v
    return(inp_sce)
  }
}


#' Correspondence analysis for a single matrix
#' 
#' This is a wrapper for \code{\link{corral_mat}} and \code{\link{corral_sce}}. See those functions for the list of possible parameters.
#'
#' @param inp matrix (any type), \code{SingleCellExperiment}, or \code{SummarizedExperiment}. If using \code{SingleCellExperiment} or \code{SummarizedExperiment}, then include the \code{whichmat} argument to specify which slot to use (defaults to \code{counts}).
#' @param ... 
#'
#' @return For matrix and \code{SummarizedExperiment} input, returns list with the correspondence analysis matrix decomposition result (u,v,d are the raw svd output; SCu and SCv are the standard coordinates; PCu and PCv are the principal coordinates)
#' For SingleCellExperiment input, returns the SCE with embeddings in the reducedDim portion under 'corral'
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
corral <- function(inp,...){
  if(is(inp,"SingleCellExperiment")){
    corral_sce(inp_sce = inp, ...)
  }
  else if(is(inp,"SummarizedExperiment")){
    if(missing(whichmat)) {whichmat <- 'counts'}
    corral_mat(inp_mat = assay(inp,whichmat),...)
  }
  else{
    corral_mat(inp_mat = inp, ...)
  }
}