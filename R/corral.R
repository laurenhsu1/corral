#' compsvd: Compute Singular Value Decomposition (SVD)
#'
#' Computes SVD.
#'
#' @param mat matrix, pre-processed input; can be sparse or full (pre-processing can be performed using \code{\link{corral_preproc}}from this package)
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irlba::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components; Default is 10
#' @param ... other parameters for the svd and irlba functions can also be utilized. See those documentation for options.
#'
#' @return SVD result - a list with the following elements:
#' \describe{
#'     \item{d}{a vector of the diagonal singular values of the input \code{mat}. Note that using \code{svd} will result in the full set of singular values, while \code{irlba} will only compute the first \code{ncomp} singular values.}
#'     \item{u}{a matrix of with the left singular vectors of \code{mat} in the columns}
#'     \item{v}{a matrix of with the right singular vectors of \code{mat} in the columns}
#' }
#' @export
#' 
#' @importFrom irlba irlba
#'
#' @examples
#' mat <- matrix(sample(0:10, 500, replace=TRUE), ncol=25)
#' compsvd(mat)
#' compsvd(mat, method = 'svd')
#' compsvd(mat, method = 'irl', ncomp = 5)
compsvd <- function(mat, method = c('irl','svd')[1], ncomp = 10, ...){
  if(method == 'irl'){
    return(irlba::irlba(mat, nv = ncomp, ...))
  }
  else if(method == 'svd'){
    return(svd(mat, nv = ncomp, nu = ncomp, ...))
  }
  else {
    print('Your provided method was not understood. Used irlba.')
    return(irlba::irlba(preproc_mat, nv = ncomp, ...))
  }
}


#' Correspondence analysis preprocessing (not intended to be called directly)
#'
#' This function performs the row and column scaling pre-processing operations, prior to SVD, for the corral methods. See \code{\link{corral}} for single matrix correspondence analysis and \code{\link{corralm}} for multi-matrix correspondence analysis.
#'
#' @param inp matrix, numeric, counts or logcounts; can be sparse Matrix or matrix
#' @param rtype character indicating what type of residual should be computed; options are "indexed" and "standardized"; defaults to 'standardized'
#'
#' @return sparse matrix, processed for input to \code{compsvd} to finish CA routine
#' @export
#' 
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' mat <- matrix(sample(0:10, 500, replace=TRUE), ncol=25)
#' mat_corral <- corral_preproc(mat)
#' corral_output <- compsvd(mat_corral, ncomp = 5)
corral_preproc <- function(inp, rtype = c('standardized','indexed')[1]){
  if(!is(inp, "dgCMatrix")) {sp_mat <- Matrix::Matrix(inp, sparse = TRUE)}  # convert to a sparse matrix
  else {sp_mat <- inp}
  N <- sum(sp_mat)
  sp_mat <- sp_mat/N
  row.w <- Matrix::rowSums(sp_mat)
  col.w <- Matrix::colSums(sp_mat)
  sp_mat <- sp_mat/row.w
  sp_mat <- sweep(sp_mat, 2, col.w, "/") - 1
  if (any(is.na(sp_mat))) {
    sp_mat <- apply(sp_mat, c(1, 2), na2zero)
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

#' corral: Correspondence analysis on a single matrix
#'
#' corral can be used for dimensionality reduction to find a set of low-dimensional embeddings for a count matrix.
#'
#' @param inp matrix (or any matrix-like object that can be coerced using Matrix), numeric raw or lognormed counts (no negative values)
#' @inheritParams compsvd
#'
#' @return When run on a matrix, a list with the correspondence analysis matrix decomposition result:
#' \describe{
#'     \item{d}{a vector of the diagonal singular values of the input \code{mat} (from SVD output)}
#'     \item{u}{a matrix of with the left singular vectors of \code{mat} in the columns (from SVD output)}
#'     \item{v}{a matrix of with the right singular vectors of \code{mat} in the columns. When cells are in the columns, these are the cell embeddings. (from SVD output)}
#'     \item{SCu and SCv}{standard coordinates, left and right, respectively}
#'     \item{PCu and PCv}{principal coordinates, left and right, respectively}
#' }
#' 
#' @rdname corral
#' @export
#'
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' mat <- matrix(sample(0:10, 5000, replace=TRUE), ncol=50)
#' result <- corral_mat(mat)
#' result <- corral_mat(mat, method = 'svd')
#' result <- corral_mat(mat, method = 'irl', ncomp = 5)
#' 
corral_mat <- function(inp, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mat <- corral_preproc(inp,...)
  result <- compsvd(preproc_mat, method, ncomp,...)
  w <- get_weights(inp)
  result[['SCu']] <- sweep(result$u,1,sqrt(w$row.w),'/') # Standard coordinates
  result[['SCv']] <- sweep(result$v,1,sqrt(w$col.w),'/')
  result[['PCu']] <- sweep(result[['SCu']],2,result$d[1:dim(result$u)[2]],'*') # Principal coordinates
  result[['PCv']] <- sweep(result[['SCv']],2,result$d[1:dim(result$v)[2]],'*')
  return(result)
}

#' corral: Correspondence analysis on a single matrix (SingleCellExperiment)
#'
#' @param inp SingleCellExperiment; raw or lognormed counts (no negative values)
#' @inheritParams compsvd
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#' @param fullout boolean; whether the function will return the full \code{corral} output as a list, or a SingleCellExperiment; defaults to SingleCellExperiment (\code{FALSE}). To get back the \code{\link{corral_mat}}-style output, set this to \code{TRUE}.
#'
#' @return When run on a \code{\link{SingleCellExperiment}}, returns a SCE with the embeddings in the \code{reducedDim} slot \code{corral} (default). Also can return the same output as \code{\link{corral_mat}}.
#' 
#' @rdname corral
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
#' sce <- sce_full_Zhengmix4eq()[1:100,1:100]
#' result_1 <- corral_sce(sce)
#' result_2 <- corral_sce(sce, method = 'svd')
#' result_3 <- corral_sce(sce, method = 'irl', ncomp = 30, whichmat = 'logcounts')
#' 
#' # example on how to add UMAP based on corral above, with 'scater' package
#' library(scater)
#' result_1 <- runUMAP(result_1, dimred = 'corral', name = 'corral_UMAP')
#' result_1 <- runTSNE(result_1, dimred = 'corral', name = 'corral_TSNE')
#' 
corral_sce <- function(inp, method = c('irl','svd')[1], ncomp = 10, whichmat = 'counts', fullout = FALSE, ...){
  inp_mat <- SummarizedExperiment::assay(inp, whichmat)
  svd_output <- corral_mat(inp_mat, ...)
  if(fullout){
    return(svd_output)
  }
  else{
    SingleCellExperiment::reducedDim(inp,'corral') <- svd_output$v
    return(inp)
  }
}


#' corral: Correspondence analysis for a single matrix
#' 
#' \code{corral} is a wrapper for \code{\link{corral_mat}} and \code{\link{corral_sce}}, and can be called on any of the acceptable input types.
#'
#' @param inp matrix (any type), \code{SingleCellExperiment}, or \code{SummarizedExperiment}. If using \code{SingleCellExperiment} or \code{SummarizedExperiment}, then include the \code{whichmat} argument to specify which slot to use (defaults to \code{counts}).
#' @param ... 
#'
#' @return For matrix and \code{SummarizedExperiment} input, returns list with the correspondence analysis matrix decomposition result (u,v,d are the raw svd output; SCu and SCv are the standard coordinates; PCu and PCv are the principal coordinates)
#' @return
#' For \code{SummarizedExperiment} input, returns the same as for a matrix.
#' @rdname corral
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
#' sce <- sce_full_Zhengmix4eq()[1:100,1:100]
#' corral_sce <- corral(sce,whichmat = 'counts')
#' 
#' mat <- matrix(sample(0:10, 500, replace=TRUE), ncol=25)
#' corral_mat <- corral(mat, ncomp=5)
#' 
corral <- function(inp,...){
  if(is(inp,"SingleCellExperiment")){
    corral_sce(inp = inp, ...)
  }
  else if(is(inp,"SummarizedExperiment")){
    if(missing(whichmat)) {whichmat <- 'counts'}
    corral_mat(inp = SummarizedExperiment::assay(inp,whichmat),...)
  }
  else{
    corral_mat(inp = inp, ...)
  }
}
