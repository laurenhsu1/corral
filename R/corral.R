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
#' @param inp_mat matrix, counts or logcounts; can be sparse or full
#' @param rtype character indicating what type of residual should be computed; options are "indexed" and "standardized"; defaults to 'standardized'
#'
#' @return sparse matrix, preprocessed, upon which compsvd can be run to complete CA routine
#' @export
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
#' @examples
#' result <- corral_mat(my_matrix)
#' result <- corral_mat(my_matrix, method = 'svd')
#' result <- corral_mat(my_matrix, method = 'irl', ncomp = 30)
corral_mat <- function(inp_mat, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mat <- corral_preproc(inp_mat)
  result <- compsvd(preproc_mat, method, ncomp)
  w <- get_w(inp_mat)
  result[['SCu']] <- sweep(result$u,1,sqrt(w$row.w),'/') # Standard coordinates
  result[['SCv']] <- sweep(result$v,1,sqrt(w$col.w),'/')
  result[['PCu']] <- sweep(result[['SCu']],2,result$d,'*') # Principal coordinates
  result[['PCv']] <- sweep(result[['SCv']],2,result$d,'*')
  # add in the correspondence coordinates
  return(result)
}

#' Correspondence analysis on a single matrix (SingleCellExperiment)
#'
#' @param inp_sce SingleCellExperiment; raw or lognormed counts (no negative values)
#' @param method character; the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric; number of components (if using irlba); Default is 10
#' @param whichmat accessor function; defaults to counts, can also use logcounts or normcounts if available
#'
#' @return SCE with the SVD output in the reducedDim portion under 'corral'
#' @export
#'
#' @examples
#' result <- corral_sce(my_sce)
#' result <- corral_sce(my_sce, method = 'svd')
#' result <- corral_sce(my_sce, method = 'irl', ncomp = 30, whichmat = logcounts)
#' 
#' # example on how to add UMAP based on corral above, with 'scater' package
#' library(scater)
#' result <- runUMAP(result, dimred = 'corral', name = 'corral_UMAP')
#' result <- runTSNE(result, dimred = 'corral', name = 'corral_TSNE')
#' 
corral_sce <- function(inp_sce, method = c('irl','svd')[1], ncomp = 10, whichmat = counts,...){
  inp_mat <- whichmat(inp_sce)
  svd_output <- corral_mat(inp_mat)
  SingleCellExperiment::reducedDim(inp_sce,'corral') <- svd_output$v
  return(inp_sce)
}
