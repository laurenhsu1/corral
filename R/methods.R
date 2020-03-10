# All the functions related to any specific step in the methods' workflows

library(Matrix)
library(SingleCellExperiment)
library(irlba)

#' Compute Singular Value Decomposition (SVD)
#'
#' @param preproc_mat matrix, pre-processed input; can be sparse or full [pre-processing can be performed using corral_preproc or pca_preproc from this package]
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irbla::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components (if using irlba); Default is 10
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
#' mat_coa <- coa_preproc(my_matrix)
#' coa_output <- compsvd(mat_coa)
corral_preproc <- function(inp_mat, rtype = c('standardized','indexed')[1]){
  if(!is(inp_mat, "dgCMatrix")) {sp_mat <- Matrix(inp_mat, sparse = TRUE)}  # convert to a sparse matrix
  else {sp_mat <- inp_mat}
  N <- sum(sp_mat)
  sp_mat <- sp_mat/N
  row.w <- rowSums(sp_mat)
  col.w <- colSums(sp_mat)
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


#' Principal components analysis (PCA) preprocessing
#'
#' @param inp_mat matrix of logcounts; can be sparse or full
#' @param type character, default, runs correlation-based (using 'corr' argument); can also run covariance-based (using 'cov' argument)
#'
#' @return preprocessed sparse matrix upon which compsvd can be run to complete PCA routine
#' @export
#'
#' @examples
#' mat_pca <- pca_preproc(my_matrix, type = 'cov')
#' pca_output <- compsvd(mat_pca)
pca_preproc <- function(inp_mat, type = 'corr'){
  if(type == 'corr'){
    scale_bool <- TRUE
    }
  else if(type == 'cov'){
    scale_bool <- FALSE
  }
  return(scale(inp_mat, scale = scale_bool, center = TRUE))
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
#' result <- coa(my_matrix)
#' result <- coa(my_matrix, method = 'svd')
#' result <- coa(my_matrix, method = 'irl', ncomp = 30)
corral <- function(inp_mat, method = c('irl','svd')[1], ncomp = 10, ...){
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
corralm <- function(mat_list, method = c('irl','svd')[1], ncomp = 10, ...){
  preproc_mats <- lapply(mat_list, corral_preproc, rtype = 'indexed')
  concatted <- list2mat(matlist = preproc_mats, direction = 'c')
  result <- compsvd(concatted, method = method, ncomp = ncomp)
  # add in the correspondence coordinates
  return(result)
}

