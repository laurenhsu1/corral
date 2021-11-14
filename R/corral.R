#' compsvd: Compute Singular Value Decomposition (SVD)
#'
#' Computes SVD.
#'
#' @param mat matrix, pre-processed input; can be sparse or full (pre-processing can be performed using \code{\link{corral_preproc}} from this package)
#' @param method character, the algorithm to be used for svd. Default is irl. Currently supports 'irl' for irlba::irlba or 'svd' for stats::svd
#' @param ncomp numeric, number of components; Default is 30
#' @param ... (additional arguments for methods)
#'
#' @return SVD result - a list with the following elements:
#' \describe{
#'     \item{\code{d}}{a vector of the diagonal singular values of the input \code{mat}. Note that using \code{svd} will result in the full set of singular values, while \code{irlba} will only compute the first \code{ncomp} singular values.}
#'     \item{\code{u}}{a matrix of with the left singular vectors of \code{mat} in the columns}
#'     \item{\code{v}}{a matrix of with the right singular vectors of \code{mat} in the columns}
#'     \item{\code{eigsum}}{sum of the eigenvalues, for calculating percent variance explained}
#' }
#' @export
#' 
#' @importFrom irlba irlba
#'
#' @examples
#' mat <- matrix(sample(0:10, 2500, replace=TRUE), ncol=50)
#' compsvd(mat, method = 'irl', ncomp = 5)
compsvd <- function(mat, method = c('irl','svd'), ncomp = 30, ...){
  method <- match.arg(method, c('irl','svd'))
  ncomp <- .check_ncomp(mat, ncomp)
  if(method == 'irl'){
    result <- irlba::irlba(mat, nv = ncomp, ...)
  }
  else if(method == 'svd'){
    result <- svd(mat, nv = ncomp, nu = ncomp, ...)
  }
  else {
    print('Provided method was not understood; used irlba.')
    result <- irlba::irlba(mat, nv = ncomp, ...)
  }
  result[['eigsum']] <- sum(mat^2)
  if(!is.null(rownames(mat))){rownames(result$u) <- rownames(mat)}
  if(!is.null(colnames(mat))){rownames(result$v) <- colnames(mat)}
  
  # add percent variance explained
  pct_var_exp <- t(data.frame('percent.Var.explained' = result$d^2 / result$eigsum))
  colnames(pct_var_exp) <- paste0(rep('PC',ncomp),seq(1,ncomp,1))
  result$pct_var_exp <- rbind(pct_var_exp,t(data.frame('cumulative.Var.explained' = cumsum(pct_var_exp[1,]))))
  
  return(result)
}

#' Apply a variance stabilizing transformation
#' 
#' Prior to running CA, there is an option to apply a variance stabilizing transformation. This function can be called explicitly or used with the `vst_mth` argument in \code{corral} and \code{corral_preproc}.
#'
#' @param inp matrix, numeric, counts or logcounts; can be sparse Matrix or matrix
#' @param transform character indicating which method should be applied. Defaults to the square root transform (`"sqrt"`). Other options include `"freemantukey"` and `"anscombe"`.
#'
#' @return variance-stabilized matrix; sparse if possible
#' @export
#'
#' @examples
#' x <- as.matrix(rpois(100, lambda = 50), ncol = 10)
#' vst_x <- var_stabilize(x)
var_stabilize <- function(inp, transform = c('sqrt','freemantukey','anscombe')){
  transform <- match.arg(transform, c('sqrt','freemantukey','anscombe'))
  .check_vals(inp)
  if(transform == 'sqrt'){
    return(inp^.5)
  }
  else if(transform == 'freemantukey'){
    return(inp^.5 + (inp + 1)^.5)
  }
  else if(transform == 'anscombe'){
    return(2 * (inp + 3/8)^.5)
  }
}


#' Preprocess a matrix for SVD to perform Correspondence Analysis (CA)
#'
#' This function performs the row and column scaling pre-processing operations, prior to SVD, for the corral methods. See \code{\link{corral}} for single matrix correspondence analysis and \code{\link{corralm}} for multi-matrix correspondence analysis.
#' 
#' In addition to standard CA (SVD of Pearson residuals), which is the default setting, this function also provides options to customize processing steps -- including the variations described in the corresponding manuscript for \code{corral}. Parameters that may be changed:
#' \describe{
#'     \item{Residual (`rtype`)}{Which analytic residual should be calculated? The default is standardized Pearson residuals. Freeman-Tukey is a more robust choice for handling overdispersed scRNAseq counts.}
#'     \item{Variance stabilization method (`vst_mth`)}{Should a variance stabilizing transformation be applied to the count matrix before computing residuals? Defaults to no.}
#'     \item{Power deflation $\alpha$ (`powdef_alpha`)}{Should a power deflation transformation be applied after computing residuals? Defaults to no. If using, set $\alpha \in (0,1)$. For a "soft" smoothing effect, we suggest $\alpha \in [0.9,0.99]$}
#'     \item{Trimming-based smoothing option (`smooth`)}{As an alternative to power deflation, this option provides a more subtle correction: the most extreme 1% of values (distributed equally between both tails) are set to the 0.5 percentile and 99.5 percentile values, respectively}
#' }
#'
#' @param inp matrix, numeric, counts or logcounts; can be sparse Matrix or matrix
#' @param rtype character indicating what type of residual should be computed; options are `"indexed"`, `"standardized"` (or `"pearson"` is equivalent), `"freemantukey"`, and `"hellinger"`; defaults to `"standardized"` for \code{\link{corral}} and `"indexed"` for \code{\link{corralm}}. `"indexed"`, `"standardized"`, and `"freemantukey"` compute the respective chi-squared residuals and are appropriate for count data. The `"hellinger"` option is appropriate for continuous data. 
#' @param vst_mth character indicating whether a variance-stabilizing transform should be applied prior to calculating chi-squared residuals; defaults to `"none"`
#' @param powdef_alpha numeric for the power that should be applied if using power deflation. Must be in (0,1), and if provided a number outside this range, will be ignored. Defaults to `NULL` which does not perform this step.
#' @param row.w numeric vector; Default is \code{NULL}, to compute row.w based on \code{inp}. Use this parameter to replace computed row weights with custom row weights
#' @param col.w numeric vector; Default is \code{NULL}, to compute col.w based on \code{inp}. Use this parameter to replace computed column weights with custom column weights
#' @param smooth logical; Whether or not to perform the additional smoothing step with `trim_matdist`. Default is \code{FALSE}. Incompatible with `powdef_alpha`, so that parameter takes precedence over this one.
#' @param ... (additional arguments for methods)
#'
#' @return matrix, processed for input to \code{compsvd} to finish CA routine
#' @export
#' 
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' mat <- matrix(sample(0:10, 500, replace=TRUE), ncol=25)
#' mat_corral <- corral_preproc(mat)
#' corral_output <- compsvd(mat_corral, ncomp = 5)
corral_preproc <- function(inp, rtype = c('standardized','indexed','hellinger','freemantukey','pearson'), vst_mth = c('none','sqrt','freemantukey','anscombe'), powdef_alpha = NULL,  row.w = NULL, col.w = NULL, smooth = FALSE, ...){
  rtype <- match.arg(rtype, c('standardized','indexed','hellinger','freemantukey','pearson'))
  vst_mth <- match.arg(vst_mth, c('none','sqrt','freemantukey','anscombe'))
  if(!is.null(powdef_alpha)){
    if(powdef_alpha > 1 | powdef_alpha < 0){
      powdef_alpha <- NULL
      cat('Invalid choice for power deflation parameter alpha (powdef_alpha) so ignoring; should be in (0,1)')
    }
  }
  if(!is(inp, "dgCMatrix")){
    x_mat <- Matrix::Matrix(inp, sparse = TRUE)
  } else {x_mat <- inp}
  N <- sum(x_mat)
  
  if(vst_mth != 'none'){
    x_mat <- var_stabilize(x_mat, transform = vst_mth)
  }
  
  p_mat <- x_mat/N
  w <- get_weights(inp)
  if(is.null(row.w)) {row.w <- w$row.w}
  if(is.null(col.w)) {col.w<- w$col.w}
  else {col.w <- col.w / sum(col.w)}
  
  res <- NULL
  
  if(rtype == 'hellinger'){
   p_mat <- p_mat / row.w
   res <- sqrt(p_mat)
  }
  if(rtype == 'freemantukey'){
    expectedp <- row.w %*% t(col.w)
    res <- p_mat^.5 + (p_mat + 1/N)^.5 - (4*expectedp + 1/N)^.5
  }
  
  p_mat <- p_mat/row.w
  p_mat <- sweep(p_mat, 2, col.w, "/") - 1
  if (any(is.na(p_mat))) {
    p_mat <- na2zero(p_mat)
  }
  
  if (rtype == 'indexed'){
    res <- p_mat
  }
  else if (rtype %in% c('standardized','pearson')){
    p_mat <- p_mat * sqrt(row.w)
    p_mat <- sweep(p_mat, 2, sqrt(col.w), "*")
    res <- p_mat
  }
  
  if(!is.null(powdef_alpha)){
    return(res^powdef_alpha)
  }
  else if (smooth){
    return(trim_matdist(res, ...))
  }
  return(res)
}


#' corral: Correspondence analysis on a single matrix
#'
#' corral can be used for dimension reduction to find a set of low-dimensional embeddings for a count matrix.
#'
#' @param inp matrix (or any matrix-like object that can be coerced using Matrix), numeric raw or lognormed counts (no negative values)
#' @param row.w numeric vector; the row weights to use in chi-squared scaling. Defaults to `NULL`, in which case row weights are computed from the input matrix.
#' @param col.w numeric vector; the column weights to use in chi-squared scaling. For instance, size factors could be given here. Defaults to `NULL`, in which case column weights are computed from the input matrix.
#' @inheritParams compsvd
#' @inheritParams corral_preproc
#'
#' @return When run on a matrix, a list with the correspondence analysis matrix decomposition result:
#' \describe{
#'     \item{\code{d}}{a vector of the diagonal singular values of the input \code{mat} (from SVD output)}
#'     \item{\code{u}}{a matrix of with the left singular vectors of \code{mat} in the columns (from SVD output)}
#'     \item{\code{v}}{a matrix of with the right singular vectors of \code{mat} in the columns. When cells are in the columns, these are the cell embeddings. (from SVD output)}
#'     \item{\code{eigsum}}{sum of the eigenvalues for calculating percent variance explained}
#'     \item{\code{SCu and SCv}}{standard coordinates, left and right, respectively}
#'     \item{\code{PCu and PCv}}{principal coordinates, left and right, respectively}
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
#' result <- corral_mat(mat, method = 'irl', ncomp = 5)
#' 
corral_mat <- function(inp, method = c('irl','svd'),ncomp = 30, row.w = NULL, col.w = NULL, rtype = c('standardized','indexed','hellinger','freemantukey','pearson'), vst_mth = c('none','sqrt','freemantukey','anscombe'), ...){
  method <- match.arg(method, c('irl','svd'))
  rtype <- match.arg(rtype, c('standardized','indexed','hellinger','freemantukey','pearson'))
  vst_mth <- match.arg(vst_mth, c('none','sqrt','freemantukey','anscombe'))
  preproc_mat <- corral_preproc(inp, row.w = row.w, col.w = col.w, rtype = rtype, vst_mth = vst_mth, ...)
  result <- compsvd(preproc_mat, method, ncomp, ...)
  w <- get_weights(inp)
  if(is.null(row.w)) {row.w <- w$row.w}
  if(is.null(col.w)) {col.w<- w$col.w}
  result[['SCu']] <- sweep(result$u,1,sqrt(row.w),'/') # Standard coordinates
  result[['SCv']] <- sweep(result$v,1,sqrt(col.w),'/')
  result[['PCu']] <- sweep(result[['SCu']],2,result$d[seq(1,ncol(result$u),1)],'*') # Principal coordinates
  result[['PCv']] <- sweep(result[['SCv']],2,result$d[seq(1,ncol(result$v),1)],'*')
  class(result) <- c(class(result),'corral')
  return(result)
}

#' corral: Correspondence analysis on a single matrix (SingleCellExperiment)
#'
#' @param inp SingleCellExperiment; raw or lognormed counts (no negative values)
#' @inheritParams compsvd
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#' @param fullout boolean; whether the function will return the full \code{corral} output as a list, or a SingleCellExperiment; defaults to SingleCellExperiment (\code{FALSE}). To get back the \code{\link{corral_mat}}-style output, set this to \code{TRUE}.
#' @param subset_row numeric, character, or boolean vector; the rows to include in corral, as indices (numeric), rownames (character), or with booleans (same length as the number of rows in the matrix). If this parameter is \code{NULL}, then all rows will be used.
#'
#' @return When run on a \code{\link{SingleCellExperiment}}, returns a SCE with the embeddings (PCv from the full corral output) in the \code{reducedDim} slot \code{corral} (default). Also can return the same output as \code{\link{corral_mat}} when \code{fullout} is set to \code{TRUE}.
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
#' 
corral_sce <- function(inp, method = c('irl','svd'), ncomp = 30, whichmat = 'counts', fullout = FALSE, subset_row = NULL, ...){
  method <- match.arg(method, c('irl','svd'))
  inp_mat <- SummarizedExperiment::assay(inp, whichmat)
  if(!is.null(subset_row)){
    inp_mat <- inp_mat[subset_row,]
  }
  svd_output <- corral_mat(inp_mat, ...)
  if(fullout){
    return(svd_output)
  }
  else{
    SingleCellExperiment::reducedDim(inp,'corral') <- svd_output$PCv
    return(inp)
  }
}


#' corral: Correspondence analysis for a single matrix
#' 
#' \code{corral} is a wrapper for \code{\link{corral_mat}} and \code{\link{corral_sce}}, and can be called on any of the acceptable input types.
#'
#' @param inp matrix (any type), \code{SingleCellExperiment}, or \code{SummarizedExperiment}. If using \code{SingleCellExperiment} or \code{SummarizedExperiment}, then include the \code{whichmat} argument to specify which slot to use (defaults to \code{counts}).
#' @param ... (additional arguments for methods)
#'
#' @return For matrix and \code{SummarizedExperiment} input, returns list with the correspondence analysis matrix decomposition result (u,v,d are the raw svd output; SCu and SCv are the standard coordinates; PCu and PCv are the principal coordinates)
#' @return
#' For \code{SummarizedExperiment} input, returns the same as for a matrix.
#' @rdname corral
#' @export
#' 
#' @importFrom irlba irlba
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
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
  } else if(is(inp,"SummarizedExperiment")){
    if(missing(whichmat)) {whichmat <- 'counts'}
    corral_mat(inp = SummarizedExperiment::assay(inp,whichmat),...)
  } else{
    corral_mat(inp = inp, ...)
  }
}


#' Print method for S3 object corral
#'
#' @param x (print method) corral object; the list output from \code{corral_mat}
#'
#' @rdname corral
#'
#' @return .
#' @export
#'
#' @examples
#' mat <- matrix(sample(1:100, 10000, replace = TRUE), ncol = 100)
#' corral(mat)
print.corral <- function(x,...){
  inp <- x
  pct_var_exp <- inp$pct_var_exp
  ncomp <- length(inp$d)
  cat('corral output summary===========================================\n')
  cat('  Output "list" includes standard coordinates (SCu, SCv),\n')
  cat('  principal coordinates (PCu, PCv), & SVD output (u, d, v)\n')
  cat('Variance explained----------------------------------------------\n')
  print(round(pct_var_exp[,seq(1,min(8,ncomp),1)],2))
  cat('\n')
  cat('Dimensions of output elements-----------------------------------\n')
  cat('  Singular values (d) :: ')
  cat(ncomp)
  cat('\n  Left singular vectors & coordinates (u, SCu, PCu) :: ')
  cat(dim(inp$u))
  cat('\n  Right singular vectors & coordinates (v, SCv, PCv) :: ')
  cat(dim(inp$v))
  cat('\n  See corral help for details on each output element.')
  cat('\n  Use plot_embedding to visualize; see docs for details.')
  cat('\n================================================================\n')
}

