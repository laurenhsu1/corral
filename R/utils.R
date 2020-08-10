#' Set na to 0
#'
#' @param x matrix of values for which na values should be changed to 0
#'
#' @return matrix, where na values are set to 0
#' @export
#'
#' @examples
#' x <- matrix(sample(0:10, 5000, replace = TRUE), ncol = 25)
#' x[sample(1:5000, 10)] <- NA
#' 
#' na2zero(x)
na2zero <- function(x) {
  func <- function(x){
    if (is.na(x)) return(0)
    else return(x)}
  
  return(apply(x, c(1,2), func))
}

#' List to Matrix
#'
#' @param matlist list of matrices to concatenate 
#' @param direction character, r or c, to indicate whether should be row-wise (i.e., rbind to match on columns) or column-wise (i.e., cbind to match on rows). Defaults to columnwise (matching on rows) to match convention of SingleCellExperiments
#'
#' @return matrix
#' @export
#'
#' @examples
#' listofmats <- list(matrix(sample(seq(0,20,1),100,replace = TRUE),nrow = 10),
#'                    matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 10))
#' newmat <- list2mat(listofmats) # to "cbind" them
#' listofmats_t <- lapply(listofmats,t)
#' newmat_t <- list2mat(listofmats_t, 'r') # to "rbind" them
list2mat <- function(matlist, direction = c('c','r')[1]){
  if (direction == 'r') concatted <- Reduce(rbind, matlist)
  if (direction == 'c') concatted <- Reduce(cbind, matlist)
  return(concatted)
}

#' Compute percent of variance explained
#'
#' @param thissvd list outputted from an svd function (svd, irlba; can also take output from \code{\link{corral_mat}} and \code{\link{corralm_matlist}})
#' @param preproc_mat matrix of pre-processed values (optional) - important to include if the svd is only partial as this is used to compute the sum of eigenvalues
#'
#' @return vector of percent variance explained values, indexed by PC
#' @export
#'
#' @examples
#' mat <- matrix(sample(seq(0,20,1),100,replace = TRUE),nrow = 10)
#' my_svd <- svd(mat)
#' get_pct_var_exp_svd(my_svd) # this works if my_svd is a full svd
#' my_irl <- irlba::irlba(mat,nv = 2)
#' get_pct_var_exp_svd(my_irl, preproc_mat = mat) # ... otherwise use this 
get_pct_var_exp_svd <- function(thissvd, preproc_mat = thissvd$d){
  denom <- sum(preproc_mat^2)
  return(thissvd$d^2 / denom)
}


#' Get weights: computes row weights and column weights
#'
#' @param inp_mat matrix for which weights should be calculated (sparse or full)
#'
#' @return list of 2 elements: 'row.w' and 'col.w' contain the row and column weights respectively
#' @export
#' 
#' @importFrom Matrix Matrix rowSums colSums
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' mat <- matrix(sample(seq(0,20,1),100,replace = TRUE),nrow = 10)
#' ws <- get_weights(mat)
get_weights <- function(inp_mat){
  w <- list()
  if(!is(inp_mat, "dgCMatrix")) {sp_mat <- Matrix(inp_mat, sparse = TRUE)}  # convert to a sparse matrix
  else {sp_mat <- inp_mat}
  N <- sum(sp_mat)
  sp_mat <- sp_mat/N
  w[['row.w']] <- Matrix::rowSums(sp_mat)
  w[['col.w']] <- Matrix::colSums(sp_mat)
  return(w)
}


#' SingleCellExperiment to list of matrices
#'
#' @param sce SingleCellExperiment that is to be separated into list of count matrices
#' @param splitby character; name of the attribute from colData that should be used to separate the SCE
#' @param to_include (optional) character vector; determines which values from the "splitby" column will be included in the outputted matlist. NULL is the default, and will result in selecting all elements
#' @param whichmat character; defaults to \code{counts}, can also use \code{logcounts} or \code{normcounts} if stored in the \code{sce} object
#'
#' @return list of matrices
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' library(DuoClustering2018)
#' sce <- sce_full_Zhengmix4eq()
#' matlist <- sce2matlist(sce = sce, splitby = 'phenoid', whichmat = 'logcounts')
sce2matlist <- function(sce,splitby,to_include = NULL,whichmat = 'counts'){
  if(is.null(to_include)){
    to_include <- unique(as.character(colData(sce)[,splitby]))
  }
  matlist <- list()
  countmatrix <- SummarizedExperiment::assay(sce, whichmat)
  for(gn in to_include){
    matlist[[gn]] <- countmatrix[,which(colData(sce)[,splitby] == gn)]
  }
  return(matlist)
}



#' @keywords internal
#' @import SingleCellExperiment
.indsbysplitby <- function(sce, splitby, to_include = NULL){
  if(is.null(to_include)){
    to_include <- unique(as.character(colData(sce)[,splitby]))
  }
  inds <- c()
  for(gn in to_include){
    inds <- c(inds, which(colData(sce)[,splitby] == gn))
  }
  return(inds)
}

#' Add embeddings to list of SCEs
#'
#' @param scelist list of SingleCellExperiments; to which the corresponding embeddings should be added
#' @param embeddings matrix; the embeddings outputted from a dimension reduction, e.g. \code{\link{corralm}}. Rows in this table correspond to columns in the SCEs in \code{scelist} (if all the SCEs were column-bound), and row indices should correspond to cells.
#' @param slotname character; name of the slot for the reduced dim embedding; defaults to \code{corral}
#'
#' @return list of SingleCellExperiments with respective embeddings stored in them
#' @export
#' 
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' library(DuoClustering2018)
#' sce <- sce_full_Zhengmix4eq()
#' scelist <- list(sce,sce)
#' embeddings <- matrix(sample(seq(0,20,1),dim(sce)[2]*6,replace = TRUE),nrow = dim(sce)[2]*2)
#' scelist <- add_embeddings2scelist(scelist, embeddings)
add_embeddings2scelist <- function(scelist,embeddings, slotname = 'corralm'){
  if(!all_are(scelist,'SingleCellExperiment')) {
    cat('Warning! You may have non-SCE elements in the list.')
  }
  dimvec <- unlist(lapply(scelist,ncol))
  for (i in seq(1,length(scelist),1)){
    if(i == 1) {start_ind <- 1}
    else{start_ind <- (sum(dimvec[seq(1,i-1,1)]) + 1)}
    
    end_ind <- start_ind + dimvec[i] - 1
    
    SingleCellExperiment::reducedDim(scelist[[i]],slotname) <- embeddings[seq(start_ind,end_ind,1),]
  }
  return(scelist)
}

#' @keywords internal
.batch_sizes <- function(matlist){
  df <- as.data.frame(do.call(rbind, lapply(matlist,dim)))
  colnames(df) <- NULL
  rownames(df) <- names(matlist)
  return(df)
}

#' @keywords internal
#' @import ggthemes
#' @importFrom grDevices colorRampPalette
.generate_palette_func <- function(ncolors, color_values){
  # adapted from ggthemes::scale_color_* functions
  if(missing(color_values)){
    values <- ggthemes::ggthemes_data$few$colors[['Medium']][['value']]
  }
  else{values <- color_values}
  values <- c(values[1],colorRampPalette(values[seq(2,length(values),1)])(1 + ncolors)) # expanding palette to fit
  max_n <- length(values) - 1L
  f <- function(n) {
    ggthemes:::check_pal_n(n, max_n)
    if (n == 1L) {
      values[[1L]]
    }
    else {
      unname(values[2L:(n + 1L)])
    }
  }
  attr(f, "max_n") <- length(values) - 1L
  f
}


#' all_are
#'
#' Checks if all elements of a list or List are of a (single) particular type \code{typechar}
#'
#' @param inplist list of List to be checked
#' @param typechar char to check for
#'
#' @return boolean, for whether the elements of \code{inplist} are all \code{typechar}
#' @export
#'
#' @examples
#' x <- list(1,2)
#' all_are(x,'numeric')
#' all_are(x,'char')
#' 
#' y <- list(1,2,'c')
#' all_are(y,'numeric')
#' all_are(y,'char')
all_are <- function(inplist,typechar){
  return(all(unlist(lapply(inplist,is,typechar))))
}

#' rv coefficient
#'
#' @param mat1 matrix (or matrix-like, e.g., df); either columns or rows should be matched with \code{mat2}
#' @param mat2 matrix (or matrix-like, e.g., df); either columns or rows should be matched with \code{mat1}
#'
#' @return numeric; RV coefficient between the matched matrices
#' @export
#'
#' @examples
#' a <- matrix(sample(1:10,100, TRUE), nrow = 10)
#' b <- matrix(sample(1:10,50, TRUE), nrow = 5)
#' 
#' rv(a, b) # matched by columns
#' rv(t(a), t(b)) # matched by rows
rv <- function(mat1, mat2) {
  if(ncol(mat1) != ncol(mat2)){
    if(nrow(mat1) == nrow(mat2)){
      mat1 <- t(mat1)
      mat2 <- t(mat2)
    }
    else{stop('Matrices must be aligned by rows or columns.')}
  }
  nscm1 <- crossprod(as.matrix(mat1))
  nscm2 <- crossprod(as.matrix(mat2))
  rv <- sum(nscm1 * nscm2)/(sum(nscm1 * nscm1) * sum(nscm2 * nscm2))^0.5
  return(rv)
}


#' Pairwise rv coefficient
#'
#' @param matlist list of matrices (or matrix-like; see \code{rv} function) for which to compute pairwise RV coefficients
#'
#' @return matrix of the pairwise coefficients
#' @export
#' 
#' @importFrom utils combn
#'
#' @examples
#' a <- matrix(sample(1:10,100,TRUE), nrow = 10)
#' b <- matrix(sample(1:10,50,TRUE), nrow = 5)
#' c <- matrix(sample(1:10,20,TRUE), nrow = 2)
#' 
#' matlist <- list(a,b,c)
#' pairwise_rv(matlist)
#' pairwise_rv(lapply(matlist, t))
pairwise_rv <- function(matlist){
  n <- length(matlist)
  a <- utils::combn(seq_len(n), 2, 
             FUN = function(x) rv(matlist[[x[1]]], matlist[[x[2]]]), simplify = TRUE)
  m <- matrix(1, n, n)
  m[lower.tri(m)] <- a
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  colnames(m) <- rownames(m) <- names(matlist)
  return(m)
}



