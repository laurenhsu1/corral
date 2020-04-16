
#' Fix NA
#'
#' @param x matrix of values for which na values should be changed to 0
#'
#' @return fixed matrix, where na values are set to 0
#' @export
#'
#' @examples
#' fix_na(x)
fix_na <- function(x) {
  if (is.na(x)) 
    return(0)
  else return(x)
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
#' list2mat(my_matrix_list)
#' list2mat(my_matrix_list, 'r')
list2mat <- function(matlist, direction = c('c','r')[1]){
  concatted <- matlist[[1]]
  for(i in 2:length(matlist)){
    if (direction == 'r') concatted <- rbind(concatted,matlist[[i]])
    if (direction == 'c') concatted <- cbind(concatted,matlist[[i]])
  }
  return(unlist(concatted))
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
#' get_pct_var_exp_svd(my_svd) # this works if my_svd is a full svd
#' get_pct_var_exp_svd(my_svd, preproc_mat = my_mat) # ... otherwise use this 
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
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' get_w(my_mat)
#' lapply(my_matlist, get_w)
get_w <- function(inp_mat){
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
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' mymatlist <- sce2matlist(sce = mysce, splitby = 'Method', whichmat = logcounts)
sce2matlist <- function(sce,splitby,to_include = NULL,whichmat = 'counts'){
  if(is.null(to_include)){
    to_include <- unique(as.character(colData(sce)[,splitby]))
  }
  matlist <- list()
  countmatrix <- assay(sce, whichmat)
  for(gn in to_include){
    matlist[[gn]] <- countmatrix[,which(colData(sce)[,splitby] == gn)]
  }
  return(matlist)
}

#' Add embeddings to SCE
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
add_embeddings2scelist <- function(scelist,embeddings,slotname = 'corral'){
  dimvec <- unlist(lapply(lapply(scelist,dim),'[',2))
  for (i in seq(1,length(scelist),1)){
    if(i == 1) {start_ind <- 1}
    else{start_ind <- (sum(dimvec[1:i-1]) + 1)}
    
    end_ind <- start_ind + dimvec[i] - 1
    
    SingleCellExperiment::reducedDim(scelist[[i]],slotname) <- embeddings[start_ind:end_ind,]
  }
  return(scelist)
}
