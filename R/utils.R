
#' Set na to 0
#'
#' @param x matrix of values for which na values should be changed to 0
#'
#' @return matrix, where na values are set to 0
#' @export
#'
#' @examples
#' x <- matrix(sample(0:10, 5000, replace=TRUE), ncol = 25)
#' x[sample(1:5000, 10)] <- NA
#' apply(x, c(1,2), na2zero)
na2zero <- function(x) {
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
#' listofmats <- list(matrix(sample(seq(0,20,1),100,replace = TRUE),nrow = 10),matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 10))
#' newmat <- list2mat(listofmats) # to "cbind" them
#' listofmats_t <- lapply(listofmats,t)
#' newmat_t <- list2mat(listofmats_t, 'r') # to "rbind" them
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
#' matlist <- sce2matlist(sce = sce, splitby = 'phenoid', whichmat = logcounts)
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
#' library(DuoClustering2018)
#' sce <- sce_full_Zhengmix4eq()
#' scelist <- list(sce,sce)
#' embeddings <- matrix(sample(seq(0,20,1),dim(sce)[2]*6,replace = TRUE),nrow = dim(sce)[2]*2)
#' scelist <- add_embeddings2scelist(scelist, embeddings)
add_embeddings2scelist <- function(scelist,embeddings, slotname = 'corralm'){
  dimvec <- unlist(lapply(lapply(scelist,dim),'[',2))
  for (i in seq(1,length(scelist),1)){
    if(i == 1) {start_ind <- 1}
    else{start_ind <- (sum(dimvec[1:i-1]) + 1)}
    
    end_ind <- start_ind + dimvec[i] - 1
    
    SingleCellExperiment::reducedDim(scelist[[i]],slotname) <- embeddings[start_ind:end_ind,]
  }
  return(scelist)
}

.batch_sizes <- function(matlist){
  df <- as.data.frame(do.call(rbind, lapply(matlist,dim)))
  colnames(df) <- NULL
  rownames(df) <- names(matlist)
  return(df)
}
