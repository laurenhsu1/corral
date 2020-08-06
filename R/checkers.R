#' @keywords internal
.check_mat <- function(mat){
  if(!is.numeric(mat)){
    stop('Matrix input should contain only numbers.')
  }
  else if (min(mat) < 0){
    stop('Matrix input should only contain positive values.')
  }
  else if (min(dim(mat)) < 2){
    stop('Matrix input must have more than 1 element and more than 1 feature.')
  }
}

#' @keywords internal
.check_dims <- function(matlist){
  dims <- unlist(lapply(matlist,nrow))
  if(length(unique(dims)) > 1) {
    stop('If performing multi-table analysis, the matrices must be matched by rows; currently the dimensions do not match. \nIf they are matched by columns, then transpose the matrices.')}
}

#' @keywords internal
.check_rw_contrib <- function(matlist, rw_contrib){
  # Verifying valid input for rw_contrib, otherwise setting to equal weight.
  matlist_len <- length(matlist)
  if(matlist_len != length(rw_contrib)){
    cat('\nThe provided weights did not match number of batches (i.e., number of matrices).\nThey will be set to equal weight.\n')
    return(rep(1, matlist_len))
  }
  else if(sum(is.na(as.numeric(rw_contrib)))){
    cat('\nNon-numeric values provided in rw_contrib.\nThey will be set to equal weight.\n')
    return(rep(1, matlist_len))
  }
  else{
    return(rw_contrib)
  }
}



