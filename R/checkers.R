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

.check_dims <- function(matlist){
  dims <- unlist(lapply(lapply(matlist,dim),'[',1))
  if(length(unique(dims)) > 1) {
    stop('If performing multi-table analysis, the matrices must be matched by rows; currently the dimensions do not match. \nIf they are matched by columns, then transpose the matrices.')}
}
