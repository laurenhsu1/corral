# Functions for evaluating performance with respect to software / computing, as well as results for dimensionality reduction and batch integration.

#' @keywords internal
.cumtotal <- function(vals, ref){
  return(sum(vals < ref))
}

#' Observations --> discrete probabilities
#' 
#' usage:
#' embedding <- matrix(sample(x = seq(0,10,.1),200, replace = TRUE))
#' disc_probs <- obs2probs(embedding)
#'
#' @param obs vector of numeric, with the observations
#' @param numbins int, the number of evenly sized bins to discretize the observations to
#' @param startbin numeric, the starting value for the smallest bin. Defaults to taking the minimum of obs
#' @param endbin numeric, the ending value for the largest bin. Defaults to taking the maximum of obs (plus a tiny decimal to ensure full range of obs is captured)
#'
#' @return dataframe, results has rows corresponding to each bin with columns for probability ('prob'), cumulative frequency ('cumfreq'), and frequency ('freq') of observations falling into that bin. The 'bins' column indicates the end of the bin (start is the preceding column)
#' @keywords internal
#'
obs2probs <- function(obs, numbins = 100, startbin = min(obs), endbin = max(obs) + .00001){
  bins <- seq(from = startbin, to = endbin, length.out = numbins)
  result <- data.frame(bins)
  result[1,'cumfreq'] <- 0
  result[1, 'freq'] <- 0
  for (ind in seq(2,numbins,1)){
    cumsum <- .cumtotal(obs, bins[ind])
    result[ind, 'cumfreq'] <- cumsum
    result[ind, 'freq'] <- cumsum - result[ind - 1, 'cumfreq']
  }
  result$probs <- result$freq / length(obs)
  return(result)
}

#' @keywords internal
.make_costmat <- function(matdim, mincost = 1, maxcost = 10){
  abval_dif <- function(x,y) {return(abs(x-y))}
  costmat <- matrix(seq(mincost, maxcost,length.out = matdim), matdim, matdim)
  costmat <- sweep(costmat, 
                   MARGIN = 2,
                   STATS = seq(mincost, maxcost, length.out = matdim),
                   FUN = abval_dif)
  return(costmat)
}


#' Earthmover distance
#' i.e., wasserstein distance with L1 (p_param = 1); can also use other penalties > 1
#' (Not technically earthmover distance if using other p_param values)
#'
#' @param batch1 matrix; subset of observations from an embedding correponding to some attribute (e.g., batch or phenotype)
#' @param batch2 matrix; subset of observations from an embedding correponding to some attribute (e.g., batch or phenotype)
#' @param whichdim int; which dimension (i.e., column) from the embeddings is used. defaults on first
#' @param numbins int; number of bins for the probability discretization (defaults to 100)
#' @param p_param int; penalty parameter for wasserstein distance. Defaults to 1, which corresonds to earthmover.
#'
#' @return num; the distance
#' @export
#' 
#' @importFrom transport wasserstein
#'
#' @examples
#' # To compare distributions of reduced dimension values to assess similarity, 
#' # e.g. as a metric for batch integration
#' embedding <- matrix(sample(x = seq(0,10,.1),1000, replace = TRUE),ncol = 5)
#' batch <- matrix(sample(c(1,2),200, replace = TRUE))
#' earthmover_dist(embedding[which(batch == 1),],embedding[which(batch == 2),]) 
earthmover_dist <- function(batch1, batch2, whichdim = 1, numbins = 100, p_param = 1){
  minval <- min(min(batch1), min(batch2))
  maxval <- max(max(batch1), max(batch2)) + .00001
  df_1 <- obs2probs(obs = batch1[,whichdim], numbins = numbins, startbin = minval, endbin = maxval)
  df_2 <- obs2probs(obs = batch2[,whichdim], numbins = numbins, startbin = minval, endbin = maxval)
  costmat <- .make_costmat(matdim = numbins)
  transport::wasserstein(df_1$probs, df_2$probs, costm = costmat, p = p_param)
}
