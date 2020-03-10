# Functions for evaluating performance of the methods

library(LaplacesDemon)
library(pryr)
library(transport)


#' Evaluate performance
#'
#' @param routine code to be evaluated
#'
#' @return list with 'time' and 'memory change' as elements
#' @export
#'
#' @examples
#' perf_report <- eval_perf(my_corral <- corral(my_matrix))
eval_perf <- function(routine){
  mem_ch <- mem_change(timing <- system.time(routine))
  perf_report <- list()
  perf_report[['time']] <- timing
  perf_report[['memory change']] <- mem_ch
  return(perf_report)
}

cumtotal <- function(vals, ref){
  return(sum(vals < ref))
}

obs2probs <- function(obs, numbins = 100, startbin = min(obs), endbin = max(obs) + .00001){
  bins <- seq(from = startbin, to = endbin, length.out = numbins)
  result <- data.frame(bins)
  result[1,'cumfreq'] <- 0
  result[1, 'freq'] <- 0
  for (ind in 2:(numbins)){
    cumsum <- cumtotal(obs, bins[ind])
    result[ind, 'cumfreq'] <- cumsum
    result[ind, 'freq'] <- cumsum - result[ind - 1, 'cumfreq']
  }
  result$probs <- result$freq / length(obs)
  return(result)
}

make_costmat <- function(matdim, mincost = 1, maxcost = 10){
  abval_dif <- function(x,y) {return(abs(x-y))}
  costmat <- matrix(seq(mincost, maxcost,length.out = matdim), matdim, matdim)
  costmat <- sweep(costmat, 
                   MARGIN = 2,
                   STATS = seq(mincost, maxcost, length.out = matdim),
                   FUN = abval_dif)
  return(costmat)
}


earthmover_dist <- function(batch1, batch2, whichdim = 1, numbins = 100, p_param = 1){
  # i.e., wasserstein distance with L1 (p_param = 1); can also use other penalties > 1
  # (Not technically earthmover distance if using other p_param values)
  minval <- min(min(batch1), min(batch2))
  maxval <- max(max(batch1), max(batch2)) + .00001
  df_1 <- obs2probs(obs = batch1[,whichdim], numbins = numbins, startbin = minval, endbin = maxval)
  df_2 <- obs2probs(obs = batch2[,whichdim], numbins = numbins, startbin = minval, endbin = maxval)
  costmat <- make_costmat(matdim = numbins)
  wasserstein(df_1$probs, df_2$probs, costm = costmat, p = p_param)
}

# write the KLD wrapper



# bring in the var ratio method from frontiers paper




