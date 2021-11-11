#' Generate a matrix of the scaled variance values
#'
#' @param inp corralm object or matrix; embedding to compute scaled variances
#' @param batchvec vector; batch labels (can be numeric or char). Defaults to `NULL`, which is appropriate for using a corralm object. If using an embedding matrix for inp, then this argument must be given and length must correspond to number of rows in `inp`.
#'
#' @return matrix of the scaled variance values by PC (batches in rows; PCs in columns)
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(5000), ncol = 50)
#' bv <- rep(seq(3),c(10,30,60))
#' scal_var_mat(dat, bv)
scal_var_mat <- function(inp, batchvec = NULL){
  if('corralm' %in% class(inp)){
    splitvec <- rep(rownames(inp$batch_sizes),inp$batch_sizes[,2])
    emb <- inp$v
  }
  else if(!is.null(batchvec)){
    if(length(batchvec) != nrow(inp)){
      stop('batchvec must be same length as number of rows in inp')
    }
    emb <- inp
    splitvec <- batchvec
  }
  else{
    stop('Please provide a valid input: either corralm object or embedding with batchvec arg')
  }
  sepemb <- apply(emb, MARGIN = 2, FUN = split, f = splitvec)
  sep_var <- lapply(sepemb, FUN = function(x) rbind(unlist(lapply(x, var))))
  sv_mat <- Reduce(cbind,lapply(sep_var, t))
  sv_mat[,1]
  
  overall_var <- apply(emb, MARGIN = 2, FUN = var)
  
  sv_mat <- sweep(sv_mat, MARGIN = 2, STATS = overall_var, FUN = '/')
  sv_mat
}

#' Generate a scaled variance plot for an integrative embedding
#'
#' @param pcs numeric; vector of which PCs should be shown. Defaults to 1:3
#' @param plot_subtitle string; the text that should show in the subtitle for the plot. defaults to NULL
#' @inheritParams scal_var_mat
#' @inheritParams plot_embedding
#'
#' @import ggplot2
#' @import ggthemes
#' @import reshape2
#' @return N/A or a ggplot object
#' @export
#' 
#' @importFrom stats quantile var
#'
#' @examples
#' dat <- matrix(rnorm(10000), ncol = 50)
#' bv <- rep(seq(4),c(10,30,60,100))
#' scal_var(dat,bv, pcs = seq(4))
scal_var <- function(inp, batchvec = NULL, pcs = seq(3), returngg = FALSE, showplot = TRUE, plot_subtitle = NULL){
  svmat <- scal_var_mat(inp, batchvec)
  svmelt <- reshape2::melt(svmat[,pcs])
  colnames(svmelt) <- c('Batch','PC','sv')
  svmelt$Batch <- as.factor(svmelt$Batch)
  
  svmelt$PC_ind <- svmelt$PC + seq(length(svmelt$PC))
  
  pc_seplines <- seq(2,max(svmelt$PC_ind))[-which(seq(2,max(svmelt$PC_ind)) %in% svmelt$PC_ind)]
  pc_seplines <- c(pc_seplines, max(svmelt$PC_ind)+1)
  pc_lablocs <- pc_seplines - nrow(svmat)/2
  
  ggobj <- ggplot(svmelt, aes(x = PC_ind, y = sv, colour = Batch, group = Batch)) + 
    geom_hline(yintercept=1, color = 'gray') + 
    geom_segment(aes(x=PC_ind, xend=PC_ind, y=1, yend=sv), color="gray") +
    geom_point(size = 2) + ggthemes::scale_color_hc() + # other color option - scale_color_few()
    theme_classic() + 
    geom_vline(xintercept = pc_seplines, color = '#555555') +
    labs(x = 'Component',
         y = 'Scaled variance by group for top PCs',
         title = 'Scaled variance by batch',
         subtitle = plot_subtitle,
         colour = 'Batch') + 
    scale_x_continuous(breaks = pc_lablocs, 
                       labels = paste0('PC',pcs)) +
    theme(axis.ticks.x = element_blank(), axis.title.x = element_text(size = rel(1.2))) +
    ylim(0,NA)
  
  if(showplot){
    show(ggobj)
  }
  if(returngg){
    return(ggobj)
  }
}


