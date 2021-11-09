#' Generate a dataframe of the scaled variance values
#'
#' @param inp corralm object or matrix; embedding to compute scaled variances
#' @param batchvec vector; batch labels (can be numeric or char). Defaults to `NULL`, which is appropriate for using a corralm object. If using an embedding matrix for inp, then this argument must be given and length must correspond to number of rows in `inp`.
#'
#' @return dataframe of the scaled variance values by PC (batches in columns)
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(5000), ncol = 50)
#' bv <- rep(seq(3),c(10,30,60))
#' scal_var_df(dat, bv)
scal_var_df <- function(inp, batchvec = NULL){
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
  sv_df <- Reduce(cbind,lapply(sep_var, t))
  sv_df[,1]
  
  overall_var <- apply(emb, MARGIN = 2, FUN = var)
  
  sv_df <- sweep(sv_df, MARGIN = 2, STATS = overall_var, FUN = '/')
  data.frame(sv_df)
}

#' Generate a scaled variance plot for an integrative embedding
#'
#' @param pcs numeric; vector of which PCs should be shown. Defaults to 1:3
#' @param plot_subtitle string; the text that should show in the subtitle for the plot. defaults to NULL
#' @inheritParams scal_var_df
#' @inheritParams plot_embedding
#'
#' @return N/A or a ggplot object
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(10000), ncol = 50)
#' bv <- rep(seq(4),c(10,30,60,100))
#' scal_var(dat,bv, pcs = seq(4))
scal_var <- function(inp, batchvec = NULL, pcs = seq(3), returngg = FALSE, showgg = TRUE, plot_subtitle = NULL){
  svdf <- scal_var_df(inp, batchvec)
  svdfmelt <- reshape2::melt(svdf[,pcs])
  colnames(svdfmelt) <- c('Batch','PC','sv')
  svdfmelt$Batch <- as.factor(svdfmelt$Batch)
  
  svdfmelt$PC_ind <- svdfmelt$PC + seq(length(svdfmelt$PC))
  
  pc_seplines <- seq(2,max(svdfmelt$PC_ind))[-which(seq(2,max(svdfmelt$PC_ind)) %in% svdfmelt$PC_ind)]
  pc_seplines <- c(pc_seplines, max(svdfmelt$PC_ind)+1)
  pc_lablocs <- pc_seplines - nrow(svdf)/2
  
  ggobj <- ggplot(svdfmelt, aes(x = PC_ind, y = sv, colour = Batch, group = Batch)) + 
    geom_segment(aes(x=PC_ind, xend=PC_ind, y=1, yend=sv), color="gray") +
    geom_point(size = 2) + ggthemes::scale_color_hc() + # other color option - scale_color_few()
    theme_classic() + 
    geom_hline(yintercept=1, color = 'gray') + 
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
  
  if(showgg){
    show(ggobj)
  }
  if(returngg){
    return(ggobj)
  }
}


