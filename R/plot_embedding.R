#' Plot selected PCs from an embedding
#'
#' @param embedding matrix or other tabular format where columns correspond to PCs and rows correspond to cells (entries)
#' @param xpc int; which PC to put on the x-axis (defaults to 1)
#' @param ypc int; which PC to put on the y-axis (defaults to the one after xpc)
#' @param plot_title char; title of plot (defaults to titling based on xpc and ypc)
#' @param color_vec vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, which will be colored the same. (e.g., this could be indicating which batch each cell is from)
#' @param color_title char; what attribute the colors represent
#' @param ellipse_vec OPTIONAL vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, and elements of the same class will be circled in an ellipse. (e.g., this could be indicating the cell type or cell line; works best for attributes intended to be compact)
#' @param saveplot boolean; whether or not to save the plot, defaults TRUE
#' @param plotfn char; what the filename is to be called. (defaults to making a name based on plot title and xpc)
#' @param showplot boolean; whether or not to show the plot, defaults TRUE
#' @param returngg boolean; whether or not to return a ggplot2 object, defaults FALSE
#'
#' @return default none; options to display plot (showplot), save plot (saveplot), and/or return ggplot2 object (returngg)
#' @export
#'
#' @examples
#' 
#' embed_mat <- corralm(list_of_mats)
#' plot_embedding(embedding = embed_mat, xpc = 1, plot_title = 'corralm plot',color_vec = cell_type_vec, color_title = 'cell type')
#' 
plot_embedding <- function(embedding, xpc = 1, ypc = xpc + 1, plot_title = paste0('PC',xpc,' by PC',ypc), color_vec, color_title, ellipse_vec = NULL, saveplot = TRUE, plotfn = paste(plot_title,xpc, sep = '_'), showplot = TRUE, returngg = FALSE){
  xvar <- as.name(paste('X',xpc,sep = ''))
  yvar <- as.name(paste('X',ypc,sep = ''))
  
  xlab = paste0('PC',xpc)
  ylab = paste0('PC',ypc)
  
  # alternate axis label, given the actual svd output with d vector. If partial svd, then also need the input matrix to compute
  # sum of eigenvalues
  # pct_var_exp <- get_pct_var_exp_svd(thissvd)
  #xlab = paste('PC ',xpc,' (', round(pct_var_exp[xpc]*100),'% of variance explained)', sep = '' )
  #ylab = paste('PC ',xpc + 1,' (', round(pct_var_exp[xpc + 1]*100),'% of variance explained)', sep = '' )
  
  df <- cbind(data.frame(embedding),color_vec)
  colvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  
  if(!is.null(ellipse_vec)){
    df <- cbind(df,ellipse_vec)
    ellvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  }
  
  gg_obj <- ggplot(df, aes(x = eval(xvar), y = eval(yvar), colour = eval(colvar))) + theme_classic() + 
    geom_hline(yintercept=0, color = 'gray') + geom_vline(xintercept=0, color = 'gray') + 
    geom_point() + scale_color_few() + 
    labs(x = xlab, 
         y = ylab, 
         title = plot_title, 
         color = color_title) + 
    theme(axis.text.x = element_text(size = rel(1.4), colour = 'black'), axis.text.y = element_text(size = rel(1.4), colour = 'black'), axis.title = element_text(size = rel(1.4), colour = 'black'))
  
  if(!is.null(ellipse_vec)){
    gg_obj <- gg_obj + stat_ellipse(aes(x = eval(xvar), y = eval(yvar), group = eval(ellvar)), type = 'norm', linetype = 2)
  }
  
  if(saveplot){
    ggsave(paste0(plotfn,'.png'))
  }
  
  if(showplot){
    plot(gg_obj)
  }
  
  if(returngg){
    return(gg_obj)
  }
}


#' Plot selected PCs from an embedding saved in a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object; contains the embedding within the 'reducedDim' slot
#' @param which_embedding character; for the embedding to plot
#' @param color_attr character; name of the attribute within colData to use for assigning colors (in lieu of "color_vec" in the plot_embedding function)
#' @param color_title character; title to use for colors legend, defaults to the same as color_attr
#' @param ellipse_attr OPTIONAL character; name of the attribute within colData to use for drawing ellipse(s) (in lieu of "ellipse_vec" in the plot_embedding function)
#' @param ... additional optional arguments - see plot_embedding function for details on other potential arguments: xpc, ypc, plot_title, color_title (if title is different from color_attr), saveplot, plotfn, showplot, returngg
#' 
#' @return default none; options to display plot (showplot), save plot (saveplot), and/or return ggplot2 object (returngg)
#' @export
#'
#' @examples
#' 
#' sce <- corralm(sce, 'Method')
#' 
#' # to plot and show only
#' plot_embedding_sce(sce = sce, which_embedding = 'corralm', xpc = 1, plot_title = 'corralm: PC1 by PC2',color_attr = "Method", ellipse_attr = 'cell_type')
#' 
#' # to return ggplot2 object and display, but not save
#' corralm_ggplot <- plot_embedding_sce(sce = sce, which_embedding = 'corralm', xpc = 1, plot_title = 'corralm: PC1 by PC2',color_attr = 'Method', ellipse_attr = 'cell_type', returngg = TRUE, saveplot = FALSE)
#' 
#' # or, to plot UMAP, using 'scater' package:
#' library(scater)
#' result <- runUMAP(sce, dimred = 'corralm', name = 'corralm_UMAP')
#' plot_embedding_sce(sce = result, which_embedding = 'corralm_UMAP', xpc = 1, plot_title = 'corralm UMAP: PC1 by PC2',color_attr = "Method", ellipse_attr = 'cell_type')
#' 
plot_embedding_sce <- function(sce, which_embedding, color_attr, color_title = color_attr, ellipse_attr = NULL, ...){
  embed_mat <- SingleCellExperiment::reducedDim(sce, which_embedding)
  
  color_vec <- SingleCellExperiment::colData(sce)[, color_attr]
  
  if(!is.null(ellipse_attr)){
    ellipse_vec <- SingleCellExperiment::colData(sce)[, ellipse_attr]
  }
  else{
    ellipse_vec <- NULL
  }
    
  plot_embedding(embed_mat, color_vec = color_vec, color_title = color_title, ellipse_vec = ellipse_vec, ...)
}