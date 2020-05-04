#' Plot selected PCs from an embedding
#'
#' @param embedding matrix or other tabular format where columns correspond to PCs and rows correspond to cells (entries)
#' @param xpc int; which PC to put on the x-axis (defaults to 1)
#' @param ypc int; which PC to put on the y-axis (defaults to the one after \code{xpc})
#' @param plot_title char; title of plot (defaults to titling based on \code{xpc} and \code{ypc})
#' @param color_vec vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, which will be colored the same. (e.g., this could be indicating which batch each cell is from)
#' @param color_title char; what attribute the colors represent
#' @param ellipse_vec OPTIONAL vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, and elements of the same class will be circled in an ellipse. (e.g., this could be indicating the cell type or cell line; works best for attributes intended to be compact)
#' @param saveplot boolean; whether or not to save the plot, defaults \code{TRUE}
#' @param plotfn char; what the filename is to be called. (defaults to making a name based on \code{plot_title} and \code{xpc})
#' @param showplot boolean; whether or not to show the plot, defaults \code{TRUE}
#' @param returngg boolean; whether or not to return a \code{\link{ggplot2}} object, defaults \code{FALSE}
#'
#' @return default none; options to display plot (\code{showplot}), save plot (\code{saveplot}), and/or return \code{\link{ggplot2}} object (\code{returngg})
#' @export
#' 
#' @import ggplot2
#' @importFrom ggthemes scale_color_few
#'
#' @examples
#' listofmats <- list(matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20),
#'                    matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20))
#' embed_mat <- corralm(listofmats, ncomp = 5)$v
#' cell_type_vec <- sample(c('type1','type2','type3'),100,replace = TRUE)
#' plot_embedding(embedding = embed_mat, 
#'                xpc = 1, 
#'                plot_title = 'corralm plot',
#'                color_vec = cell_type_vec, 
#'                color_title = 'cell type',
#'                saveplot = FALSE)
#' 
plot_embedding <- function(embedding, xpc = 1, ypc = xpc + 1, plot_title = paste0('PC',xpc,' by PC',ypc), color_vec, color_title, ellipse_vec = NULL, saveplot = TRUE, plotfn = paste(plot_title,xpc, sep = '_'), showplot = TRUE, returngg = FALSE){
  xvar <- as.name(paste('X',xpc,sep = ''))
  yvar <- as.name(paste('X',ypc,sep = ''))
  
  xlab = paste0('PC',xpc)
  ylab = paste0('PC',ypc)

  colnames(embedding) <- NULL
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
#' @param sce \code{\link{SingleCellExperiment}} object; contains the embedding within the \code{reducedDim} slot
#' @param which_embedding character; for the embedding to plot
#' @param color_attr character; name of the attribute within \code{colData} to use for assigning colors (in lieu of \code{color_vec} in the \code{\link{plot_embedding}} function)
#' @param color_title character; title to use for colors legend, defaults to the same as \code{color_attr}
#' @param ellipse_attr OPTIONAL character; name of the attribute within \code{colData} to use for drawing ellipse(s) (in lieu of \code{ellipse_vec} in the \code{\link{plot_embedding}} function)
#' @param ... additional optional arguments - see \code{\link{plot_embedding}} function for details on other potential arguments: \code{xpc}, \code{ypc}, \code{plot_title}, \code{color_title} (if title is different from \code{color_attr}), \code{saveplot}, \code{plotfn}, \code{showplot}, \code{returngg}
#' 
#' @return default none; options to display plot (\code{showplot}), save plot (\code{saveplot}), and/or return \code{\link{ggplot2}} object (\code{returngg})
#' @export
#' 
#' @import ggplot2
#' @importFrom ggthemes scale_color_few
#' @importFrom SingleCellExperiment colData reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' library(DuoClustering2018)
#' library(SingleCellExperiment)
#' sce <- sce_full_Zhengmix4eq()[1:100,sample(1:3500,100,replace = FALSE)]
#' colData(sce)$Method <- matrix(sample(c('Method1','Method2'),100,replace = TRUE))
#' sce <- corralm(sce, splitby = 'Method')
#' 
#' # to plot and show only
#' plot_embedding_sce(sce = sce, 
#'                    which_embedding = 'corralm', 
#'                    xpc = 1, 
#'                    plot_title = 'corralm: PC1 by PC2',
#'                    color_attr = "Method", 
#'                    ellipse_attr = 'phenoid',
#'                    saveplot = FALSE)
#' 
#' # to return ggplot2 object and display, but not save
#' corralm_ggplot <- plot_embedding_sce(sce = sce, 
#'                                      which_embedding = 'corralm', 
#'                                      xpc = 1, 
#'                                      plot_title = 'corralm: PC1 by PC2',
#'                                      color_attr = 'Method', 
#'                                      ellipse_attr = 'phenoid', 
#'                                      returngg = TRUE, 
#'                                      saveplot = FALSE)
#' 
#' # or, to plot UMAP, using 'scater' package:
#' library(scater)
#' result <- runUMAP(sce, 
#'                   dimred = 'corralm', 
#'                   name = 'corralm_UMAP')
#' plot_embedding_sce(sce = result, 
#'                    which_embedding = 'corralm_UMAP', 
#'                    xpc = 1, 
#'                    plot_title = 'corralm UMAP: PC1 by PC2',
#'                    color_attr = "Method", 
#'                    ellipse_attr = 'phenoid',
#'                    saveplot = FALSE)
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
