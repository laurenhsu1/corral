#' Plot selected PCs from an embedding
#'
#' @param embedding matrix or other tabular format where columns correspond to PCs and rows correspond to cells (entries). \code{corral} and \code{corralm}  objects are also accepted.
#' @param xpc int; which PC to put on the x-axis (defaults to 1)
#' @param ypc int; which PC to put on the y-axis (defaults to the one after \code{xpc})
#' @param plot_title char; title of plot (defaults to titling based on \code{xpc} and \code{ypc})
#' @param color_vec vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, which will be colored the same. (e.g., this could be indicating which batch each cell is from)
#' @param color_title char; what attribute the colors represent
#' @param ellipse_vec vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class, and elements of the same class will be circled in an ellipse. (e.g., this could be indicating the cell type or cell line; works best for attributes intended to be compact)
#' @param facet_vec vector; length should correspond to the number of rows in embedding, and each element of the vector classifies that cell (entry) in the embedding to that particular class. Plot will be faceted by this attribute.
#' @param ptsize numeric; the size of the points as passed to \code{geom_point()}. Defaults to 0.8.
#' @param saveplot boolean; whether or not to save the plot, defaults \code{FALSE}
#' @param plotfn char; what the filename is to be called. (defaults to making a name based on \code{plot_title} and \code{xpc})
#' @param showplot boolean; whether or not to show the plot, defaults \code{TRUE}
#' @param returngg boolean; whether or not to return a \code{\link{ggplot2}} object, defaults \code{FALSE}
#' @param color_pal_vec char; hex codes for the color palette to be used. Default is to use the ggthemes few for plots with less than 9 colors, and to use/"stretch" pals polychrome if more colors are needed.
#' @param dimname char; the name of the dimensions. defaults to "Dim"
#'
#' @return default none; options to display plot (\code{showplot}), save plot (\code{saveplot}), and/or return \code{\link{ggplot2}} object (\code{returngg})
#' @export
#' 
#' @import ggplot2
#' @import gridExtra
#' @import ggthemes
#' @import pals
#'
#' @examples
#' listofmats <- list(matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20),
#'                    matrix(sample(seq(0,20,1),1000,replace = TRUE),nrow = 20))
#' corralm_obj <- corralm(listofmats, ncomp = 5)
#' embed_mat <- corralm_obj$v
#' cell_type_vec <- sample(c('type1','type2','type3'),100,replace = TRUE)
#' plot_embedding(embedding = embed_mat, 
#'                xpc = 1, 
#'                plot_title = 'corralm plot',
#'                color_vec = cell_type_vec, 
#'                color_title = 'cell type',
#'                saveplot = FALSE)
#' 
#' # or, call directly on the corralm object              
#' plot_embedding(corralm_obj)
#' 
plot_embedding <- function(embedding, xpc = 1, ypc = xpc + 1, plot_title = paste0('Dim',xpc,' by Dim',ypc), color_vec = NULL, color_title = NULL, ellipse_vec = NULL, facet_vec = NULL, ptsize = 0.8, saveplot = FALSE, plotfn = paste(plot_title,xpc, sep = '_'), showplot = TRUE, returngg = FALSE, color_pal_vec = NULL, dimname = 'Dim'){
  if('corral' %in% class(embedding)){
    embedding <- embedding$PCv
  }
  if('corralm' %in% class(embedding)){
    corralm_obj <- embedding
    embedding <- corralm_obj$v
    if(is.null(color_vec)){
      color_vec <- rep(rownames(corralm_obj$batch_sizes), corralm_obj$batch_sizes[,2])
      color_title <- 'Batch'
    }
  }
  if(is.null(color_vec)){
    color_vec <- rep('emb', nrow(embedding))
  }
  
  xlab <- paste0(dimname, xpc)
  ylab <- paste0(dimname, ypc)

  df <- cbind(data.frame(embedding[,c(xpc, ypc)]), color_vec)
  
  colnames(df) <- c('Xdim','Ydim', 'color_vec')
  
  if(!is.null(ellipse_vec)){
    df <- cbind(df, ellipse_vec)
  }
  
  if(!is.null(facet_vec)){
    df <- cbind(df, facet_vec)
  }
  
  # Setting up the colors
  .colscale <- function(palette){ggplot2::discrete_scale('colour','colscale',palette)}
  
  if(!is.null(color_pal_vec)){
    palette_func <- .generate_palette_func(ncolors = length(unique(color_vec)), color_pal_vec)
  } else if(length(unique(color_vec)) < 9){
    palette_func <- ggthemes::few_pal('Medium')
  } else{
    palette_func <- .generate_palette_func(ncolors = length(unique(color_vec)), pals::alphabet2())
  }
  
  gg_obj <- ggplot(df, aes(x = Xdim, y = Ydim, colour = color_vec)) + 
    geom_point(size = ptsize) +
    theme_classic() + .colscale(palette_func) + 
    geom_hline(yintercept = 0, color = 'gray') + 
    geom_vline(xintercept = 0, color = 'gray') + 
    labs(x = xlab, y = ylab, 
         title = plot_title, 
         color = color_title) + 
    theme(axis.text.x = element_text(size = rel(1.4), colour = 'black'), 
          axis.text.y = element_text(size = rel(1.4), colour = 'black'), 
          axis.title = element_text(size = rel(1.4), colour = 'black')) + 
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  if(!is.null(ellipse_vec)){
    gg_obj <- gg_obj + stat_ellipse(aes(x = Xdim, y = Ydim, group = ellipse_vec), type = 'norm', linetype = 2)
  }
  
  if(!is.null(facet_vec)){
    gg_obj <- gg_obj + facet_wrap(facets = vars(facet_vec))
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
#' @param ellipse_attr character; name of the attribute within \code{colData} to use for drawing ellipse(s) (in lieu of \code{ellipse_vec} in the \code{\link{plot_embedding}} function)
#' @param facet_attr character; name of the attribute within \code{colData} to use for faceting (in lieu of \code{facet_vec} in the \code{\link{plot_embedding}} function)
#' @param ... additional optional arguments - see \code{\link{plot_embedding}} function for details on other potential arguments: \code{xpc}, \code{ypc}, \code{plot_title}, \code{color_title} (if title is different from \code{color_attr}), \code{ptsize}, \code{saveplot}, \code{plotfn}, \code{showplot}, \code{returngg}, \code{color_pal_vec}, \code{dimname}
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
plot_embedding_sce <- function(sce, which_embedding, color_attr = NULL, color_title = color_attr, ellipse_attr = NULL, facet_attr = NULL, ...){
  embed_mat <- SingleCellExperiment::reducedDim(sce, which_embedding)
  
  if(is.null(color_attr)){
    color_vec <- NULL
  }
  else{
    color_vec <- SingleCellExperiment::colData(sce)[, color_attr]
  }
  
  if(!is.null(ellipse_attr)){
    ellipse_vec <- SingleCellExperiment::colData(sce)[, ellipse_attr]
  }
  else{
    ellipse_vec <- NULL
  }
  
  if(!is.null(facet_attr)){
    facet_vec <- SingleCellExperiment::colData(sce)[, facet_attr]
  }
  else{
    facet_vec <- NULL
  }
    
  plot_embedding(embed_mat, color_vec = color_vec, color_title = color_title, ellipse_vec = ellipse_vec, facet_vec = facet_vec, ...)
}
