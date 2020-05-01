# corral
*Dimensionality reduction and batch integration methods for single cell data*

corral is a package to perform correspondence analysis on count data, particularly for single-cell (e.g., RNA-seq).

It can be used on a single table for dimensionality reduction, or it can be used for integration across batches, samples, and modes.

To install the dev version: \code{devtools::install_github('laurenhsu1/corral')}


# How to

The package contains two main function calls:

1. \code{corral} for dimensionality reduction on a single table 
2. \code{corralm} for integration across tables.

See the vignettes for details on outputs and how to plot.
