# CompositionnalAnalysisR
R functions for analysis and projection of compositionnal data I develloped for my work.

## Requirements
Requires the following R packages : ggplot2, plotly and ggfortify. 

## How to use
Load the the functions file in your R project with source("path/to/compositionnal_fonctions.R"), then call directly the functions from you R project !

The functions are briefly presented here, the arguments for each are explained in the main file. Feel free to contact me for more informations !


### multiplicative_replacement
Function to replace missing values or rounded zeros, as described in Martín-Fernández, J. A., Barceló-Vidal, C., and Pawlowsky-Glahn, V., 2003, Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation, Mathematical Geology, 35(3), 253–78.

### normalise

Function to close a compositionnal dataset to a specified value.

### ALR_biplot

Function to plot compositionnal data, using additive log-ratio transformations, as described by Aitchison, J., 1982, The Statistical Analysis of Compositional Data, Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139–77. This function takes up to three datasets, and plots them in a single interactive plotly graph, wich can then be exported as html and used outside of R !

### CLR_biplot

Function to plot compositionnal data, using centered log-ratio transformations, as described by Aitchison, J., 1982, The Statistical Analysis of Compositional Data, Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139–77. This function takes up to three datasets, and plots them in a single interactive plotly graph, wich can then be exported as html and used outside of R !

### PCA_plot

Function to plot compositionnal data, using principal componant analysis (prcomp method). This function takes up to three datasets, and plots them in a single interactive plotly graph, wich can then be exported as html and used outside of R !

### PCA_plot_3D
Experimental function, identical to PCA_plot but create and interractive plot in 3 dimensions.
