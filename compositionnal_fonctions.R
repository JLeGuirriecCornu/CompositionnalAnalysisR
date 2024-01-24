library("ggplot2")
library("plotly")
library("ggfortify")


### Function for the analysis of compositional data, by Julien Le Guirriec-Cornu

# How to use : load the the functions in your R project with source("path/to/this/file")
# Call directly the functions from you R project, with appropriate arguments.


#data <- read.csv("simulation_standards", row.names = 1)

### Multiplicative replacement of missing values
# Method to replace missing values or rounded zeros, as described in Martín-Fernández, J. A., Barceló-Vidal, C., and Pawlowsky-Glahn, V., 2003, Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation, Mathematical Geology, 35(3), 253–78.

## Inputs :
# df : a data frame containing compositional data
# index (optionnal): a vector containing the index of columns of the data, if none supplied the function will assume the whole dataframe is compositionnal data to transform
# factor (default = 0.65) : 

multiplicative_replacement <- function(df, index, factor = 0.65){
  replace <- function(data, factor){
    minimumdetected <- apply(data[],2,min, na.rm=TRUE) * factor #Facteur arbitraitre, valeur par default = 0.65
    
    for (i in 1:nrow(data)){
      sum_values <- sum(data[i,], na.rm = TRUE)
      sum_missing <- sum(minimumdetected[which(is.na(data[i,]))])
      normfactor <- ((1000000 - sum_missing)/(sum_values))
      
      for (k in 1:ncol(data)){
        if (is.na(data[i,k])) {
          data[i,k] <- minimumdetected[k]
        } else{
          data[i,k] <- data[i,k] * normfactor
        }
      }
    }
    return(data)
  }
  if(missing(index)) {
    df <- replace(df, factor)
    return(df)
  }  else{
    df[,index] <- replace(df[,index], factor)
    return(df)
  }
  
}


### Normalisation

# Function to close a compositionnal dataset to a specified value
# data : a data frame containing only compositional data
# closure : closure value (numeric)

normalise <- function(data, closure) {
  for (i in 1:nrow(data)){
    sum_values <- sum(data[i,])
    normfactor <- (closure  / sum_values)
    
    for (k in 1:ncol(data)){
      data[i,k] <- data[i,k] * normfactor
    }
  }
  return(data)
}

### Additive log-ratio biplot
# Method to plot compositionnal data. 

## Inputs :

#Mandatory arguments for the reference dataset
# ref_data : dataframe containing compositional data only 
# ref_data_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# ref_data_group : vector defining the group each sample belongs to, corresponds to geological sources for example
# ref_data_name : character string defining the dataset, for the plot's title

#Mandatory arguments for the projection
# var1 : character string containing the name of the x variable
# var2 : character string containing the name of the y variable
# normalvar : character string containing the name of the normalisation variable

#Optionnal arguments for the first projected dataset (all must be supplied for the dataset to be projected)
# proj_data_1 : dataframe containing compositional data only 
# proj_data_1_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# proj_data_1_name : character string defining the dataset, for the plot's title and legend

#Optionnal arguments for the second projected dataset (all must be supplied for the dataset to be projected)
# proj_data_2 : dataframe containing compositional data only 
# proj_data_2_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# proj_data_2_name : character string defining the dataset, for the plot's title and legend

ALR_biplot <- function(ref_data, ref_data_id, ref_data_group, ref_data_name, var1, var2, normalvar, proj_data_1, proj_data_1_id, proj_data_1_name, proj_data_2, proj_data_2_id, proj_data_2_name) {
  
  ref_data <- data.frame(log(ref_data[,which(grepl(var1, colnames(ref_data)))]/ref_data[,which(grepl(normalvar, colnames(ref_data)))]),
                         log(ref_data[,which(grepl(var2, colnames(ref_data)))]/ref_data[,which(grepl(normalvar, colnames(ref_data)))]),
                         ref_data_id,
                         ref_data_group
  )
  
  colnames(ref_data) <- c("V1",
                           "V2",
                           "Info",
                           "Group"
  )
  
  
  if(missing(proj_data_1) | missing(proj_data_1_name) | missing(proj_data_1_id)){
    ref_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
      xlab(paste0('log(',var1,'/',normalvar,')')) +
      ylab(paste0('log(',var2,'/',normalvar,')')) + 
      stat_ellipse() + ggtitle(paste("ALR projection of ", ref_data_name))
    
    ggplotly(ref_plot)  
  } else {
    proj_data_1 <- data.frame(log(proj_data_1[,which(grepl(var1, colnames(proj_data_1)))]/proj_data_1[,which(grepl(normalvar, colnames(proj_data_1)))]),
                           log(proj_data_1[,which(grepl(var2, colnames(proj_data_1)))]/proj_data_1[,which(grepl(normalvar, colnames(proj_data_1)))]),
                           proj_data_1_id,
                           proj_data_1_name
    )
    
    colnames(proj_data_1) <- c("V1",
                            "V2",
                            "Info",
                            "Group"
    )
    
    if(missing(proj_data_2) | missing(proj_data_2_name) | missing(proj_data_2_id)){
      proj1_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
        xlab(paste0('log(',var1,'/',normalvar,')')) +
        ylab(paste0('log(',var2,'/',normalvar,')')) + 
        stat_ellipse() +  geom_point(data = proj_data_1, aes(V1, V2), shape = 18, size = 1) + ggtitle(paste("ALR projection of ", ref_data_name, "and", proj_data_1_name))
      
      ggplotly(proj1_plot)
      
    }else{
      proj_data_2 <- data.frame(log(proj_data_2[,which(grepl(var1, colnames(proj_data_2)))]/proj_data_2[,which(grepl(normalvar, colnames(proj_data_2)))]),
                                log(proj_data_2[,which(grepl(var2, colnames(proj_data_2)))]/proj_data_2[,which(grepl(normalvar, colnames(proj_data_2)))]),
                                proj_data_2_id,
                                proj_data_2_name
      )
      
      colnames(proj_data_2) <- c("V1",
                              "V2",
                              "Info",
                              "Group"
      )
      
      
      proj2_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
        xlab(paste0('log(',var1,'/',normalvar,')')) +
        ylab(paste0('log(',var2,'/',normalvar,')')) + 
        stat_ellipse() +  geom_point(data = proj_data_1, aes(V1, V2), shape = 18, size = 1) +  geom_point(data = proj_data_2, aes(V1, V2), shape = 15, size = 1) + ggtitle(paste("ALR projection of ", ref_data_name, proj_data_1_name ,"and", proj_data_2_name))
      
      ggplotly(proj2_plot)
    }
    
  }
  
}


### Centered log-ratio projection

## Inputs :

#Mandatory arguments for the reference dataset
# ref_data : dataframe containing compositional data only 
# ref_data_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# ref_data_group : vector defining the group each sample belongs to, corresponds to geological sources for example
# ref_data_name : character string defining the dataset, for the plot's title

#Mandatory arguments for the projection
# var1 : character string containing the name of the x variable
# var2 : character string containing the name of the y variable


#Optionnal arguments for the first projected dataset (all must be supplied for the dataset to be projected)
# proj_data_1 : dataframe containing compositional data only 
# ref_data_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# ref_data_group : vector defining the group each sample belongs to, corresponds to geological sources for example
# ref_data_name : character string defining the dataset, for the plot's title

CLR_biplot <- function(ref_data, ref_data_id, ref_data_group, ref_data_name, var1, var2, proj_data_1, proj_data_1_id, proj_data_1_name, proj_data_2, proj_data_2_id, proj_data_2_name) {
  
  geommean <- function(x) {exp(mean(log(x), na.rm = TRUE))}
  
  if(missing(proj_data_1) | missing(proj_data_1_name) | missing(proj_data_1_id)){
    
    geommean_ref_data <- apply(ref_data, 1, geommean)
    
    ref_data <- data.frame(log(ref_data[,which(grepl(var1, colnames(ref_data)))]/geommean_ref_data),
                           log(ref_data[,which(grepl(var2, colnames(ref_data)))]/geommean_ref_data),
                           ref_data_id,
                           ref_data_group
    )
    
    colnames(ref_data) <- c("V1",
                            "V2",
                            "Info",
                            "Group"
    )
    
    ref_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
      xlab(paste0('log(',var1,'/ g(x) )')) +
      ylab(paste0('log(',var2,'/ g(x) )')) + 
      stat_ellipse() + ggtitle(paste("CLR of ", ref_data_name))
    
    ggplotly(ref_plot)  
  } else {
    
    if(missing(proj_data_2) | missing(proj_data_2_name) | missing(proj_data_2_id)){
      
      common_elements <- intersect(colnames(ref_data), colnames(proj_data_1))
      
      geommean_ref_data <- apply(ref_data[,which(colnames(ref_data) %in% common_elements)],1,geommean)
      geommean_proj_1_data <- apply(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1,geommean)
      
      ref_data <- data.frame(log(ref_data[,which(grepl(var1, colnames(ref_data)))]/geommean_ref_data),
                             log(ref_data[,which(grepl(var2, colnames(ref_data)))]/geommean_ref_data),
                             ref_data_id,
                             ref_data_group
      )
      
      colnames(ref_data) <- c("V1",
                              "V2",
                              "Info",
                              "Group"
      )
      
      proj_data_1 <- data.frame(log(proj_data_1[,which(grepl(var1, colnames(proj_data_1)))]/geommean_proj_1_data),
                                log(proj_data_1[,which(grepl(var2, colnames(proj_data_1)))]/geommean_proj_1_data),
                                proj_data_1_id,
                                proj_data_1_name
      )
      
      colnames(proj_data_1) <- c("V1",
                                 "V2",
                                 "Info",
                                 "Group"
      )
      
      proj1_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
        xlab(paste0('log(',var1,'/ g(x) )')) +
        ylab(paste0('log(',var2,'/ g(x) )')) + 
        stat_ellipse() +  geom_point(data = proj_data_1, aes(V1, V2), shape = 18, size = 1) + ggtitle(paste("CLR projection of ", ref_data_name, "and", proj_data_1_name))
      
      ggplotly(proj1_plot)
      
    }else{
      
      common_elements <- Reduce(intersect, list(colnames(ref_data),colnames(proj_data_1), colnames(proj_data_2)))
      
      geommean_ref_data <- apply(ref_data[,which(colnames(ref_data) %in% common_elements)],1,geommean)
      geommean_proj_1_data <- apply(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1,geommean)
      geommean_proj_2_data <- apply(proj_data_2[,which(colnames(proj_data_2) %in% common_elements)],1,geommean)
      
      ref_data <- data.frame(log(ref_data[,which(grepl(var1, colnames(ref_data)))]/geommean_ref_data),
                             log(ref_data[,which(grepl(var2, colnames(ref_data)))]/geommean_ref_data),
                             ref_data_id,
                             ref_data_group
      )
      
      colnames(ref_data) <- c("V1",
                              "V2",
                              "Info",
                              "Group"
      )
      
      proj_data_1 <- data.frame(log(proj_data_1[,which(grepl(var1, colnames(proj_data_1)))]/geommean_proj_1_data),
                                log(proj_data_1[,which(grepl(var2, colnames(proj_data_1)))]/geommean_proj_1_data),
                                proj_data_1_id,
                                proj_data_1_name
      )
      
      colnames(proj_data_1) <- c("V1",
                                 "V2",
                                 "Info",
                                 "Group"
      )
      
      proj_data_2 <- data.frame(log(proj_data_2[,which(grepl(var1, colnames(proj_data_2)))]/geommean_proj_2_data),
                                log(proj_data_2[,which(grepl(var2, colnames(proj_data_2)))]/geommean_proj_2_data),
                                proj_data_2_id,
                                proj_data_2_name
      )
      
      colnames(proj_data_2) <- c("V1", 
                                 "V2",
                                 "Info",
                                 "Group"
      )
      
      
      
      proj2_plot <- ggplot(ref_data, aes(x=V1, V2, color=Group, label=Info)) + geom_point() + ggplot2::theme_bw() + 
        xlab(paste0('log(',var1,'/ g(x) )')) +
        ylab(paste0('log(',var2,'/ g(x) )')) +
        stat_ellipse() +  geom_point(data = proj_data_1, aes(V1, V2), shape = 18, size = 1) +  geom_point(data = proj_data_2, aes(V1, V2), shape = 15, size = 1) + ggtitle(paste("Log-ratio projection of ", ref_data_name, proj_data_1_name ,"and", proj_data_2_name))
      
      ggplotly(proj2_plot)
    }
    
  }
  
}

### PCA_plot

## Inputs :

#Mandatory arguments for the reference dataset
# ref_data : dataframe containing compositional data only 
# ref_data_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# ref_data_group : vector defining the group each sample belongs to, corresponds to geological sources for example
# ref_data_name : character string defining the dataset, for the plot's title


#Optionnal arguments for the first projected dataset (all must be supplied for the dataset to be projected)
# proj_data_1 : dataframe containing compositional data only 
# proj_data_1_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# proj_data_1_name : character string defining the dataset, for the plot's title and legend

#Optionnal arguments for the second projected dataset (all must be supplied for the dataset to be projected)
# proj_data_2 : dataframe containing compositional data only 
# proj_data_2_id : vector containing information about each compositions in the data frame, such as the sample reference for example
# proj_data_2_name : character string defining the dataset, for the plot's title and legend

# Parametric arguments
# xpca : number of the component to be plotted on the x axis (default = 1)
# ypca : number of the component to be plotted on the y axis (default = 2)
# plot_mode : boolean, TRUE by default. If true, the function will return a ggplot object, which can then be plotted with ggplotly(object), if false the function will return the the list of class 'prcomp' returned by the pca function, as well as the results for each dataset, useful for debugging and to analyse the eigenvectors


PCA_plot <- function(ref_data, ref_data_id, ref_data_group, ref_data_name, proj_data_1, proj_data_1_id, proj_data_1_name, proj_data_2, proj_data_2_id, proj_data_2_name, plot_mode = TRUE, xpca = 1, ypca = 2) {
  
  xpca_lab <- paste("PC",xpca, sep = "")
  ypca_lab <- paste("PC",ypca, sep = "")
  
  calib_pca <- function(df){
    pca_results <- prcomp(df, center = TRUE, scale. = TRUE )
    
    
    var_pca = pca_results$sdev^2 / sum(pca_results$sdev^2)
    
    return(list(pca_results, var_pca))
  }
  
  if(missing(proj_data_1) | missing(proj_data_1_name) | missing(proj_data_1_id)){
    
    pca_results <- calib_pca(normalise(ref_data,1))
    pca <- pca_results[[1]]
    var_pca <- pca_results[[2]]
    
    pca_x <- pca$x
    pca_ref_results <- cbind(pca_x[,c(xpca,ypca)], ref_data_id, ref_data_group)
    
    pca_ref_results <- as.data.frame(pca_ref_results)
    colnames(pca_ref_results) <- c("PC1", "PC2", "Info", "Group")
    
    pca_plot <- ggplot(pca_ref_results, aes(x = as.numeric(PC1), y = as.numeric(PC2), color=Group, label=Info)) + 
      xlab(paste(xpca_lab ,"(", round(var_pca[xpca], 2)*100,") %", sep="")) + 
      ylab(paste(ypca_lab ,"(", round(var_pca[ypca], 2)*100,") %", sep="")) +
      geom_point() + stat_ellipse() + theme_bw() + ggtitle(paste("PCA of", ref_data_name))
    
    if (plot_mode == TRUE){
      return(pca_plot)
    }
    if (plot_mode == FALSE){
      return(list(pca, pca_ref_results))
    }
  } else {
    
    if(missing(proj_data_2) | missing(proj_data_2_name) | missing(proj_data_2_id)){
      
      common_elements <- intersect(colnames(ref_data), colnames(proj_data_1))
      
      pca_results <- calib_pca(normalise(ref_data[,which(colnames(ref_data) %in% common_elements)],1))
      
      pca <- pca_results[[1]]
      var_pca <- pca_results[[2]]
      
      pca_x <- pca$x
      pca_ref_results <- cbind(pca_x[,c(xpca,ypca)], ref_data_id, ref_data_group)
      
      pca_proj_1_results <- predict(pca, normalise(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1))
      
      pca_proj_1_results <- cbind(pca_proj_1_results[,c(xpca,ypca)], proj_data_1_id, proj_data_1_name)
      
      #print(pca_proj_1_results)
      
      pca_ref_results <- as.data.frame(pca_ref_results)
      pca_proj_1_results <- as.data.frame(pca_proj_1_results)
      colnames(pca_ref_results) <- c("PC1", "PC2", "Info", "Group")
      colnames(pca_proj_1_results) <- c("PC1", "PC2", "Info", "Group")
      
      pca_combined_results <- rbind(pca_ref_results, pca_proj_1_results)
      
      pca_plot <- ggplot(pca_ref_results, aes(x = as.numeric(PC1), y = as.numeric(PC2), color=Group, label=Info)) + 
        xlab(paste(xpca_lab ,"(", round(var_pca[xpca], 2)*100,") %", sep="")) + 
        ylab(paste(ypca_lab ,"(", round(var_pca[ypca], 2)*100,") %", sep="")) +
        geom_point() + stat_ellipse() + theme_bw() + geom_point(data = pca_proj_1_results, aes(x = as.numeric(PC1), y = as.numeric(PC2)), shape = 18, size = 1) + ggtitle(paste("PCA of", ref_data_name, "and projection of", proj_data_1_name))
      
      
      if (plot_mode == TRUE){
        return(pca_plot)
      }
      if (plot_mode == FALSE){
        return(list(pca, pca_ref_results, pca_proj_1_results))
      }  
      
      
      
    }else{
      
      common_elements <- Reduce(intersect, list(colnames(ref_data),colnames(proj_data_1), colnames(proj_data_2)))
      
      pca_results <- calib_pca(normalise(ref_data[,which(colnames(ref_data) %in% common_elements)],1))
      
      pca <- pca_results[[1]]
      var_pca <- pca_results[[2]]
      
      pca_x <- pca$x
      pca_ref_results <- cbind(pca_x[,c(xpca,ypca)], ref_data_id, ref_data_group)
      
      pca_proj_1_results <- predict(pca, normalise(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1))
      pca_proj_1_results <- cbind(pca_proj_1_results[,c(xpca,ypca)], proj_data_1_id, proj_data_1_name)
      
      pca_proj_2_results <- predict(pca, normalise(proj_data_2[,which(colnames(proj_data_1) %in% common_elements)],1))
      pca_proj_2_results <- cbind(pca_proj_2_results[,c(xpca,ypca)], proj_data_2_id, proj_data_2_name)
      
      pca_ref_results <- as.data.frame(pca_ref_results)
      pca_proj_1_results <- as.data.frame(pca_proj_1_results)
      pca_proj_2_results <- as.data.frame(pca_proj_2_results)
      colnames(pca_ref_results) <- c("PC1", "PC2", "Info", "Group")
      colnames(pca_proj_1_results) <- c("PC1", "PC2", "Info", "Group")
      colnames(pca_proj_2_results) <- c("PC1", "PC2", "Info", "Group")
      
      
      pca_plot <- ggplot(pca_ref_results, aes(x = as.numeric(PC1), y = as.numeric(PC2), color=Group, label=Info)) + 
        xlab(paste(xpca_lab ,"(", round(var_pca[xpca], 2)*100,") %", sep="")) + 
        ylab(paste(ypca_lab ,"(", round(var_pca[ypca], 2)*100,") %", sep="")) +
        geom_point() + stat_ellipse() + theme_bw() + geom_point(data = pca_proj_1_results, aes(x = as.numeric(PC1), y = as.numeric(PC2)), shape = 18, size = 1) + geom_point(data = pca_proj_2_results, aes(x = as.numeric(PC1), y = as.numeric(PC2)), shape = 15, size = 1) + ggtitle(paste("PCA of", ref_data_name, "and projection of", proj_data_1_name, "and", proj_data_2_name))
      
      if (plot_mode == TRUE){
        return(pca_plot)
      }
      if (plot_mode == FALSE){
        return(list(pca, pca_ref_results, pca_proj_1_results, pca_proj_2_results))
      }
      
    }
    
  }
  
}


### PCA_plot_3D (experimental)

# Makes a 3D plot of first 3 components of the pca, same arguments as PCA plot with the exception of xpca and ypca not available for this function.

PCA_plot_3D <- function(ref_data, ref_data_id, ref_data_group, ref_data_name, proj_data_1, proj_data_1_id, proj_data_1_name, proj_data_2, proj_data_2_id, proj_data_2_name, plot_mode = TRUE) {
  
  calib_pca <- function(df){
    pca_results <- prcomp(df, center = TRUE, scale. = TRUE )
    
    
    var_pca = pca_results$sdev^2 / sum(pca_results$sdev^2)
    
    return(list(pca_results, var_pca))
  }
  
  if(missing(proj_data_1) | missing(proj_data_1_name) | missing(proj_data_1_id)){
    
    pca_results <- calib_pca(normalise(ref_data,1))
    pca <- pca_results[[1]]
    var_pca <- pca_results[[2]]
    
    pca_x <- pca$x
    
    pca_ref_results <- cbind(pca_x[,c(1,2,3)], ref_data_id, ref_data_group)
    
    pca_ref_results <- as.data.frame(pca_ref_results)
    colnames(pca_ref_results) <- c("PC1", "PC2","PC3", "Info", "Group")
    
    pca_plot <- plot_ly(pca_ref_results, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group, marker = list(size=2), text = ~paste("Info: ", Info))
    
    pca_plot <- pca_plot %>% add_markers()
    
    pca_plot <- pca_plot %>% layout(scene = list(xaxis = list(title = paste("PC1 (", round(var_pca[1], 2)*100,") %", sep="")),
                                                 
                                                 yaxis = list(title = paste("PC2 (", round(var_pca[2], 2)*100,") %", sep="")),
                                                 
                                                 zaxis = list(title = paste("PC3 (", round(var_pca[3], 2)*100,") %", sep=""))))
    
    if (plot_mode == TRUE){
      return(pca_plot)
    }
    if (plot_mode == FALSE){
      return(list(pca, pca_ref_results))
    }
  } else {
    
    if(missing(proj_data_2) | missing(proj_data_2_name) | missing(proj_data_2_id)){
      
      common_elements <- intersect(colnames(ref_data), colnames(proj_data_1))
      
      pca_results <- calib_pca(normalise(ref_data[,which(colnames(ref_data) %in% common_elements)],1))
      
      pca <- pca_results[[1]]
      var_pca <- pca_results[[2]]
      
      pca_x <- pca$x
      
      pca_ref_results <- cbind(pca_x[,c(1,2,3)], ref_data_id, ref_data_group)
      
      pca_ref_results <- as.data.frame(pca_ref_results)
      
      
      pca_proj_1_results <- predict(pca, normalise(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1))
      
      pca_proj_1_results <- cbind(pca_proj_1_results[,c(1,2,3)], proj_data_1_id, proj_data_1_name)
      
      #print(pca_proj_1_results)
      
      pca_ref_results <- as.data.frame(pca_ref_results)
      pca_proj_1_results <- as.data.frame(pca_proj_1_results)
      colnames(pca_ref_results) <- c("PC1", "PC2","PC3", "Info", "Group")
      colnames(pca_proj_1_results) <- c("PC1", "PC2","PC3", "Info", "Group")
      
      pca_plot <- plot_ly(pca_ref_results, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group, marker = list(size=2), text = ~paste("Info: ", Info))
      
      pca_plot <- pca_plot %>% add_markers()
      
      pca_plot <- pca_plot %>% add_trace(data = pca_proj_1_results, type='scatter3d', mode = 'markers')
     
      pca_plot <- pca_plot %>% layout(scene = list(xaxis = list(title = paste("PC1 (", round(var_pca[1], 2)*100,") %", sep="")),
                                                   
                                                   yaxis = list(title = paste("PC2 (", round(var_pca[2], 2)*100,") %", sep="")),
                                                   
                                                   zaxis = list(title = paste("PC3 (", round(var_pca[3], 2)*100,") %", sep=""))))
      
      if (plot_mode == TRUE){
        return(pca_plot)
      }
      if (plot_mode == FALSE){
        return(list(pca, pca_ref_results, pca_proj_1_results))
      }  
      
      
      
    }else{
      
      common_elements <- Reduce(intersect, list(colnames(ref_data),colnames(proj_data_1), colnames(proj_data_2)))
      
      pca_results <- calib_pca(normalise(ref_data[,which(colnames(ref_data) %in% common_elements)],1))
      
      pca <- pca_results[[1]]
      var_pca <- pca_results[[2]]
      
      pca_x <- pca$x
      pca_ref_results <- cbind(pca_x[,c(1,2,3)], ref_data_id, ref_data_group)
      
      pca_proj_1_results <- predict(pca, normalise(proj_data_1[,which(colnames(proj_data_1) %in% common_elements)],1))
      pca_proj_1_results <- cbind(pca_proj_1_results[,c(1,2,3)], proj_data_1_id, proj_data_1_name)
      
      pca_proj_2_results <- predict(pca, normalise(proj_data_2[,which(colnames(proj_data_1) %in% common_elements)],1))
      pca_proj_2_results <- cbind(pca_proj_2_results[,c(1,2,3)], proj_data_2_id, proj_data_2_name)
      
      pca_ref_results <- as.data.frame(pca_ref_results)
      pca_proj_1_results <- as.data.frame(pca_proj_1_results)
      pca_proj_2_results <- as.data.frame(pca_proj_2_results)
      colnames(pca_ref_results) <- c("PC1", "PC2","PC3", "Info", "Group")
      colnames(pca_proj_1_results) <- c("PC1", "PC2","PC3", "Info", "Group")
      colnames(pca_proj_2_results) <- c("PC1", "PC2","PC3", "Info", "Group")
      
      
      pca_plot <- plot_ly(pca_ref_results, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group, marker = list(size=2), text = ~paste("Info: ", Info))
      
      pca_plot <- pca_plot %>% add_markers()
      
      pca_plot <- pca_plot %>% add_trace(data = pca_proj_1_results, type='scatter3d', mode = 'markers')
      pca_plot <- pca_plot %>% add_trace(data = pca_proj_2_results, type='scatter3d', mode = 'markers')
      
      pca_plot <- pca_plot %>% layout(scene = list(xaxis = list(title = paste("PC1 (", round(var_pca[1], 2)*100,") %", sep="")),
                                         
                                         yaxis = list(title = paste("PC2 (", round(var_pca[2], 2)*100,") %", sep="")),
                                         
                                         zaxis = list(title = paste("PC3 (", round(var_pca[3], 2)*100,") %", sep=""))))
      
      
      if (plot_mode == TRUE){
        return(pca_plot)
      }
      if (plot_mode == FALSE){
        return(list(pca, pca_ref_results, pca_proj_1_results, pca_proj_2_results))
      }
      
    }
    
  }
  
}





