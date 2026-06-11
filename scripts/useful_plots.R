library(tidyverse)

pca_plot <- function(pca_data, coloring = pca_data[[1]], title = 'PCA'){
  X <- scale(pca_data[, -1])
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  var_explained <- summary(pca)$importance[2,]
  
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    coloring = coloring
  )
  
  p <- ggplot(pca_df, aes(PC1, PC2, color = coloring)) +
        geom_point(alpha = 0.6) +
        theme_bw() +
        labs(x = paste0('PC1 (', round(var_explained[1]*100,1), '%)'),
             y = paste0('PC2 (', round(var_explained[2]*100,1), '%)'),
             title = title)
  
  invisible(list(plot = p, pca = pca))
}


pca_biplot <- function(pca_data, coloring = pca_data[[1]], scale_arrows = 10, compnames = NULL, title = 'PCA biplot'){
  X <- scale(pca_data[, -1])
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  var_explained <- summary(pca)$importance[2,]
  
  # sample coordinates
  scores <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    coloring = coloring
  )
  
  # motif vectors
  loadings <- data.frame(
    PC1 = pca$rotation[,1],
    PC2 = pca$rotation[,2],
    comp_names = if (is.null(compnames)) rownames(pca$rotation) else compnames)
  
  # scale the loadings (smaller than sample coordinates) to make arrows more visible
  loadings$PC1 <- loadings$PC1 * scale_arrows
  loadings$PC2 <- loadings$PC2 * scale_arrows
  
  library(ggrepel)
  
  p <- ggplot(scores, aes(PC1, PC2, color = coloring)) +
          geom_point(alpha = 0.7, size = 2) +
          geom_segment(data = loadings,
                       aes(x = 0, y = 0, xend = PC1, yend = PC2),
                       arrow = arrow(length = unit(0.25,'cm')),
                       inherit.aes = FALSE,
                       color = 'black') +
          geom_text_repel(data = loadings,
                          aes(x = PC1, y = PC2, label = comp_names),
                          inherit.aes = FALSE,
                          size = 3.5,
          ) +
          theme_bw() +
          labs(
            title = title,
            x = paste0('PC1 (', round(var_explained[1]*100,1), '%)'),
            y = paste0('PC2 (', round(var_explained[2]*100,1), '%)')
          )
  
  invisible(list(plot = p, pca = pca))
}
