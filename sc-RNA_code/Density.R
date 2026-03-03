##density 
library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
set.seed(1)
plot_umap_density <- function(data, type, base_size = 12, base_family = "") {
  ggplot(data = data, aes(x = UMAP_1, y = UMAP_2)) +
    theme_black(base_size, base_family) +
    xlim(c(min(data$UMAP_1)-2, max(data$UMAP_1)+2)) +
    ylim(c(min(data$UMAP_2)-2, max(data$UMAP_2)+2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha = 1) +
    geom_point(color = 'grey90', size = 0.05, alpha = 0.4) +
    scale_fill_viridis(option = "magma", alpha = 1) +
    labs(title = paste0("UMAP Density Plot for ", type)) +
    theme(legend.position = "none")
}
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(0, 4), "lines")
    )
}
sub$type <- sub$class
meta <- as.data.frame(sub@meta.data)
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
coord <- Embeddings(object = sub, reduction = "umap")
coord <- coord[, c(1, 2)]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord = data.frame(ID = rownames(coord), coord)
meta <- left_join(meta, coord, by = "ID")
dir.create("density/", showWarnings = FALSE)
unique_types <- unique(meta$type)
lapply(unique_types, function(tissue_type) {
  tissue_data <- meta %>% filter(type == tissue_type)
  p <- plot_umap_density(tissue_data, tissue_type)
  ggsave(filename = paste0("density/", tissue_type, "_umap_density.pdf"),
         width = 10, height = 10, plot = p)
})