library(tidyverse)
library(cowplot)

data<-AverageExpression(merged,assays = 'RNA',add.ident = 'class',features = gene$int)
data<-data[["RNA"]]
t_mat_PCA <- t(data)

# Running PCA
PCA <- prcomp(t_mat_PCA[], center = TRUE,scale. = TRUE)

# Getting PCA variance
summary(PCA)$importance %>%
  as.data.frame() -> Variance_PCA

# Getting 10 first components
Variance_PCA[2,] %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  dplyr::rename(components = rowname) %>% 
  dplyr::rename(Proportion_of_Variance = "Proportion of Variance") %>% 
  filter(components == "PC1" |
           components == "PC2" |
           components == "PC3" |
           components == "PC4" |
           components == "PC5" |
           components == "PC6" |
           components == "PC7" |
           components == "PC8" |
           components == "PC9" |
           components == "PC10" ) -> Variance_PCA

# Ordering PCs
Variance_PCA$components <- factor( Variance_PCA$components , levels=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10") )

# Plotting variance
Variance_PCA %>% 
  ggplot() +
  geom_bar( aes(components,Proportion_of_Variance*100),stat = "identity" )+
  geom_text(aes(x=components,y=Proportion_of_Variance*100,label=sprintf("%0.1f", round(Proportion_of_Variance*100, digits = 1))),vjust=-0.5) +
  theme_cowplot()



PCA$x %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  mutate(Species = case_when( grepl("Mammalia$",rowname) ~ "Mammalia",
                              grepl("Aves$",rowname) ~ "Aves",
                              grepl("Reptilia$",rowname) ~ "Reptilia",
                              grepl("Amphibia$",rowname) ~ "Amphibia",
                              grepl("Pisces$",rowname) ~ "Pisces"))%>% 
  mutate(Cell_type = case_when( grepl("Differentiated-KC",rowname) ~ "Differentiated_KC",
                                grepl("Endothelial",rowname) ~ "Endothelial",
                                grepl("Fibroblast",rowname) ~ "Fibroblast",
                                grepl("Goblet_cell",rowname) ~ "Goblet_cell",
                                grepl("Hair-follicle-matrix-cells",rowname) ~ "Hair_follicle_matrix_cells",
                                grepl("Hair-follicle-stem-cells",rowname) ~ "Hair_follicle_stem_cells",
                                grepl("Immune",rowname) ~ "Immune",
                                grepl("Ionocytes",rowname) ~ "Ionocytes",
                                grepl("Melanocytes",rowname) ~ "Melanocytes",
                                grepl("Merkel-cell",rowname) ~ "Merkel_cell",
                                grepl("Mucous-cell",rowname) ~ "Mucous_cell",
                                grepl("Proliferating",rowname) ~ "Proliferating",
                                grepl("Schwann",rowname) ~ "Schwann",
                                grepl("Sebaceous-gland-channel",rowname) ~ "Sebaceous_gland_channel",
                                grepl("Sebaceous-gland-secretion",rowname) ~ "Sebaceous_gland_secretion",
                                grepl("SMC-Pericytes",rowname) ~ "SMC_Pericytes",
                                grepl("Undifferentiated-KC",rowname) ~ "Undifferentiated_KC")) -> PCA_PC

PCA_PC %>% 
  ggplot() +
  geom_point( aes(PC1,PC2,col=Species),size=6,stroke = 3) +
  theme_cowplot()

with(PCA_PC, scatter3D(x = PCA_PC$PC2, y = PCA_PC$PC3, z = PCA_PC$PC1,
                       pch = 21, cex = 1.5,
                       col = "white",
                       xlab = "PC2",
                       ylab = "PC3",
                       zlab = "PC1", 
                       bg = PCA_PC$col,
                       ticktype = "detailed",
                       bty = "f", box = T,
                       theta = -20, phi = 0, d=3,
                       colkey = FALSE))