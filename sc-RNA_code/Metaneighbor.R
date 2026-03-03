####MetaNeighbor
library(MetaNeighbor)
kc<-as.SingleCellExperiment(kc)
var_genes = variableGenes(dat = kc, exp_labels = kc$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = pbmc,
                             study_id = pbmc$orig.ident,
                             cell_type = pbmc$ident,
                             fast_version = TRUE)

df<-as.data.frame(table(merged$species))
get_species <- function(names) sapply(strsplit(names, "\\|"), `[`, 1)
row_species <- get_species(rownames(MetaNeighbor_species))
col_species <- get_species(colnames(MetaNeighbor_species))
row_order <- order(match(row_species, df$Var1))
col_order <- order(match(col_species, df$Var1))
ordered_matrix <- MetaNeighbor_species[row_order, col_order]

row_names <- rownames(ordered_matrix)
cell_types <- sapply(strsplit(row_names, "\\|"), `[`, 2)
cell_types<-as.data.frame(cell_types)
df<-as.data.frame(table(merged$celltype))
cell_types$cell_types<-factor(cell_types$cell_types,levels =df$Var1)

order <- order(cell_types$cell_types)
sorted_matrix <- ordered_matrix[order, order]
plotHeatmap(sorted_matrix,cex=0.5,Rowv = FALSE, Colv = FALSE, cexRow = 0.8, cexCol = 0.8)