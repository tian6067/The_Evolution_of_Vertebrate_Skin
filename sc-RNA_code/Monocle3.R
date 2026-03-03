################monocle3
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)

Idents(kc)<-kc$species
cell<-list()
anno<-list()
df<-as.data.frame(table(kc$species))
df<-as.character(df$Var1)
all<-data.frame()
for (i in 1:length(df)) {
  cell[[i]]<-subset(kc,idents = df[[i]])
  Idents(cell[[i]])<-cell[[i]]$finalcelltype
  write.csv(t(as.matrix(cell[[i]]@assays$RNA@counts)),file=paste(df[[i]],".csv",sep = "_"))
}


kc$sp_finalcell<-all
table(kc$sp_finalcell)
Idents(kc)<-kc$sp_finalcell

Idents(kc)<-kc$celltype
sub<-subset(kc,idents=c('Basel_cell',  'Spinous_cell', 'Granular_cell',    'Lipid_cell'))



Idents(sub)<-sub$sp_finalcell
sce<-subset(sub,downsample=1000)

expression_matrix <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds,reduction_method='UMAP',
                        preprocess_method = 'PCA')
cds.embed <- cds@int_colData$reducedDims$UMAP#monocle3中UMAP信息
int.embed <- Embeddings(sce, reduction = "umap")#seurat中UMAP降维信息
int.embed <- int.embed[rownames(cds.embed),]#替换
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

plot_cells(cds, 
           color_cells_by = 'celltype',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,trajectory_graph_segment_size = 1)+scale_color_distiller(palette = "Spectral")
saveRDS(cds,file = 'cds.rds')
cds_fish <- cds[, colData(cds)$sp == "fish"]
cds_amp <- cds[, colData(cds)$sp == "amp"]
cds_rept <- cds[, colData(cds)$sp == "rept"]
cds_aves <- cds[, colData(cds)$sp == "aves"]
cds_mam <- cds[, colData(cds)$sp == "mam"]

cds_fish <- graph_test(cds_fish, neighbor_graph="principal_graph", cores=8)
cds_amp <- graph_test(cds_amp, neighbor_graph="principal_graph", cores=8)
cds_rept <- graph_test(cds_rept, neighbor_graph="principal_graph", cores=8)
cds_aves <- graph_test(cds_aves, neighbor_graph="principal_graph", cores=8)
cds_mam <- graph_test(cds_mam, neighbor_graph="principal_graph", cores=8)

fish_genes <- row.names(subset(cds_fish_gene, q_value< 0.05))
amp_genes <- row.names(subset(cds_amp_gene, q_value< 0.05))
rept_genes <- row.names(subset(cds_rept_gene, q_value< 0.05))
aves_genes <- row.names(subset(cds_aves_gene, q_value< 0.05))
mam_genes <- row.names(subset(cds_mam_gene, q_value< 0.05))

cutColumn_Means <- function(data_exp,#需要分割数据
                            cut#需要分割的列数
){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
  
  
}

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#####fish
plot_matrix <- exprs(cds_fish)[match(all$gene,
                                     rownames(rowData(cds_fish))),
                               order(pseudotime(cds_fish))]
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- all$gene;
dim(plot_matrix)
plot_test <- cutColumn_Means(plot_matrix,cut = 25)
p_fish<-pheatmap::pheatmap(plot_test, 
                           useRaster = T,
                           cluster_cols=FALSE, 
                           cluster_rows=T, 
                           show_rownames=F, 
                           show_colnames=F, 
                           clustering_method = "ward.D2",
                           cutree_rows=7,
                           filename=NA,
                           border_color = NA,
                           fontsize_row = 8,
                           color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                           clustering_callback = callback)
fish_module <- as.data.frame(cutree(p_fish$tree_row, k=7))
colnames(fish_module) <- "Module"
fish_module$gene <- rownames(fish_module)

#####amp
plot_matrix <- exprs(cds_amp)[match(all$gene,
                                    rownames(rowData(cds_amp))),
                              order(pseudotime(cds_amp))]
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- all$gene;
dim(plot_matrix)
plot_test <- cutColumn_Means(plot_matrix,cut = 25)
p_amp<-pheatmap::pheatmap(plot_test, 
                          useRaster = T,
                          cluster_cols=FALSE, 
                          cluster_rows=T, 
                          show_rownames=F, 
                          show_colnames=F, 
                          clustering_method = "ward.D2",
                          cutree_rows=4,
                          filename=NA,
                          border_color = NA,
                          fontsize_row = 8,
                          color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                          clustering_callback = callback)
amp_module <- as.data.frame(cutree(p_amp$tree_row, k=4))
colnames(amp_module) <- "Module"
amp_module$gene <- rownames(amp_module)


#####rept
plot_matrix <- exprs(cds_rept)[match(all$gene,
                                     rownames(rowData(cds_rept))),
                               order(pseudotime(cds_rept))]
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- all$gene;
dim(plot_matrix)
plot_test <- cutColumn_Means(plot_matrix,cut = 25)
p_rept<-pheatmap::pheatmap(plot_test, 
                           useRaster = T,
                           cluster_cols=FALSE, 
                           cluster_rows=T, 
                           show_rownames=F, 
                           show_colnames=F, 
                           clustering_method = "ward.D2",
                           cutree_rows=4,
                           filename=NA,
                           border_color = NA,
                           fontsize_row = 8,
                           color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                           clustering_callback = callback)
rept_module <- as.data.frame(cutree(p_rept$tree_row, k=4))
colnames(rept_module) <- "Module"
rept_module$gene <- rownames(rept_module)


#####aves
plot_matrix <- exprs(cds_aves)[match(sub$gene,
                                     rownames(rowData(cds_aves))),
                               order(pseudotime(cds_aves))]
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- sub$gene;
dim(plot_matrix)
plot_test <- cutColumn_Means(plot_matrix,cut = 25)

pdf('aves.pdf',height=10,width=10)
plot_test <- na.omit(plot_test)
p_aves<-pheatmap::pheatmap(plot_test, 
                           useRaster = T,
                           cluster_cols=FALSE, 
                           cluster_rows=T, 
                           show_rownames=F, 
                           show_colnames=F, 
                           clustering_method = "ward.D2",
                           cutree_rows=6,
                           filename=NA,
                           border_color = NA,
                           fontsize_row = 8,
                           color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                           clustering_callback = callback)
dev.off()

aves_module <- as.data.frame(cutree(p_aves$tree_row, k=4))
colnames(aves_module) <- "Module"
aves_module$gene <- rownames(aves_module)


#####mam
plot_matrix <- exprs(cds_mam)[match(all$gene,
                                    rownames(rowData(cds_mam))),
                              order(pseudotime(cds_mam))]
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- all$gene;
dim(plot_matrix)
plot_test <- cutColumn_Means(plot_matrix,cut = 25)

pdf('mam.pdf',height=10,width=10)
plot_test <- na.omit(plot_test)
p_mam<-pheatmap::pheatmap(plot_test, 
                          useRaster = T,
                          cluster_cols=FALSE, 
                          cluster_rows=T, 
                          show_rownames=F, 
                          show_colnames=F, 
                          clustering_method = "ward.D2",
                          cutree_rows=6,
                          filename=NA,
                          border_color = NA,
                          fontsize_row = 8,
                          color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                          clustering_callback = callback)
dev.off()

mam_module <- as.data.frame(cutree(p_mam$tree_row, k=6))
colnames(mam_module) <- "Module"
mam_module$gene <- rownames(mam_module)





plot_test <- na.omit(plot_test)
pheatmap::pheatmap(plot_test, 
                   cluster_cols=FALSE, 
                   cluster_rows=T, 
                   show_rownames=F, 
                   show_colnames=F, 
                   clustering_method = "ward.D2",
                   filename=NA,
                   border_color = NA,
                   fontsize_row = 8,
                   color=colorRampPalette(c('#5E91CB','white',"#E5839B"))(100),
                   annotation_colors=ann_colors,
                   clustering_callback = callback,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   main="Pseudotime")

gene <- c('ELOVL1','ELOVL2','ELOVL3','ELOVL4','ELOVL6','ELOVL7')
plot_genes_in_pseudotime(cds_rept[gene,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2,cell_size = 0)