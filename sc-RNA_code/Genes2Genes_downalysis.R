####Slingshot
hf <- RunSlingshot(srt = hf, group.by = "celltype", reduction = "UMAP",start = 'Hair_follicle_stem_cell_KRT15')
lineages_layers <- SCP::LineagePlot(hf, lineages = c("Lineage1","Lineage1","Lineage1","Lineage4","Lineage5"),return_layer = T)
lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer")]
SCP::CellDimPlot(hf, group.by = "celltype",pt.size = 1,
                 reduction = "UMAP",show_stat = T,
                 label = F, label_insitu = F,,palcolor = col[1:14])+
  guides(color = guide_legend(
    title = "Groups",
    override.aes = list(color= col[1:14],size=4)))+
  ggnewscale::new_scale_colour()+lineages_layers


###g2g genes
Idents(hf)<-hf$class
pr<-subset(hf,idents='no')
sec<-subset(hf,idents='yes')

Idents(pr)<-pr$Lineage1
df<-as.data.frame(pr@active.ident)
df$cell<-rownames(df)
colnames(df)<-c('cc','cell')
df <- df[!is.na(df$cc), ]
pr_sub <- subset(pr, cells = df$cell)
non_zero_genes <- rownames(pr_sub)[rowSums(GetAssayData(pr_sub, assay = "RNA", slot = "counts")) > 0]
pr_sub <- subset(pr_sub, features = non_zero_genes)

Idents(sec)<-sec$Lineage1
df<-as.data.frame(sec@active.ident)
df$cell<-rownames(df)
colnames(df)<-c('cc','cell')
df <- df[!is.na(df$cc), ]
sec_sub <- subset(sec, cells = df$cell)
non_zero_genes <- rownames(sec_sub)[rowSums(GetAssayData(sec_sub, assay = "RNA", slot = "counts")) > 0]
sec_sub <- subset(sec_sub, features = non_zero_genes)

pr_sub[["RNA"]] <- as(pr_sub[["RNA"]], "Assay")
sec_sub[["RNA"]] <- as(sec_sub[["RNA"]], "Assay")

sceasy::convertFormat('CT.h5ad', from="anndata", to="seurat", outFile='CT.rds')
sceasy::convertFormat(sec_sub, from="seurat", to="anndata", outFile='sub_pr_Lineage1.h5ad')


ncbi<-str_subset(rownames(pr_sub),pattern = '^LOC')
en<-str_subset(rownames(pr_sub),pattern = '^ENS')
all<-setdiff(c(ncbi,en),'ENSA')
filter<-setdiff(rownames(pr_sub),all)%>%as.data.frame()
colnames(filter)<-'gene'
filter1 <- filter[!grepl("-|\\.", filter$gene), ]%>%as.data.frame()
colnames(filter1)<-'gene'

ncbi<-str_subset(rownames(sec_sub),pattern = '^LOC')
en<-str_subset(rownames(sec_sub),pattern = '^ENS')
all<-setdiff(c(ncbi,en),'ENSA')
filter<-setdiff(rownames(sec_sub),all)%>%as.data.frame()
colnames(filter)<-'gene'
filter2 <- filter[!grepl("-|\\.", filter$gene), ]%>%as.data.frame()
colnames(filter2)<-'gene'


all_gene<-intersect(filter1$gene,filter2$gene)%>%as.data.frame()
colnames(all_gene)<-'gene'

#all_gene<-intersect(all_gene$gene,deg$gene)%>%as.data.frame()

colnames(all_gene)<-'gene'
write.csv(all_gene,file = 'Lineage1_gene.csv')


###g2g genes GO terms
Idents(hf)<-hf$class
hf<-subset(hf,idents=c('yes','no'))
Idents(hf)<-hf$celltype
df<-as.data.frame(table(hf$celltype))
df<-as.character(df$Var1)
markers<-list()
cell<-list()
deg<-data.frame()
for (i in 1:length(df)) {
  #cell[[i]]<-subset(hf,idents = df[[i]])
  #Idents(cell[[i]])<-cell[[i]]$class
  #markers[[i]] <- FindMarkers(cell[[i]],ident.1 = 'yes',min.pct = 0.25, logfc.threshold = 0.25)
  #markers[[i]]$cell<-df[[i]]
  #markers[[i]]$gene<-rownames(markers[[i]])
  deg<-rbind(deg,markers[[i]])
}


deg_up<-subset(deg,deg$avg_log2FC>0)
deg_down<-subset(deg,deg$avg_log2FC<0)

gene =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = deg_up$gene,
                                                 columns = 'ENTREZID',
                                                 keytype = 'SYMBOL')[,2]))
go_sec <- enrichGO(gene, OrgDb = "org.Hs.eg.db", ont="BP")
go_sec=DOSE::setReadable(go_sec, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
go_sec@result$group<-'sec'


gene =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = deg_down$gene,
                                                 columns = 'ENTREZID',
                                                 keytype = 'SYMBOL')[,2]))
go_pr <- enrichGO(gene, OrgDb = "org.Hs.eg.db", ont="BP")
go_pr=DOSE::setReadable(go_pr, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
go_pr@result$group<-'pr'


library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(aPEAR)

deg_up<-subset(deg,deg$avg_log2FC>0)
df<-as.data.frame(table(deg_up$cell))
df<-as.character(df$Var1)
t<-list()
gene<-list()
go_sec<-list()
go_sec_all<-data.frame()
kegg_sec<-list()
kegg_sec_all<-data.frame()
for (i in 1:length(df)) {
  t[[i]]<-subset(deg_up,deg_up$cell%in%df[[i]])
  gene[[i]] =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                        keys = t[[i]]$gene,
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))
  go_sec[[i]] <- enrichGO(gene[[i]], OrgDb = "org.Hs.eg.db", ont="BP")
  go_sec[[i]]=DOSE::setReadable(go_sec[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  go_sec[[i]]@result$celltype<-df[[i]]
  go_sec_all<-rbind(go_sec_all,go_sec[[i]]@result)
  
  kegg_sec[[i]] <- enrichKEGG(gene[[i]])
  kegg_sec[[i]]=DOSE::setReadable(kegg_sec[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  kegg_sec[[i]]@result$celltype<-df[[i]]
  kegg_sec_all<-rbind(kegg_sec_all,kegg_sec[[i]]@result)
}

deg_down<-subset(deg,deg$avg_log2FC<0)
df<-as.data.frame(table(deg_down$cell))
df<-as.character(df$Var1)
t<-list()
gene<-list()
go_pr<-list()
go_pr_all<-data.frame()
kegg_pr<-list()
kegg_pr_all<-data.frame()
for (i in 1:length(df)) {
  t[[i]]<-subset(deg_down,deg_down$cell%in%df[[i]])
  gene[[i]] =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                        keys = t[[i]]$gene,
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))
  go_pr[[i]] <- enrichGO(gene[[i]], OrgDb = "org.Hs.eg.db", ont="BP")
  go_pr[[i]]=DOSE::setReadable(go_pr[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  go_pr[[i]]@result$celltype<-df[[i]]
  go_pr_all<-rbind(go_pr_all,go_pr[[i]]@result)
  
  kegg_pr[[i]] <- enrichKEGG(gene[[i]])
  kegg_pr[[i]]=DOSE::setReadable(kegg_pr[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  kegg_pr[[i]]@result$celltype<-df[[i]]
  kegg_pr_all<-rbind(kegg_pr_all,kegg_pr[[i]]@result)
}


go_pr_all<-subset(go_pr_all,go_pr_all$p.adjust<0.05)
go_sec_all<-subset(go_sec_all,go_sec_all$p.adjust<0.05)

go_pr_all_top5 <- go_pr_all %>% group_by(celltype) %>% top_n(n = 5, wt = Count)
go_sec_all_top5 <- go_sec_all %>% group_by(celltype) %>% top_n(n = 5, wt = Count)


go_sec_all_sub <- subset(go_sec_all, go_sec_all$celltype %in%  c('Dermal_papilla','Duct_cells','Hair_follicle_stem_cell_KRT15',
                                                                 'Hair_follicle_stem_cell_KRT15_LGR5','Hair_follicle_stem_cell_LGR5','Hair_Matrix',
                                                                 'Inner_root_sheath', 'Outer_root_sheath', 'Proliferating_cells_matrix','Sebaceous_gland',
                                                                 'Sweat_cells', 'Transit_amplifying_cells'))
go_pr_all_sub <- subset(go_pr_all, go_pr_all$celltype %in% c('Interfollicle_epidermis','Proliferating_cells_basal'))
go_pr_all_top5 <- go_pr_all_sub %>% group_by(celltype) %>% top_n(n = 100, wt = Count)

go_sec_all_sub<-as.data.frame(table(go_sec_all_sub$Description))
go_pr_all_sub<-as.data.frame(table(go_pr_all_sub$Description))


deg_up<-subset(deg,deg$avg_log2FC>0)
df<-as.data.frame(table(deg_up$cell))
df<-as.character(df$Var1)
t<-list()
gene<-list()
kegg_sec<-list()
kegg_sec_all<-data.frame()
for (i in 1:length(df)) {
  t[[i]]<-subset(deg_up,deg_up$cell%in%df[[i]])
  gene[[i]] =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                        keys = t[[i]]$gene,
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))
  
  kegg_sec[[i]] <- enrichKEGG(gene[[i]])
  kegg_sec[[i]]=DOSE::setReadable(kegg_sec[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  kegg_sec[[i]]@result$celltype<-df[[i]]
  kegg_sec_all<-rbind(kegg_sec_all,kegg_sec[[i]]@result)
}


deg_down<-subset(deg,deg$avg_log2FC<0)
df<-as.data.frame(table(deg_down$cell))
df<-as.character(df$Var1)
t<-list()
gene<-list()
kegg_pr<-list()
kegg_pr_all<-data.frame()
for (i in 1:length(df)) {
  t[[i]]<-subset(deg_down,deg_down$cell%in%df[[i]])
  gene[[i]] =as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                        keys = t[[i]]$gene,
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))
  kegg_pr[[i]] <- enrichKEGG(gene[[i]])
  kegg_pr[[i]]=DOSE::setReadable(kegg_pr[[i]], OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  kegg_pr[[i]]@result$celltype<-df[[i]]
  kegg_pr_all<-rbind(kegg_pr_all,kegg_pr[[i]]@result)
}


kegg_pr_all<-subset(kegg_pr_all,kegg_pr_all$category%in%c('Cellular Processes','Environmental Information Processing','Genetic Information Processing','Metabolism','Organismal Systems'))
kegg_sec_all<-subset(kegg_sec_all,kegg_sec_all$category%in%c('Cellular Processes','Environmental Information Processing','Genetic Information Processing','Metabolism','Organismal Systems'))

kegg_pr_all<-subset(kegg_pr_all,kegg_pr_all$pvalue<0.05)
kegg_sec_all<-subset(kegg_sec_all,kegg_sec_all$pvalue<0.05)

kegg_pr_all_top5 <- kegg_pr_all %>% group_by(celltype) %>% top_n(n = 5, wt = Count)
kegg_sec_all_top5 <- kegg_sec_all %>% group_by(celltype) %>% top_n(n = 5, wt = Count)


kegg_sec_all_sub <- subset(kegg_sec_all, kegg_sec_all$celltype %in%  c('Dermal_papilla','Duct_cells','Hair_follicle_stem_cell_KRT15',
                                                                       'Hair_follicle_stem_cell_KRT15_LGR5','Hair_follicle_stem_cell_LGR5','Hair_Matrix',
                                                                       'Inner_root_sheath', 'Outer_root_sheath', 'Proliferating_cells_matrix','Sebaceous_gland',
                                                                       'Sweat_cells', 'Transit_amplifying_cells'))
kegg_pr_all_sub <- subset(kegg_pr_all, kegg_pr_all$celltype %in% c('Interfollicle_epidermis','Proliferating_cells_basal'))
kegg_pr_all_top5 <- kegg_pr_all_sub %>% group_by(celltype) %>% top_n(n = 10, wt = Count)

kegg_sec_all_sub<-as.data.frame(table(kegg_sec_all_sub$Description))
kegg_pr_all_sub<-as.data.frame(table(kegg_pr_all_sub$Description))










