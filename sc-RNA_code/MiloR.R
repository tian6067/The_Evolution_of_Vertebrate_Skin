library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
sce<-readRDS("/public/home/2020304010114/A_SKIN/re_KC/KC_final_12_11.rds")
sce <- as.SingleCellExperiment(sce)
sce_milo <- Milo(sce)  

sce_milo <- buildGraph(sce_milo, k = 500, d = 30,reduced.dim="INTEGRATED.RPCA")

sce_milo <- makeNhoods(sce_milo, prop = 0.05, k =500, d=30, refined = TRUE,reduced_dims="INTEGRATED.RPCA",refinement_scheme="graph") 

#plotNhoodSizeHist(scRNA_pre) 

sce_milo <- countCells(sce_milo, meta.data = as.data.frame(colData(sce_milo)), samples="orig.ident") 

thy.design <- data.frame(colData(sce_milo))[,c("orig.ident", "sp")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$orig.ident
## Reorder rownames to match columns of nhoodCounts(milo)
thy.design <- thy.design[colnames(nhoodCounts(sce_milo)), , drop=FALSE]
table(thy.design$sp)

model <- model.matrix(~ 0 + sp, data=thy.design)
ave.contrast <- c("spmam-( spamp + spfish + 3*sprept +10*spaves)/15 ")
mod.constrast <- makeContrasts(contrasts=ave.contrast, levels=model)
mod.constrast

da_results <- testNhoods(sce_milo, design = ~ 0 + sp, design.df = thy.design, reduced.dim = "INTEGRATED.RPCA",model.contrasts = ave.contrast, fdr.weighting="graph-overlap")

#table(da_results$SpatialFDR < 0.1)

sce_milo <- buildNhoodGraph(sce_milo)

pdf("/public/home/2020304010114/picture/adjust_miloR_pdf/k500_umap_mam.pdf")
plotUMAP(sce_milo, colour_by="finalcelltype") + plotNhoodGraphDA(sce_milo, da_results, alpha=0.1) +
  plot_layout(guides="auto" )
dev.off()

da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "finalcelltype")
pdf("/public/home/2020304010114/picture/adjust_miloR_pdf/k500_dot_mam.pdf")
plotDAbeeswarm(da_results, group.by = "finalcelltype")
dev.off()