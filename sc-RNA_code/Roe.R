library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(Seurat)
library(tidyverse)
library(readr)
library(qs)
data<-as.data.frame(sc@meta.data)

R_oe<-calTissueDist(data,
                    byPatient=F,
                    colname.cluster='seurat_clusters',
                    colname.patient="orig.ident",
                    colname.tissue="sp",
                    method="chisq",
                    min.rowSum=0)

col_fun <- colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), 
                      c("blue", "white", "red"))

Heatmap(as.matrix(R_oe),
        show_heatmap_legend = TRUE, 
        cluster_rows = T, 
        cluster_columns = F,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",  # 自定义图注名称
          at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
          labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)
