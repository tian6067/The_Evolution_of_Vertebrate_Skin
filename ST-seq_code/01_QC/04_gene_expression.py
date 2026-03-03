import os
import scanpy as sc
import scanpy.external as sce

dir = 'D:/Research/skin/01_read_ref/'
os.chdir(dir)


########## 02 逐个切片展示基因表达量#############
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

colors = ["#E6E6E6", "#FF0000"]  # gray → red
custom_cmap = LinearSegmentedColormap.from_list("gray_to_red", colors)

#CLDN1
file_names=[
 'h_goat10_cellbin_QC.h5ad',
 'h_goat29_cellbin_QC.h5ad',
 'h_sheep12_cellbin_QC.h5ad',
 'h_sheep9_cellbin_QC.h5ad',
 'rg_10_cellbin_QC.h5ad',
 'rg_12_cellbin_QC.h5ad',
 's_x2_cellbin_QC.h5ad',
 's_x_33_cellbin_QC.h5ad']


 
for i in file_names:
    adata=sc.read('D:/BT-KYFW-2023-8028-分析结果\\ST\\cellbin\\02_QC\\'+i)
    adata.layers["counts"] = adata.X.copy()

    #process
    sc.pp.normalize_total(adata, inplace=True) # 所有细胞counts的中位数作为target_sum
    sc.pp.log1p(adata)

    adata.layers["log"] = adata.X.copy()
    
    sc.pl.spatial(
                    adata,
                    color='SPINK5',  # 需要展示的基因
                    #groups=['ELOVL3'],
                    spot_size=50,
                    size=0.8,
                    ncols=3,
                    vmin=0,      # 最小值（colorbar 下限）
                    #vmax=10,     # 最大值（colorbar 上限）
                    #save='snake10.pdf'
                    #layer='counts',
                    #title='ELOVL3',
                    #palette= celltype_colors,
                    save=f'_{i}.pdf' ,
                    cmap=custom_cmap
                    )
        