import os
import cv2
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
import numpy as np
import SOAPy_st as sp
import matplotlib.pyplot as plt
from anndata import AnnData
from tqdm import tqdm

sample = 'horse_R1_zw'
outdir = '/public/home/s20223040710/rumen/11_stereo_seq/06_velocity/'  # 根据需要设置路径
adata = sc.read_h5ad("/public/home/s20223040710/rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
#outdir = 'D:/Research/Rumen/11_stereo_seq/06_velocity/'  # 根据需要设置路径
#adata = sc.read_h5ad("D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
os.chdir(outdir)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace= True, subset = False)#  subset = True

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver="arpack")
# sce.pp.harmony_integrate(adata_concat, 'sample')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30) #use_rep="X_pca_harmony",
sc.tl.leiden(adata, resolution=0.6, key_added="leiden")#horse_R2_zw horse_R1_zw 0.6
# sc.tl.umap(adata)
# sc.pl.umap(adata)


sc.pl.spatial(adata, img_key="hires", color='leiden', spot_size=50)
plt.savefig(f'{sample}.spatial.niche.png', bbox_inches='tight', dpi=500)

# mask = sp.tl.get_mask_from_domain(adata, clusters='ctniche_5', KSize=35, cluster_key='celltype_niche')
mask = sp.tl.get_mask_from_domain(adata, clusters='4', KSize=51, cluster_key='leiden')
# plt.imshow(mask, cmap='gray')
plt.imshow(mask, cmap='gray', vmin=0, vmax=1)  # plt只能显示0到1之间的值，将0-255进行压缩
plt.savefig(f'{sample}.spatial.mask.png', bbox_inches='tight', dpi=500)

# ###Statstical testing
# # wilcoxon_res = sp.tl.wilcoxon_test(adata,mask,radius=1000,location='out',cut=500)
# # 使用 tqdm 包装计算过程
# with tqdm(total=adata.n_vars) as pbar:
#     wilcoxon_res = sp.tl.wilcoxon_test(adata,mask,radius=1000,location='out',cut=500)
#     pbar.update(adata.n_vars)  # 更新进度条
# wilcoxon_res.to_csv('wilcoxon.gene.csv')
# # spearman_res = sp.tl.spearman_correlation(adata,mask,radius=1000,num=5)
# with tqdm(total=adata.n_vars) as pbar:
#     spearman_res = sp.tl.spearman_correlation(adata,mask,radius=1000,num=5)
#     pbar.update(adata.n_vars)  # 更新进度条
# spearman_res.to_csv('spearman.gene.csv')

###Regression
sp.tl.spatial_tendency(adata,mask,method= 'poly', radius=1000,location='out',frac=5)
# sp.tl.spatial_tendency(adata,mask,method= 'poly',gene_name = 'LAP3',radius=10,location='out')
sp.pl.show_tendency(adata, gene_name = ['KRT6A','TSPYL4','LUC7L'], show=True)
plt.savefig(f'{sample}.spatial.gene.trajectory.pdf', bbox_inches='tight', dpi=500)

# ###Clustering genes based on regression curves
# sp.tl.gene_cluster(adata=adata, k=10, range_min=0.03, fdr=True, pvalue=0.05)
# sp.pl.show_curves_cluster(adata)
# plt.savefig('spatial.gene.cluster.trajectory.png',bbox_inches = 'tight',dpi = 500)

####细胞轨迹
adata.obsm['q05_cell_abundance_w_sf'] = adata.obsm['q05_cell_abundance_w_sf'].rename(
    columns=lambda col: col.replace('q05cell_abundance_w_sf_', ''))
Spatial = AnnData(X = adata.obsm['q05_cell_abundance_w_sf'],obs = adata.obs,obsm = adata.obsm ,uns = adata.uns)
# mask1 = sp.tl.get_mask_from_domain(Spatial, clusters='ctniche_5', KSize=35, cluster_key='celltype_niche')
# mask1 = sp.tl.get_mask_from_domain(adata, clusters='0', KSize=35, cluster_key='leiden')
# wilcoxon_cell = sp.tl.wilcoxon_test(Spatial,mask1,radius=1000,location='out',cut=500)
# wilcoxon_cell.to_csv('wilcoxon.cell.csv')
# spearman_cell = sp.tl.spearman_correlation(Spatial,mask1,radius=1000,num=5)
# spearman_cell.to_csv('spearman.cell.csv')
sp.tl.spatial_tendency(Spatial,mask,method='poly', radius=1000,location='out',frac=5)  #mask1

for cell in Spatial.var.index :
    sp.pl.show_tendency(adata, gene_name = cell, show=True)
    plt.savefig(f'{sample}.spatial.{cell}.trajectory.png', bbox_inches='tight', dpi=500)
    
    
sp.pl.show_tendency(adata, gene_name = Spatial.var.index, show=True)
plt.legend(title='celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('cell trajectory')
plt.ylabel('cell compositions')
plt.tight_layout()
plt.savefig(f'{sample}.spatial.cell.trajectory.pdf', bbox_inches='tight', dpi=500)

