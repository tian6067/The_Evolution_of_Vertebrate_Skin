# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 22:24:04 2025
@author: Yahui Zhang
"""
import os
import scniche as sn
import scanpy as sc
import pandas as pd
import numpy as np
import scanpy.external as sce
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
print("Last run with scNiche version:", sn.__version__)
# set seed
sn.pp.set_seed()

outdir = '/public/home/s20223040710/skin/04_scNiche/'  # 根据需要设置路径
# outdir = 'D:/Research/skin/04_scNiche/'  # 根据需要设置路径
os.chdir(outdir)

sample = 'adata_8'
adata = sc.read_h5ad(sample +".h5ad")#merge

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
adata.raw = adata.copy()
sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True)

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

sce.pp.harmony_integrate(adata, 'sample')
# sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=10, n_pcs=30)

#设置细胞类型
cell_type_predictions = adata.obsm['q05_cell_abundance_w_sf']
# 使用 np.argmax() 获取每一行最大值的索引
predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
# 获取列名列表，即细胞类型的名称
cell_type_names = cell_type_predictions.columns
# 根据最大值的索引选择细胞类型的名称
predicted_cell_types = cell_type_names[predicted_cell_types_idx]
# 将预测的细胞类型添加到 adata_sp.obs 中
adata.obs['predicted_cell_type'] = predicted_cell_types

columns_to_keep = ['predicted_cell_type', 'sample']
adata.obs = adata.obs[columns_to_keep]
adata.obsm = {key: adata.obsm[key] for key in ['q05_cell_abundance_w_sf','X_pca_harmony','X_pca', 'spatial']}
adata.obsm['q05_cell_abundance_w_sf'] = adata.obsm['q05_cell_abundance_w_sf'].fillna(0)

celltype_key='predicted_cell_type' # not use, just for `sn.pp.process_multi_slices`
sample_key = 'sample'
use_rep = 'X_pca_harmony'
k_cutoff = 15
lr = 0.01
epochs = 100
batch_num = 8
# prepare multi slices
adata = sn.pp.process_multi_slices(
    adata=adata,
    celltype_key=celltype_key,
    sample_key=sample_key,
    mode='KNN',
    k_cutoff=k_cutoff,
    is_pca=False,
    verbose=True,
    layer_key=use_rep
)

# choose the features of the three views `X_C2L`, `X_data`, and `X_data_nbr` to run scNiche     
choose_views = ['q05_cell_abundance_w_sf', 'X_data', 'X_data_nbr']
# 强制转换所有特征矩阵为numpy数组
for view in choose_views:
    if view in adata.obsm:
        original_type = type(adata.obsm[view])
        if hasattr(adata.obsm[view], 'values'):
            adata.obsm[view] = adata.obsm[view].values
        else:
            adata.obsm[view] = np.array(adata.obsm[view])
        
adata = sn.pp.prepare_data_batch(adata=adata, verbose=True, batch_num=batch_num, choose_views=choose_views)
# training
model = sn.tr.Runner_batch(adata=adata, device='cuda:0', verbose=False, choose_views=choose_views)
adata = model.fit(lr=0.01, epochs=epochs)
adata.write(f'{sample}_model_scNiche.h5ad')



