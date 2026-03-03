# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 18:22:02 2025
"""
import scanpy as sc
import os
import pandas as pd
import numpy as np
indir = 'D:/Research/skin/02_cell2location/'  # 根据需要设置路径
outdir = 'D:/Research/skin/03_niche_NMF'  # 根据需要设置路径
os.chdir(outdir)

sample_names = ["h_goat10_cellbin_subtype", "h_goat29_cellbin_subtype", "h_sheep9_cellbin_subtype", "h_sheep12_cellbin_subtype",
           "rg_12_cellbin_subtype", "rg_10_cellbin_subtype", "s_x_33_cellbin_subtype", "s_x2_cellbin_subtype", 
          ]
#abundance_list = [sc.read(indir + sample + '_cell2loc_anno.h5ad').obsm['q05_cell_abundance_w_sf'] for sample in sample_names ]
abundance_list = [pd.read_csv(indir  + sample + '_predicted_cell_types.csv', index_col=0) for sample in sample_names ]
all_columns = sorted(set().union(*[df.columns for df in abundance_list]))

# 一次性合并所有DataFrame
X_data = pd.concat(
    [df.reindex(columns=all_columns).fillna(0) for df in abundance_list], 
    axis=0)
X_data.to_csv('X_data.csv', index=True)


adata_list = []
for sample in sample_names:
    adata_sp = sc.read(indir + sample + '_cell2loc_anno.h5ad')
    adata_sp.layers['counts'] = adata_sp.X.copy()
    adata_list.append(adata_sp)
    
adata_combined = sc.concat(adata_list, label='sample', merge='first', join = 'inner', keys=sample_names)

adata_combined.obsm['q05_cell_abundance_w_sf'] = X_data

cell_type_predictions = adata_combined.obsm['q05_cell_abundance_w_sf']
# 使用 np.argmax() 获取每一行最大值的索引
predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
# 获取列名列表，即细胞类型的名称
cell_type_names = cell_type_predictions.columns
# 根据最大值的索引选择细胞类型的名称
predicted_cell_types = cell_type_names[predicted_cell_types_idx]
# 将预测的细胞类型添加到 adata_sp.obs 中
adata_combined.obs['predicted_cell_type'] = predicted_cell_types
adata_combined.obs['predicted_cell_type'] = adata_combined.obs['predicted_cell_type'].str.replace('^q05cell_abundance_w_sf_', '', regex=True)


adata_combined.write('adata_8.h5ad')

