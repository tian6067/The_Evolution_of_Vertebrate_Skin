import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy.external as sce


dir = 'D:/Research/skin/01_read_ref/'
os.chdir(dir)

sample_names = ["h_sheep9_cellbin", "h_sheep12_cellbin", "h_goat10_cellbin", "h_goat29_cellbin", 
          "rg_10_cellbin", "rg_12_cellbin", "s_x_33_cellbin", "s_x2_cellbin"] #对cellbin进行质控

sample_names = ["chickenL_r1_cellbin", "chickenL_r2_cellbin",
                "axolotl_r1_cellbin",  "axolotl_r2_cellbin",
                "goldfish_r1_cellbin", "goldfish_r2_cellbin", 
                "snake_r1_cellbin", "snake_r2_cellbin", "snake_r3_cellbin", 
                "cy2_cellbin","cy4_cellbin"
                ] #对分片后的cellbin进行质控

MTgenes_dict = {
    "h_sheep9_cellbin":  ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "h_sheep12_cellbin": ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "s_x_33_cellbin":    ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "s_x2_cellbin":      ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "s_x_33_cellbin":    ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "s_x2_cellbin":      ["ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],

    "h_goat10_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP6"],
    "h_goat29_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP6"],
    "rg_10_cellbin":    ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP6"],
    "rg_12_cellbin":    ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP6"],
    
    "goldfish_cellbin": ["ND6","ND5","ND4L","ND4L","MT-ND4","MT-ND3","MT-ND2","ND1"],
    "chickenL_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "cy2_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "cy4_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "axolotl_cellbin": [],
    "snake_cellbin": [],
    
    "chickenL_r1_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "chickenL_r2_cellbin": ["ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6"],
    "goldfish_r1_cellbin": ["ND6","ND5","ND4L","ND4L","MT-ND4","MT-ND3","MT-ND2","ND1"],
    "goldfish_r2_cellbin": ["ND6","ND5","ND4L","ND4L","MT-ND4","MT-ND3","MT-ND2","ND1"],
    "axolotl_r1_cellbin": [],
    "axolotl_r2_cellbin": [],
    "snake_r1_cellbin": [],
    "snake_r2_cellbin": [],
    "snake_r3_cellbin": [],
}

# 循环读取每一行数据
for sample in sample_names:
    # 获取质控标准
    tresh = {
        'pct_counts_mt': 10,
        'total_counts_low': 100,
        'n_genes_by_counts_low': 50,
    } 
    adata = sc.read_h5ad(sample +'.h5ad')
    n0 = adata.shape[1]
    adata.obs['orig.ident'] = sample
    adata.obs_names = [sample + '_' + x for x in adata.obs_names]
    ####2.质控#############
    #2.2 过滤低质量读数细胞（每个条形码的计数数量（计数深度）、每个条形码的基因数量、每个条形码的线粒体基因计数比例）
    adata.var["mt"] = adata.var_names.isin(MTgenes_dict[sample])
    #2.2.1.1 计算 QC 协变量或度量
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[300], log1p=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save='_' + sample +'_plot.png')

    #2.2.3.1 手动过滤低质量读数的细胞  我们首先定义一个过滤字典 
    n1 = adata.shape[0]
    adata.obs['passing_mt'] = adata.obs['pct_counts_mt'] < tresh['pct_counts_mt']
    adata.obs['passing_nUMIs'] = (adata.obs['total_counts'] > tresh['total_counts_low']) 
    adata.obs['passing_ngenes'] = (adata.obs['n_genes_by_counts'] > tresh['n_genes_by_counts_low'])
    print(f'Lower treshold, total_counts: {tresh["total_counts_low"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_nUMIs"])}')
    print(f'Lower treshold, n genes: {tresh["n_genes_by_counts_low"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_ngenes"])}')
    print(f'Lower treshold, mito %: {tresh["pct_counts_mt"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_mt"])}')
    #我们需要对数据取保留的细胞的交集
    QC_test = (adata.obs['passing_mt']) & (adata.obs['passing_nUMIs']) & (adata.obs['passing_ngenes'])
    removed = QC_test.loc[lambda x : x == False]
    adata = adata[QC_test, :].copy()
    n2 = adata.shape[0]
    print(f'Cells retained after filtering of low quality cells: {n2}, {n1-n2} removed.')

    n3 = adata.shape[1]
    sc.pp.filter_genes(adata, min_cells=3)  # 只有在至少有3个细胞表达的基因才会被保留
    n4 = adata.shape[1]
    print(f'Genes retained after filtering min_cells: {sample}, {n4}, {n3-n4} removed,{n3} before.\n\n')
    sc.pp.filter_cells(adata, min_counts=10)  # 每个细胞至少需要表达10个基因

    # 确保这些对象是数值类型
    sc.pl.spatial(adata, color=['n_genes_by_counts','total_counts','pct_counts_mt'],spot_size=50, size=1.3, cmap="tab20")
    adata.obs['n_genes_by_counts'] = pd.to_numeric(adata.obs['n_genes_by_counts'], errors='coerce')
    adata.obs['total_counts'] = pd.to_numeric(adata.obs['total_counts'], errors='coerce')
    adata.obs['pct_counts_mt'] = pd.to_numeric(adata.obs['pct_counts_mt'], errors='coerce')
   
    adata.write(dir+ sample + '_QC.h5ad')
