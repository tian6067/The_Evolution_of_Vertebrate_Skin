# -*- coding: utf-8 -*-
"""
Created on Sat Nov  1 15:56:54 2025
@author: Yahui Zhang
"""
import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors
import scniche as sn
from scniche.pl import *
# set seed
sn.pp.set_seed()
warnings.filterwarnings('ignore')
print("Last run with scNiche version:", sn.__version__)

outdir = '/public/home/s20223040710/skin/03_niche_NMF/'  # 根据需要设置路径
# outdir = 'D:/Research/skin/03_niche_NMF/'  # 根据需要设置路径
os.chdir(outdir)


sample = 'adata_8'
adata = sc.read_h5ad(f'{outdir}{sample}_res0.8_k12_Niche.h5ad')


###############03 计算每个niche不同分组的差异表达基因###############################
sample_to_group = {
    'h_sheep9_cellbin_subtype': 'h_sheep_goat',
    'h_sheep12_cellbin_subtype': 'h_sheep_goat',
    'h_goat10_cellbin_subtype': 'h_sheep_goat',
    'h_goat29_cellbin_subtype': 'h_sheep_goat',
    'rg_12_cellbin_subtype': 'rg',
    'rg_10_cellbin_subtype': 'rg',
    's_x_33_cellbin_subtype': 's_x',
    's_x2_cellbin_subtype': 's_x'
}
adata.obs['group'] = adata.obs['sample'].map(sample_to_group)
niche_type = 'scNiche' #molNiche
niche_names = adata.obs[niche_type].cat.categories.tolist()

all_deg = []   # 用来收集所有 DEG 结果

for niche in niche_names:
    sub_adata = adata[adata.obs[niche_type] == niche].copy()
    sc.tl.rank_genes_groups(
        sub_adata,
        groupby="group",
        groups=["s_x"],      # 你关心的组
        reference="h_sheep_goat",     # 对照组
        method="wilcoxon"
    )
    res = sub_adata.uns["rank_genes_groups"]
    df1 = pd.DataFrame({
        "gene": res["names"]["s_x"],
        "logFC": res["logfoldchanges"]["s_x"],
        "pvals": res["pvals"]["s_x"],
        "pvals_adj": res["pvals_adj"]["s_x"],
        "niche": niche,
        "group": "s_x",
        "reference": "h_sheep_goat",
        "comparison": "s_x_vs_h_sheep_goat"
    })
    all_deg.append(df1)

    sc.tl.rank_genes_groups(
        sub_adata,
        groupby="group",
        groups=["rg"],      # 你关心的组
        reference="s_x",     # 对照组
        method="wilcoxon"
    )
    res = sub_adata.uns["rank_genes_groups"]
    df2 = pd.DataFrame({
        "gene": res["names"]["rg"],
        "logFC": res["logfoldchanges"]["rg"],
        "pvals": res["pvals"]["rg"],
        "pvals_adj": res["pvals_adj"]["rg"],
        "niche": niche,
        "group": "rg",
        "reference": "s_x",
        "comparison": "rg_vs_s_x"
    })
    all_deg.append(df2)

df_all_deg = pd.concat(all_deg, ignore_index=True)
df_all_deg_sig = df_all_deg[ (df_all_deg["pvals_adj"] < 0.05) &  (df_all_deg["logFC"] > 0)]
df_all_deg_sig.to_csv(f"DEG_by_{niche_type}.csv", index=False)




###############03 计算每个niche不同分组的差异表达基因###############################
sample_to_group = {
    'h_sheep9_cellbin_subtype': 'h_sheep',
    'h_sheep12_cellbin_subtype': 'h_sheep',
    'h_goat10_cellbin_subtype': 'h_goat',
    'h_goat29_cellbin_subtype': 'h_goat',
    'rg_12_cellbin_subtype': 'rg',
    'rg_10_cellbin_subtype': 'rg',
    's_x_33_cellbin_subtype': 's_x',
    's_x2_cellbin_subtype': 's_x'
}
adata.obs['group'] = adata.obs['sample'].map(sample_to_group)
niche_type = 'scNiche' #molNiche
niche_type = 'molNiche' #
niche_names = adata.obs[niche_type].cat.categories.tolist()

all_deg = []   # 用来收集所有 DEG 结果

for niche in niche_names:
    sub_adata = adata[adata.obs[niche_type] == niche].copy()
    sc.tl.rank_genes_groups(
        sub_adata,
        groupby="group",
        groups=["s_x"],      # 你关心的组
        reference="h_sheep",     # 对照组
        method="wilcoxon"
    )
    res = sub_adata.uns["rank_genes_groups"]
    df1 = pd.DataFrame({
        "gene": res["names"]["s_x"],
        "logFC": res["logfoldchanges"]["s_x"],
        "pvals": res["pvals"]["s_x"],
        "pvals_adj": res["pvals_adj"]["s_x"],
        "niche": niche,
        "group": "s_x",
        "reference": "h_sheep",
        "comparison": "s_x_vs_h_sheep"
    })
    all_deg.append(df1)

    sc.tl.rank_genes_groups(
        sub_adata,
        groupby="group",
        groups=["rg"],      # 你关心的组
        reference="h_goat",     # 对照组
        method="wilcoxon"
    )
    res = sub_adata.uns["rank_genes_groups"]
    df2 = pd.DataFrame({
        "gene": res["names"]["rg"],
        "logFC": res["logfoldchanges"]["rg"],
        "pvals": res["pvals"]["rg"],
        "pvals_adj": res["pvals_adj"]["rg"],
        "niche": niche,
        "group": "rg",
        "reference": "h_goat",
        "comparison": "rg_vs_h_goat"
    })
    all_deg.append(df2)

df_all_deg = pd.concat(all_deg, ignore_index=True)
df_all_deg_sig = df_all_deg[ (df_all_deg["pvals_adj"] < 0.05) & (abs(df_all_deg["logFC"]) > 0)]
df_all_deg_sig.to_csv(f"DEG_species_by_{niche_type}.csv", index=False)
