import os
import scniche as sn
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
from scniche.pl import *
import warnings
warnings.filterwarnings('ignore')
print("Last run with scNiche version:", sn.__version__)

outdir = '/public/home/s20223040710/skin/04_scNiche/'  # 根据需要设置路径
# outdir = 'D:/Research/skin/04_scNiche/'  # 根据需要设置路径
os.chdir(outdir)
# set seed
sn.pp.set_seed()
sample = 'adata_8'


samples = ["h_sheep9_cellbin_subtype", "h_sheep12_cellbin_subtype", "h_goat10_cellbin_subtype", "h_goat29_cellbin_subtype", 
           "rg_12_cellbin_subtype", "rg_10_cellbin_subtype", "s_x_33_cellbin_subtype", "s_x2_cellbin_subtype"]
for target_num in range(12, 16):
    print(f"\n=== k={target_num} ===")
    adata= sc.read_h5ad(f'{sample}_model_scNiche.h5ad')
    adata = sn.tr.clustering(adata=adata, target_k=target_num)
    palette = {
        "Hair_follicle_stem_cell_KRT15": "#b8b8d8",
        "Hair_follicle_stem_cell_KRT15_LGR5": "#e0b6d5",
        "Hair_follicle_stem_cell_LGR5": "#b58bb9",
        "Hair_Matrix": "#dfd8a9",
        "Inner_root_sheath": "#7dd1e7",
        "Interfollicle_epidermis": "#f4a796",
        "Dermal_papilla": "#f2bc95",
        "Duct_cells": "#91b7e1",
        "Outer_root_sheath": "#dfd8a9",
        "Proliferating_cells_basal": "#f9c3c6",
        "Proliferating_cells_matrix": "#f8dac5",
        "Sebaceous_gland": "#a6dad4",
        "Sweat_cells": "#b9d1b3",
        "Proliferating": "#f8dac5",  # Added for completeness, this was likely meant to be Proliferating_cells_matrix
        "Transit_amplifying_cells": "#f49896",
        "Angiogenic_EC": "#e0cc8c",
        "Arteriole_EC": "#da9d75",
        "capillary_EC": "#ce95be",
        "LEC": "#78b2a9",
        "Post_capillary_venule": "#e7b7d3",
        "Venule_EC": "#66a6d1",
        "SMC_Pericytes": "#f2af2b",
        "Mesenchymal": "#e6898a",
        "Pro_inflammatory": "#9c77ad",
        "Secretory_papillary": "#6f9cb7",
        "Secretory_reticular": "#96c780",
        "Langerhans_cell": "#c0c7e2",
        "T_cell": "#c9e1c8",
        "Macrophages": "#f9c8b8",
        "Mast_cell": "#b1c6e7",
        "NKT": "#d3c9c1",
        "Dendritic_cells": "#f1d3e1",
        
    ########空间NNiche的颜色定义#########
        "Niche0": "#94ffb5",
        "Niche1": "#000000",
        "Niche2": "#8f7c00",
        "Niche3": "#dfd8a9",
        "Niche4": "#ffcc99",
        "Niche5": "#2bce48",
        "Niche6": "#993f00",
        "Niche7": "#005c31",
        "Niche8": "#0075dc",
        "Niche9": "#f0a0ff",
        "Niche10": "#c20088",
        "Niche11": "#4c005c",
        "Niche12": "#9dcc00",
        "Niche13": "#f49896",
        "Niche14": "#e0cc8c",
    }
    for i in samples:
        sc.pl.spatial(adata[adata.obs['sample'] ==  i], color=['scNiche', 'predicted_cell_type'],
                      spot_size=25,title=f"{i}", size=1.3, palette=palette, save=f"_k{target_num}_{i}_scniche.pdf")
    # 计算ARI
    res = adata.obs.copy()
    for i in samples:
        res_tmp = res.loc[res['sample'] == i]
        ari_tmp = adjusted_rand_score(res_tmp['predicted_cell_type'], res_tmp['scNiche'])
        print(f"{i}: {ari_tmp:.4f}")
    adata.write(f'{sample}_k{target_num}_scNiche.h5ad')


#####################读取文件################
sample = 'adata_8'
adata= sc.read_h5ad(f'{sample}_scNiche.h5ad')
adata= sc.read_h5ad(f'{sample}_model_scNiche.h5ad')




####################堆叠柱状图##################
cell_counts = pd.crosstab(adata.obs['sample'], adata.obs['scNiche'])
cell_ratios = cell_counts.div(cell_counts.sum(axis=1), axis=0)
cell_ratios.to_csv('scNiche_ratios.csv')

cell_counts = pd.crosstab(adata.obs['sample'], adata.obs['molNiche'])
cell_ratios = cell_counts.div(cell_counts.sum(axis=1), axis=0)
cell_ratios.to_csv('molNiche_ratios.csv')


####################堆叠柱状图 按照三个进行分组##################
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
cell_counts = pd.crosstab(adata.obs['group'], adata.obs['molNiche'])
cell_ratios = cell_counts.div(cell_counts.sum(axis=1), axis=0)
cell_ratios.to_csv('molNiche_group3_ratios.csv')


cell_counts = pd.crosstab(adata.obs['group'], adata.obs['scNiche'])
cell_ratios = cell_counts.div(cell_counts.sum(axis=1), axis=0)

cell_ratios.to_csv('scNiche_group3_ratios.csv')
