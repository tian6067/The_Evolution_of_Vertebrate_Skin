import scanpy as sc
import os
import pandas as pd
import seaborn as sns
import pandas as pd
import scanpy as sc
import scanpy.external as sce

outdir = 'D:/Research/skin/03_niche_NMF'  # 根据需要设置路径
os.chdir(outdir)

adata = sc.read('adata_8.h5ad')

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
adata.raw = adata.copy()
sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True)

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

sce.pp.harmony_integrate(adata, 'sample')
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=10, n_pcs=30)
# sc.pp.neighbors(adata,use_rep = 'X_pca')
# sc.external.pp.bbknn(adata, batch_key='sample')

resolutions = [0.6] #选择分辨率为0.6的

for res in resolutions:
    print(f"正在处理分辨率: {res}")
    sc.tl.leiden(adata, resolution=res, key_added="leiden")#horse_R2_zw horse_R1_zw 0.6
    sc.tl.umap(adata)
    # adata.write(f"{sample}_res{res}.h5ad")
    sc.pl.umap(adata, color=["leiden"], palette="tab20",save=f"_mol_niche_res{res}.pdf")
    
    for sample in adata.obs['sample'].unique():
        adata_sample = adata[adata.obs["sample"] == sample]
        adata_sample.obs['molNiche'] = 'Niche' + adata_sample.obs['leiden'].astype(str)
        palette = {
        ########分子Niche的颜色定义#########
            "Niche0": "#2235b5",
            "Niche1": "#6f774b",
            "Niche2": "#d10040",
            "Niche3": "#84cd5d",
            "Niche4": "#8e76db",
            "Niche5": "#71c8a5",
            "Niche6": "#f6a26d",
            "Niche7": "#924373",
            "Niche8": "#a6513c",
            "Niche9": "#000000",
            "Niche10": "#72c7ff",
            "Niche11": "#cabd00",
            "Niche12": "#f0a0ff",
            "Niche13": "#00a34d",
            "Niche14": "#c8aecb",
            "Niche15": "#de75ad",
            
        }
        sc.pl.spatial(adata_sample, color="molNiche",spot_size=25,title=f"{sample}_mol_niche", size=1.3, palette=palette, 
                      save=f"_mol_niche_{sample}_res{res}.pdf")
    
    ######################03 保存到csv#############
    leiden_clusters = adata.obs['leiden']  # 'leiden' 是聚类结果所在的列名
    cluster_info = adata.obs[['leiden']].copy()  # 只保留 leiden 聚类结果
    cluster_info['row_id'] = adata.obs_names  # 获取细胞的标识符
    cluster_info['mol_niche'] = cluster_info['leiden']
    cluster_info = cluster_info[['row_id', 'mol_niche']]
    cluster_info.to_csv(f'leiden_clusters_{res}.csv', index=False)
