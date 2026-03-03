import os
import scniche as sn
from scniche.pl import *
print("Last run with scNiche version:", sn.__version__)
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors
warnings.filterwarnings('ignore')

outdir = '/public/home/s20223040710/skin/04_scNiche/'  # 根据需要设置路径
# outdir = 'D:/Research/skin/04_scNiche/'  # 根据需要设置路径
os.chdir(outdir)

# set seed
sn.pp.set_seed()
sample = 'adata_8'
for target_num in range(11, 12):  #5, 16
    print(f"\n=== k={target_num} ===")
    adata=sc.read_h5ad(f'{outdir}{sample}_k{target_num}_scNiche.h5ad')
    adata.X = adata.layers['counts']
    # normalize first
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    kwargs = {'figsize': (5, 3)}
    palette_use = ["#94ffb5", "#000000", "#8f7c00", "#dfd8a9", "#ffcc99", "#2bce48", "#993f00", "#005c31", "#0075dc", "#f0a0ff", "#c20088", "#4c005c", "#9dcc00", "#f49896", "#e0cc8c"]
    sn.pl.stacked_barplot(adata, x_axis='scNiche', y_axis='predicted_cell_type', mode='proportion', palette=palettes.default_20, kwargs=kwargs)
    plt.savefig(f'{outdir}figures/stacked_barplot_celltype_k{target_num}.pdf', bbox_inches='tight', dpi=300)
    
    # sn.pl.stacked_barplot(adata, x_axis='scNiche', y_axis='sample', mode='proportion', palette=palette_use, kwargs=kwargs)
    # plt.savefig(f'{outdir}figures/stacked_barplot_sample_k{target_num}.pdf', bbox_inches='tight', dpi=300)
    
    sn.pl.stacked_barplot(adata, x_axis='sample', y_axis='scNiche', mode='proportion', palette=palette_use, kwargs=kwargs)
    plt.savefig(f'{outdir}figures/stacked_barplot_sample_niche_k{target_num}.pdf', bbox_inches='tight', dpi=300)
    def simplify_breed(x):
        x = str(x)
        if 'h_goat' in x: return 'h_goat'
        elif 'h_sheep' in x: return 'h_sheep'
        elif 'rg_' in x: return 'rg'
        elif 's_x' in x: return 's_x'
        else: return 'unknown'
        
    adata.obs['breed'] = adata.obs['sample'].apply(simplify_breed).astype('category')
    sn.pl.stacked_barplot(adata, x_axis='breed', y_axis='scNiche', mode='proportion', palette=palette_use, kwargs=kwargs)
    plt.savefig(f'{outdir}figures/stacked_barplot_breed_niche_k{target_num}.pdf', bbox_inches='tight', dpi=300)
    
    def simplify_usage(x):
        x = str(x)
        if 'h_' in x: return 'Hair'
        elif 'rg_' in x: return 'Cashmere'
        elif 's_x' in x: return 'Wool'
        else: return 'unknown'
        
    adata.obs['usage'] = adata.obs['sample'].apply(simplify_usage).astype('category')
    sn.pl.stacked_barplot(adata, x_axis='usage', y_axis='scNiche', mode='proportion', palette=palette_use, kwargs=kwargs)
    plt.savefig(f'{outdir}figures/stacked_barplot_usage_niche_k{target_num}.pdf', bbox_inches='tight', dpi=300)
    
    # cell type enrichment
    sn.al.enrichment(adata, id_key='predicted_cell_type', val_key='scNiche', library_key='sample')
    
    # colors = ["white", "#D73027"]
    # cmap_custom = mcolors.LinearSegmentedColormap.from_list("white_to_red", colors)
    
    # plot
    kwargs = {'figsize': (11, 8), 'vmax': 4, 'cmap': 'Reds', 'linewidths': 0, 'linecolor': 'white'}
    # row_order =[ 'Niche0', 'Niche1', 'Niche2', 'Niche3', 'Niche4', 'Niche5', 'Niche6', 'Niche7',  'Niche8', 'Niche9', 'Niche10', 'Niche11', 'Niche12', 'Niche13', 'Niche14', ]
    row_order = [f'Niche{i}' for i in range(len(sorted(adata.obs['scNiche'].cat.categories)))]
    col_order = [
    "Hair_follicle_stem_cell_KRT15",
    "Hair_follicle_stem_cell_KRT15_LGR5",
    "Hair_follicle_stem_cell_LGR5",
    "Hair_Matrix",
    "Inner_root_sheath",
    "Interfollicle_epidermis",
    "Dermal_papilla",
    "Duct_cells",
    "Outer_root_sheath",
    "Proliferating_cells_basal",
    "Proliferating_cells_matrix",
    "Sebaceous_gland",
    "Sweat_cells",
    "Transit_amplifying_cells",
    "Angiogenic_EC",
    "Arteriole_EC",
    "capillary_EC",
    "LEC",
    "Post_capillary_venule",
    "Venule_EC",
    "SMC_Pericytes",
    "Mesenchymal",
    "Pro_inflammatory",
    "Secretory_papillary",
    "Secretory_reticular",
    "Langerhans_cell",
    "T_cell",
    "Macrophages",
    "Mast_cell",
    "NKT",
    "Dendritic_cells"
    ]
    sn.pl.enrichment_heatmap(adata, id_key='predicted_cell_type', val_key='scNiche', binarized=False, show_pval=True, row_order=row_order,col_order=col_order, anno_key=None, anno_palette=palette_use, kwargs=kwargs)    
    plt.savefig(f'{outdir}figures/stacked_enrichment_heatmap_k{target_num}.pdf', bbox_inches='tight', dpi=300)

#################保存到csv文件里面##############
columns_to_save = ['scNiche', 'predicted_cell_type', 'sample', 'breed', 'usage']
adata.obs[columns_to_save].to_csv(f'{outdir}{sample}_k{target_num}_obs.csv')


###############把mol的niche保存到scniche里面###########
columns_to_copy = ['leiden','molNiche']
sample = 'adata_8'
target_num = 12
adata_target=sc.read_h5ad(f'{outdir}{sample}_k{target_num}_scNiche.h5ad')

for col in columns_to_copy:
    adata_target.obs[col] = adata.obs[col]

adata_target.write(f'{sample}_res0.8_k12_Niche.h5ad')
