import scanpy as sc
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import numpy as np
import matplotlib.image as mpimg
from matplotlib.patches import Patch
outdir = 'D:/Research/skin/02_cell2location/'  # 根据需要设置路径
os.chdir(outdir)
from plot_cell2loc import plot_spatial
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


sample = 's_x2_bin50'
sample = 's_x_33_bin50'
sample_names = ['s_x_33_bin50_5_50', 's_x_33_bin50_10_50', 's_x_33_bin50_5_100', 's_x_33_bin50_10_100',
                's_x_33_bin50_5_200', 's_x_33_bin50_10_200'] #s_x_33_bin50

sample_names = ['s_x_33_cellbin_subtype_5_20', 's_x_33_cellbin_subtype_5_50', 's_x_33_cellbin_subtype_2_50']

sample_names = ['s_x_33_cellbin_KC_subtype', 's_x2_cellbin_KC_subtype',
                's_x_33_cellbin_withoutKC_subtype', 's_x2_cellbin_withoutKC_subtype'] #s_x_33_bin50

sample_names = ["h_sheep9_cellbin_subtype", "h_sheep12_cellbin_subtype", "h_goat10_cellbin_subtype", "h_goat29_cellbin_subtype", 
           "rg_12_cellbin_subtype", "rg_10_cellbin_subtype", "s_x_33_cellbin_subtype", "s_x2_cellbin_subtype",
          ] #"goldfish_cellbin", "chickenL_cellbin", "axolotl_cellbin", "snake_cellbin"

sample_names = ["chickenL_r1_cellbin", "chickenL_r2_cellbin",]

sample_names = ["chickenL_r1_cellbin", "chickenL_r2_cellbin",
                "axolotl_r1_cellbin",  "axolotl_r2_cellbin",
                "goldfish_r1_cellbin", "goldfish_r2_cellbin", 
                "snake_r1_cellbin", "snake_r2_cellbin", "snake_r3_cellbin"]

sample_names = ["adata_8"]

for sample in sample_names:
    adata_sp = sc.read(outdir + sample + '_cell2loc_anno.h5ad')
    ####0401得分值最高的为细胞类型####
    cell_type_predictions = adata_sp.obsm['q05_cell_abundance_w_sf']
    cell_type_predictions.to_csv(outdir +sample+'_predicted_cell_types.csv', index=True)
    cell_type_predictions = pd.read_csv(outdir  + sample + '_predicted_cell_types.csv', index_col=0)
    # 提取每个细胞的最大预测概率列名
    # 使用 np.argmax() 获取每一行最大值的索引
    predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
    # 获取列名列表，即细胞类型的名称
    cell_type_names = cell_type_predictions.columns
    # 根据最大值的索引选择细胞类型的名称
    predicted_cell_types = cell_type_names[predicted_cell_types_idx]
    # 将预测的细胞类型添加到 adata_sp.obs 中
    adata_sp.obs['predicted_cell_type'] = predicted_cell_types
    adata_sp.obs['predicted_cell_type'] = adata_sp.obs['predicted_cell_type'].str.replace('^q05cell_abundance_w_sf_', '', regex=True)
    adata_sp.write(outdir + sample + '_cell2loc_anno_cell_type.h5ad')
    # 使用 sc.pl.spatial 可视化每个点的细胞类型
    palette = {"Duct_cells": "#91b7e1",
        "Sebaceous_gland": "#a6dad4",
        "Sweat_cells": "#b9d1b3",
        "Outer_root_sheath": "#dfd8a9",
        "Inner_root_sheath": "#7dd1e7",
        "Dermal_papilla": "#f2bc95",
        "Interfollicle_epidermis": "#f4a796",
        "Hair_Matrix": "#d04a53",
        "Proliferating_cells_matrix": "#f8dac5",
        "Proliferating_cells_basal": "#f9c3c6",
        "Transit_amplifying_cells": "#f49896",
        "Hair_follicle_stem_cell_KRT15_LGR5": "#e0b6d5",
        "Hair_follicle_stem_cell_LGR5": "#b58bb9",
        "Hair_follicle_stem_cell_KRT15": "#b8b8d8",
        
        "Pro_inflammatory": "#9c77ad",
        "Secretory_reticular": "#96c780",
        "Secretory_papillary": "#6f9cb7",
        "Mesenchymal": "#e6898a",

        "Venule_EC": "#66a6d1",
        "Post_capillary_venule": "#e7b7d3",
        "Arteriole_EC": "#da9d75",
        "Angiogenic_EC": "#e0cc8c",
        "capillary_EC": "#ce95be",
        "LEC": "#78b2a9",
        'SMC_Pericytes': "#f2af2b",
        
        "T_cell": "#c9e1c8",
        "NKT": "#d3c9c1",
        "Langerhans_cell": "#c0c7e2",
        "Macrophages": "#f9c8b8",
        "Dendritic_cells": "#f1d3e1",
        "Mast_cell": "#b1c6e7",
    }#颜色板
    sc.pl.spatial(adata_sp, color='predicted_cell_type',spot_size=25,palette=palette, title=f"Spatial Cell Type Distribution\n{sample}",save ='_' + sample + '_cell2loc_anno.pdf')
    #cells_list = list(palette.keys())
    cells_list = adata_sp.obs['predicted_cell_type'].unique()
    for cells in cells_list: #adata_sp.obs['predicted_cell_type'].unique():
        sc.pl.spatial(adata_sp, color='predicted_cell_type',spot_size=25,palette=palette, groups=cells,
                      title=f"Spatial Celltype\n{sample}_{cells}",
                      save ='_' + sample + '_' + cells + '_cell2loc_anno.png')
    #  groups=['Angiogenic_EC']
    #sc.pl.spatial(adata_sp, color='predicted_cell_type', spot_size=40,palette="tab20", title=f"Spatial Cell Type Distribution\nSample: {sample}")
    #plt.savefig(outdir + '02_cell2location/'+ sample + '_cell2loc_anno.png', dpi=300, bbox_inches = 'tight')
    # tab20  tab20c  Paired
    

###############
# 定义你希望单独展示的细胞类型
target_cell_types = ['Dermal_papilla', 'Fibroblast', 'SMC_Pericytes']

for sample in sample_names:
    # 读取文件
    adata_sp = sc.read(outdir + sample + '_cell2loc_anno.h5ad')
    cell_type_predictions = adata_sp.obsm['q05_cell_abundance_w_sf']
    cell_type_predictions.to_csv(outdir + sample + '_predicted_cell_types.csv', index=True)
    cell_type_predictions = pd.read_csv(outdir + sample + '_predicted_cell_types.csv', index_col=0)
    
    # 提取每个细胞的最大预测概率列名
    predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
    
    # 获取列名列表，即细胞类型的名称
    cell_type_names = cell_type_predictions.columns
    
    # 根据最大值的索引选择细胞类型的名称
    predicted_cell_types = cell_type_names[predicted_cell_types_idx]
    
    # 将预测的细胞类型添加到 adata_sp.obs 中
    adata_sp.obs['predicted_cell_type'] = predicted_cell_types
    adata_sp.obs['predicted_cell_type'] = adata_sp.obs['predicted_cell_type'].str.replace('^q05cell_abundance_w_sf_', '', regex=True)

    # 使用 sc.pl.spatial 可视化每个目标细胞类型
    for target_cell_type in target_cell_types:
        if target_cell_type in adata_sp.obs['predicted_cell_type'].unique():
            # 只对目标细胞类型进行可视化
            sc.pl.spatial(adata_sp, 
                          color='predicted_cell_type', 
                          group='Inner_root_sheath',
                          spot_size=25, 
                          title=f"Spatial Distribution of {target_cell_type} in {sample}",
                          color_map=palette.get(target_cell_type, 'tab20c'), 
                          # save=f'_{sample}_{target_cell_type}_distribution.png'
                          )




# 0402 plot in spatial coordinates 每一个图片展示一个细胞类型
# cells = list(set(adata_ref.obs['celltype']))
cells = ['Dermal_papilla', 'Fibroblast', 'SMC_Pericytes']
# 筛选出在 adata_sp.obs 中存在的细胞类型列
# cells = [cell for cell in set(adata_ref.obs['celltype']) if cell in adata_sp.obs.columns]
sc.pl.spatial(adata_sp, cmap='bwr',color=cells,ncols=4, size=1.3,img_key='hires',vmin=0, vmax='p99.2',spot_size=50,save ='_' + sample + '_cell2loc_cell_anno.pdf')


# 040301 在一个面板中可视化多种细胞类型，用plot_spatial函数
sample_names = ['s_x_33_bin50_5_50', 's_x_33_bin50_10_50', 's_x_33_bin50_5_100', 's_x_33_bin50_10_100',
                's_x_33_bin50_5_200', 's_x_33_bin50_10_200'] #s_x_33_bin50
for sample in sample_names:
    adata_sp= sc.read(outdir + sample + '_cell2loc_anno.h5ad')
    clust_col = adata_sp.obsm['q05_cell_abundance_w_sf'].columns[:7] #默认选择前7个细胞类型，可以自己设置7个细胞类型
    clust_col = clust_col.str.replace('q05cell_abundance_w_sf_', '', regex=False)
    # clust_col = ["Granular_cell","Basal_cell", "Interfollicle_epidermis", "Transit_amplifying_cells" ,"Proliferating_cells_basal"]
    clust_labels = clust_col
    if 'spatial' in adata_sp.obsm:
        with plt.rc_context({"figure.figsize": (15, 15)}):
            # 绘制空间图
            fig = plot_spatial(
                adata=adata_sp,
                color=clust_col,
                labels=clust_labels,
                max_color_quantile=0.992,
                circle_diameter=3,
                show_img=False,  # 禁用图像显示
                colorbar_position="right",
                colorbar_shape={"horizontal_gaps": 0.2},
            )
            plt.savefig(outdir + 'figures/'+ sample + "_cell2loc_celltype_anno.png", dpi=300)  # 保存为 PNG 文件，dpi=300 用于高质量保存
    else:
        print("缺少空间坐标或图像数据。")



# 040302 在一个面板中可视化一种细胞类型，用plot_spatial函数
sample_names = ["sheep_R1_bw", "sheep_R1_lw", "sheep_R1_ww", "sheep_R2_lw", "sheep_R2_bw", "sheep_R2_ww"]  # 样本名称
clust_col = ["Basal cell","Capillary endothelial","Proliferative T cell","Spinous cell","Marcphages","Fibroblast",  'Smooth muscle cell']
for sample in sample_names:
    adata_sp= sc.read("D:/Research/Rumen/11_stereo_seq/02_cell2location/" + sample + "_cell2loc_anno.h5ad")
    if 'spatial' not in adata_sp.obsm:
        print(f"{sample} 缺少空间数据，跳过")
        continue
    for ct in clust_col:
        fig = plot_spatial(
            adata=adata_sp,
            color=[ct],                    # 单独传递一个列名
            labels=[ct],                    # 或者 labels=[ct]
            max_color_quantile=0.992,
            circle_diameter=3,
            show_img=False,
            colorbar_position="right",
            colorbar_shape={"horizontal_gaps": 0.2},
        )
        plt.savefig(f"{outdir}{sample}_{ct.replace(' ','_')}_cell2loc_single_celltype_anno.pdf", dpi=300, bbox_inches='tight')



##############循环绘制饼图#################
import seaborn as sns
# sample_names = ['s_x_33_bin50', 's_x2_bin50']  # 样本名称

color_dict = palette
# 循环绘图
for sample in sample_names:
    print(f"Processing {sample}...")
    # cell_type_predictions = pd.read_csv(f"{sample}_predicted_cell_types.csv", index_col=0)

    adata_sp = sc.read(f"{sample}_cell2loc_anno.h5ad")
    cell_type_predictions = adata_sp.obsm['q05_cell_abundance_w_sf']
    coords = adata_sp.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]
    # 计算图像比例
    x_span = x.max() - x.min()
    y_span = y.max() - y.min()
    aspect_ratio = y_span / x_span if x_span > 0 else 1.0
    fig, ax = plt.subplots(figsize=(12, 10 * aspect_ratio))
    
    # 自动获取所有细胞类型并创建颜色映射
    all_cell_types = []
    for i, spot in enumerate(cell_type_predictions.index):
        proportions = cell_type_predictions.loc[spot]
        proportions = proportions[proportions > 0]
        all_cell_types.extend(proportions.index.tolist())
    unique_cell_types = sorted(set(all_cell_types))
    color_palette = sns.color_palette("tab20", n_colors=len(unique_cell_types))
    color_dict = {ct: color_palette[i] for i, ct in enumerate(unique_cell_types)}
    #    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # color_dict = {ct: color_cycle[i % len(color_cycle)] for i, ct in enumerate(unique_cell_types)}
    
    
    #####绘制饼图
    for i, spot in enumerate(cell_type_predictions.index):
        proportions = cell_type_predictions.loc[spot]
        proportions = proportions[proportions > 0]
        # 只保留你手动定义了颜色的类型
        proportions = proportions[proportions.index.isin(color_dict.keys())]
        if proportions.sum() == 0:
            continue
        # 归一化 + 最多保留 4 类
        # proportions = proportions.sort_values(ascending=False).head(4) # 最多保留4类
        proportions = proportions.sort_values(ascending=False) # 保留所有的细胞类型
        proportions = proportions / proportions.sum()
        spot_colors = [color_dict[ct] for ct in proportions.index]
        axins = inset_axes(ax, width=0.03, height=0.03, loc='center',
                           bbox_to_anchor=(x[i], y[i]),
                           bbox_transform=ax.transData,
                           borderpad=0)
        axins.pie(proportions, labels=None, colors=spot_colors)
        axins.set_aspect('equal')
        axins.axis('off')
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.max(), y.min())
    ax.set_aspect('equal')
    ax.axis('off')
    # 图例（只显示你手动设定的类型）
    legend_elements = [Patch(facecolor=color_dict[ct], label=ct) for ct in color_dict.keys()]
    plt.subplots_adjust(right=0.82)
    ax.legend(handles=legend_elements, loc='center left',
              bbox_to_anchor=(1.02, 0.5), title="Cell Types")
    plt.savefig(f"{sample}_cell2location_piecharts.pdf", format="pdf", bbox_inches="tight")
    plt.close()
print("所有样本绘图完成。")



#################单个样本绘制饼图##################
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

adata_sp = sc.read(f"{sample}_cell2loc_anno.h5ad")
cell_type_predictions = adata_sp.obsm['q05_cell_abundance_w_sf']
   
coords = adata_sp.obsm['spatial']
x, y = coords[:, 0], coords[:, 1]

# 计算图像比例
x_span = x.max() - x.min()
y_span = y.max() - y.min()
aspect_ratio = y_span / x_span if x_span > 0 else 1.0

fig, ax = plt.subplots(figsize=(12, 10 * aspect_ratio))

# 自动获取所有细胞类型并创建颜色映射
all_cell_types = []
for i, spot in enumerate(cell_type_predictions.index):
    proportions = cell_type_predictions.loc[spot]
    proportions = proportions[proportions > 0]
    all_cell_types.extend(proportions.index.tolist())
unique_cell_types = sorted(set(all_cell_types))
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
color_dict = {ct: color_cycle[i % len(color_cycle)] for i, ct in enumerate(unique_cell_types)}

for i, spot in enumerate(cell_type_predictions.index):
    proportions = cell_type_predictions.loc[spot]
    proportions = proportions[proportions > 0]

    if proportions.sum() == 0:
        continue

    # 归一化,计算比例
    # proportions = proportions.sort_values(ascending=False).head(4) # 最多保留 4 类
    proportions = proportions.sort_values(ascending=False)
    
    proportions = proportions / proportions.sum()
    spot_colors = [color_dict[ct] for ct in proportions.index]

    axins = inset_axes(ax, width=0.06, height=0.06, loc='center',
                       bbox_to_anchor=(x[i], y[i]),
                       bbox_transform=ax.transData,
                       borderpad=0)
    axins.pie(proportions, labels=None, colors=spot_colors)
    axins.set_aspect('equal')
    axins.axis('off')

# 设置主图
ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.max(), y.min())
ax.set_aspect('equal')
ax.axis('off')

# 图例
legend_elements = [Patch(facecolor=color_dict[ct], label=ct) for ct in unique_cell_types]
plt.subplots_adjust(right=0.82)
ax.legend(handles=legend_elements, loc='center left',
          bbox_to_anchor=(1.02, 0.5), title="Cell Types", fontsize=8)

plt.savefig(f"{sample}_cell2location_piecharts.pdf", format="pdf", bbox_inches="tight")
plt.close()


######替换成9大细胞类型##########
cell_type_proportions = []
for sample in sample_names:
    adata_sp= sc.read("D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
    adata_sp.layers['counts'] = adata_sp.X.copy()
    sc.pp.normalize_total(adata_sp, target_sum=1e4)
    # sc.pp.log1p(adata_sp)
    # adata_sp.layers["logcounts"] = adata_sp.X.copy()
    # 读取文件
    cell_type_predictions = pd.read_csv(outdir  + sample + '_predicted_cell_types.csv', index_col=0)
    # 提取每个细胞的最大预测概率列名
    # 使用 np.argmax() 获取每一行最大值的索引
    predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
    # 获取列名列表，即细胞类型的名称
    cell_type_names = cell_type_predictions.columns
    # 根据最大值的索引选择细胞类型的名称
    predicted_cell_types = cell_type_names[predicted_cell_types_idx]
    # 将预测的细胞类型添加到 adata_sp.obs 中
    adata_sp.obs['predicted_cell_type'] = predicted_cell_types
    adata_sp.obs['predicted_cell_type'] = adata_sp.obs['predicted_cell_type'].str.replace('^q05cell_abundance_w_sf_', '', regex=True)

    # 计算每个细胞类型的比例
    cell_type_counts = adata_sp.obs['subcell'].value_counts(normalize=True)
    # 将样本名称与细胞类型比例合并到一个 DataFrame
    sample_proportion = cell_type_counts.reset_index()
    sample_proportion.columns = ['Cell_Type', sample]
    
    # 将每个样本的比例合并到一个列表中
    cell_type_proportions.append(sample_proportion)

