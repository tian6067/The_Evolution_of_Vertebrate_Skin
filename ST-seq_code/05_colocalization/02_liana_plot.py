# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from functools import reduce
####################01 单细胞共定位分析作图 R1和R2求均值############
# outdir = '/public/home/s20223040710/skin/05_colocalization/'  # 根据需要设置路径
outdir = 'D:/Research/skin/05_colocalization/'  # 根据需要设置路径
os.chdir(outdir)

# 样本列表和分组方式
groupings = {
    "four": {
        "h_sheep": ['h_sheep9_cellbin_subtype', 'h_sheep12_cellbin_subtype'],  # 对应 ['h_sheep9_cellbin_subtype', 'h_sheep12_cellbin_subtype']
        "h_goat": ['h_goat10_cellbin_subtype', 'h_goat29_cellbin_subtype'],    # 对应 ['h_goat10_cellbin_subtype', 'h_goat29_cellbin_subtype']
        "rg": ['rg_12_cellbin_subtype', 'rg_10_cellbin_subtype'],              # 对应 ['rg_12_cellbin_subtype', 'rg_10_cellbin_subtype']
        "s_x": ['s_x_33_cellbin_subtype', 's_x2_cellbin_subtype']              # 对应 ['s_x_33_cellbin_subtype', 's_x2_cellbin_subtype']
    },
    "three": {
        "h_sheep_goat": ['h_sheep9_cellbin_subtype', 'h_sheep12_cellbin_subtype', 'h_goat10_cellbin_subtype', 'h_goat29_cellbin_subtype'],  # 对应 ['h_sheep9_cellbin_subtype', 'h_sheep12_cellbin_subtype', 'h_goat10_cellbin_subtype', 'h_goat29_cellbin_subtype']
        "rg": ['rg_12_cellbin_subtype', 'rg_10_cellbin_subtype'],              # 对应 ['rg_12_cellbin_subtype', 'rg_10_cellbin_subtype']
        "s_x": ['s_x_33_cellbin_subtype', 's_x2_cellbin_subtype']              # 对应 ['s_x_33_cellbin_subtype', 's_x2_cellbin_subtype']
    },
    "one": {
        "all_samples": ['h_sheep9_cellbin_subtype', 'h_sheep12_cellbin_subtype', 'h_goat10_cellbin_subtype', 'h_goat29_cellbin_subtype', 
                        'rg_12_cellbin_subtype', 'rg_10_cellbin_subtype', 's_x_33_cellbin_subtype', 's_x2_cellbin_subtype']  # 对应所有样本 ['h_sheep9_cellbin_subtype', ..., 's_x2_cellbin_subtype']
    }
}

views = ["intra", "juxta", "para"]
# view = "para"  # intra / juxta / para

for view in views:
    for group_name, groups in groupings.items():
        group_pivots = []
        all_values = []
        for breed_name, breeds in groups.items():
            # 读取并合并组内数据
            group_data = []
            for sample in breeds:
                df = pd.read_csv(outdir + sample + '_misty_data_' + view + '.csv', index_col=0)
                df = df.rename(columns={'importances': f'importances_{sample}'})
                group_data.append(df)
            
            # 合并组内样本
            merged = reduce(lambda left, right: pd.merge(left, right, on=['target', 'predictor', 'view'], how='outer'), group_data)
                
            # 计算组内平均
            importance_cols = [col for col in merged.columns if 'importances' in col]
            merged['mean_importance'] = merged[importance_cols].mean(axis=1)
            
            # 创建透视表
            pivot = merged.pivot_table(index='target', columns='predictor', values='mean_importance')
            group_pivots.append(pivot)
            
            values = pivot.values.flatten()
            valid_values = values[~np.isnan(values)]  # 移除NaN值
            all_values.extend(valid_values)
            vmin, vmax = np.min(all_values), np.max(all_values)
            
            #给行和列进行排序
            order = [
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
            pivot = pivot.reindex(index=order[::-1], columns=order)
            plt.figure(figsize=(10, 12))
            sns.heatmap(pivot, annot=False, cmap='Greens', 
               vmin=vmin, vmax=vmax,
               cbar_kws={'label': 'Mean Importance', 'shrink': 0.4, 'orientation': 'horizontal'}, 
               linewidths=0)
            plt.title(f'Colocalization - {breed_name} - {view}  ')
            plt.xlabel('Predictor')
            plt.ylabel('Target')
            plt.xticks(rotation=90)
            plt.tight_layout()        
            plt.savefig(f'figures/{view}/{group_name}/colocalization_{breed_name}_{view}_cluster.pdf', format='pdf')
            plt.show()
        
    
