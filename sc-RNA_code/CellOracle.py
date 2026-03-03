import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co

adata = sc.read('rat.h5ad')
#vg=pd.read_csv("human_hvg.csv")
#adata=adata[:,vg['gene']]
#adata

sc.pp.neighbors(adata)
sc.tl.diffmap(adata)
adata.uns["iroot"] = np.flatnonzero(adata.obs["celltype"] == "Hair_follicle_stem_cell_KRT15")[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color=["dpt_pseudotime"])


#base_GRN = co.data.load_human_promoter_base_GRN()
base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata,cluster_column_name="celltype",embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_GRN)

oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

oracle.to_hdf5("cc.celloracle.oracle")

links = oracle.get_links(cluster_name_for_GRN_unit="celltype", alpha=10,verbose_level=10)

links.to_hdf5(file_path="cc.celloracle.links")
links.filter_links(p=0.05, weight="coef_abs", threshold_number=2000)


plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,save=f"{save_folder}/degree_distribution/",)
plt.rcParams["figure.figsize"] = [6, 6]

links.get_network_score()
links.merged_score.head()
links.to_hdf5(file_path="cc.celloracle.links")

import os
import sys
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co

######CLDN1 gene OE
oracle = co.load_hdf5("cc.celloracle.oracle")
links = co.load_hdf5(file_path="cc.celloracle.links")

links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, 
                              use_cluster_specific_TFdict=True)

goi = "Fos"


####KO
oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)
####OE
oracle.simulate_shift(perturb_condition={goi: 2},
                      n_propagation=3)
                      
                      
                      
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)




import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale = 25
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.savefig('quiver_plot.png',dpi=300, bbox_inches='tight')


###n_grid 和 min_mass 的参数
n_grid = 80
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('mass_thresholds_plot.png', dpi=300, bbox_inches='tight')

min_mass = 2.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('mass_filter_plot1.png', dpi=300, bbox_inches='tight')


##绘制矢量场
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
scale_simulation = 15###选择合适阈值用来调整箭头的大小,数值越小，箭头越大
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.savefig('Simulated cell identity shift vector155.png', dpi=300, bbox_inches='tight')


fig, ax = plt.subplots(figsize=[8, 8])
oracle.plot_cluster_whole(ax=ax, s=10)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
plt.savefig('scale_simulation.png', dpi=300, bbox_inches='tight')
plt.savefig('scale_simulation.pdf', bbox_inches='tight')

###比较模拟向量与开发向量
from celloracle.applications import Gradient_calculator

# Instantiate Gradient calculator object
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="dpt_pseudotime")

gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('gradient.calculate_mass_filter.png', dpi=300, bbox_inches='tight')

###伪时间数据转换为网格点
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('n_poly.png', dpi=300, bbox_inches='tight')

gradient.calculate_gradient()
scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('gradient.visualize_results.png', dpi=300, bbox_inches='tight')
plt.savefig('gradient.visualize_results.pdf',  bbox_inches='tight')

fig, ax = plt.subplots(figsize=[6, 6])
gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)
plt.savefig('gradient.plot_dev_flow_on_grid.png', dpi=300, bbox_inches='tight')

###保存结果
gradient.to_hdf5("mouse_Fosb_OE.celloracle.gradient")

##扰动分数
from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)


###调整vm使得右侧面板无法看到颜色
vm = 1###如果你在随机结果（右侧面板）中看到颜色，这意味着 vm 参数太小
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax[0])
ax[0].set_title(f"PS")
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
plt.savefig('KO.png', dpi=300, bbox_inches='tight')

###用扰动模拟矢量场显示扰动分数
scale_simulation = 15
fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
plt.savefig('FOSB_OE_scale_simulation.png',dpi=300, bbox_inches='tight')
plt.savefig('FOSB_OE_scale_simulation.pdf', bbox_inches='tight')

###绘制所有结果
dev.visualize_development_module_layout_0(s=5,
                                          scale_for_simulation=scale_simulation,
                                          s_grid=50,
                                          scale_for_pseudotime=scale_dev,
                                          vm=vm)
plt.savefig('FOS_OE_all_result.png',dpi=300, bbox_inches='tight')
plt.savefig('FOS_OE_all_result.pdf',bbox_inches='tight')