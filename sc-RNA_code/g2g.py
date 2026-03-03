import pandas as pd
import scanpy as sc
import pandas as pd
import anndata
import numpy as np
import seaborn as sb
import numpy as np
import warnings
import pickle
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

from genes2genes import Main
from genes2genes import ClusterUtils
from genes2genes import TimeSeriesPreprocessor
from genes2genes import PathwayAnalyser
from genes2genes import VisualUtils
from genes2genes import MyFunctions
from optbinning import ContinuousOptimalBinning

adata_ref = anndata.read_h5ad('sub_pr_Lineage1.h5ad') # Reference dataset
adata_query = anndata.read_h5ad('sub_sec_Lineage1.h5ad') # Query dataset

adata_ref.obs.columns

sc.pp.neighbors(adata_ref)
sc.pp.neighbors(adata_query)

sc.pl.umap(adata_ref, color=["celltype"])
sc.pl.umap(adata_query, color=["celltype"])

adata_ref.obs['time']=adata_ref.obs['Lineage1']
adata_query.obs['time']=adata_query.obs['Lineage1']


print(min(adata_ref.obs['time']), max(adata_ref.obs['time']))
print(min(adata_query.obs['time']), max(adata_query.obs['time']))

adata_ref.obs['time'] = TimeSeriesPreprocessor.Utils.minmax_normalise(np.asarray(adata_ref.obs['time']))
adata_query.obs['time'] = TimeSeriesPreprocessor.Utils.minmax_normalise(np.asarray(adata_query.obs['time']))

sb.kdeplot(adata_ref.obs['time'], fill=True, label='pr', color='forestgreen') 
sb.kdeplot(adata_query.obs['time'], fill=True, label='sec', color='midnightblue'); 
plt.xlabel('pseudotime'); plt.legend(); plt.show()

x = np.asarray(adata_ref.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

x = np.asarray(adata_query.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

n_bins =15

annotation_colname = 'celltype' 
joint_cmap = ['#BD8FBA','#C4C1DE','#E6BFDA','#F19C98',
              '#F8CCCE','#F9E1CD','#E5E0AD','#F3AF9A',
              '#F3C499','#8FD2E8','#E5E0AD','#C3D8B6','#B4DDD7','#A0BEE4'
]
VisualUtils.plot_celltype_barplot(adata_ref, n_bins, annotation_colname, joint_cmap)
VisualUtils.plot_celltype_barplot(adata_query, n_bins, annotation_colname, joint_cmap)


#####G2G trajectory alignment


# 读取 CSV 
df = pd.read_csv('Lineage1_gene.csv')
gene_list = pd.Index(df['Gene'].values)

deg = pd.read_csv('g2g_tf.csv')
deg = pd.Index(deg['gene'].values)
gene_list = np.intersect1d(gene_list, deg)

aligner = Main.RefQueryAligner(adata_ref, adata_query, gene_list, n_bins)
aligner.align_all_pairs()

##########Aggregate (average) cell-level alignment across all aligned genes
aligner.get_aggregate_alignment() 

def plot_alignment_path_on_given_matrix(mat, paths, cmap='viridis'):
        fig,ax = plt.subplots(1,1, figsize=(7,7))
        sb.heatmap(mat, square=True,  cmap='viridis', ax=ax, cbar=True)  
        for path in paths: 
            path_x = [p[0]+0.5 for p in path]
            path_y = [p[1]+0.5 for p in path]
            ax.plot(path_y, path_x, color='white', linewidth=6)
        plt.xlabel("Primary_Follicle",fontweight='bold')
        plt.ylabel("Secondary_Follicle",fontweight='bold')
        ax.xaxis.tick_top() # x axis on top
        ax.xaxis.set_label_position('top')
        plt.savefig('heatmap_Lineage1.pdf')

average_alignment, alignment_path =  ClusterUtils.get_cluster_average_alignments(aligner, aligner.gene_list)
mat = ClusterUtils.get_pairwise_match_count_mat(aligner, aligner.gene_list )
print('Average Alignment: ', average_alignment)
plot_alignment_path_on_given_matrix(paths = [alignment_path], mat=mat)

# 保存对象到文件
with open("aligner_Lineage1.pkl", "wb") as f:
    pickle.dump(aligner, f)

# ordered genes according to alignment similarity statistics
df = aligner.get_stat_df()  
df

VisualUtils.plot_alignmentSim_vs_l2fc(df)
earliest_match_sorted_genes_list = aligner.show_ordered_alignments() 

#######Gene-set overrepresentation analysis on the top dissimilar genes
df = aligner.get_stat_df() 
df.sort_values(['opt_alignment_cost','alignment_similarity_percentage'], ascending=[False,False])
df[np.logical_and(list(df['alignment_similarity_percentage'] <=0.5), list(df['opt_alignment_cost'] >=30)) ].Gene
topDEgenes = df[np.logical_and(list(df['alignment_similarity_percentage'] <=0.5), list(df['opt_alignment_cost'] >=30)) ].Gene
topDEgenes



# 存储为 CSV（不带索引）
df = pd.DataFrame(df)
df.to_csv("Lineage1_g2g_gene.csv", index=False)



###绘制热图
def compute_mmldist(gex1,gex2):

        μ_S = np.mean(gex1); σ_S = np.std(gex1); 
        μ_T = np.mean(gex2); σ_T = np.std(gex2); 

        ref_data = gex1 
        query_data = gex2

        I_ref_model, I_refdata_g_ref_model = MyFunctions.run_dist_compute_v3(ref_data, μ_S, σ_S) 
        I_query_model, I_querydata_g_query_model = MyFunctions.run_dist_compute_v3(query_data, μ_T, σ_T) 
        I_ref_model, I_querydata_g_ref_model = MyFunctions.run_dist_compute_v3(query_data, μ_S, σ_S) 
        I_query_model, I_refdata_g_query_model = MyFunctions.run_dist_compute_v3(ref_data, μ_T, σ_T) 

        match_encoding_len1 = I_ref_model + I_querydata_g_ref_model + I_refdata_g_ref_model
        match_encoding_len1 = match_encoding_len1/(len(query_data)+len(ref_data))
        match_encoding_len2 = I_query_model + I_refdata_g_query_model + I_querydata_g_query_model
        match_encoding_len2 = match_encoding_len2/(len(query_data)+len(ref_data))
        match_encoding_len = (match_encoding_len1 + match_encoding_len2 )/2.0 

        null = (I_ref_model + I_refdata_g_ref_model + I_query_model + I_querydata_g_query_model)/(len(query_data)+len(ref_data))
        match_compression =   match_encoding_len - null 

        return round(float(match_compression.numpy()),4) 


genes = aligner.gene_list
match_costs_across_time = {}
ranked_genes = []

for t in range(0,13):
    match_costs = []
    for g in genes:
        al = aligner.results_map[g]
        Mmat = al.fwd_DP.DP_M_matrix[1:14,1:14] 
        c = compute_mmldist(al.S.data_bins[t], al.T.data_bins[t])
        match_costs.append(c)
        #match_costs.append(Mmat[t,t])
    d = pd.DataFrame( [genes, match_costs] ) .transpose() 
    d.columns = ['gene','mml_opt_cost_in_prefix_alignment']
    d['rank'] = d.sort_values('mml_opt_cost_in_prefix_alignment', ascending=False).index
    match_costs_across_time[t] = d

gene_ordered = aligner.gene_list
mat = [] 
for t in range(0,13):
    d = match_costs_across_time[t] 
    d = d.set_index('gene')
    mat.append(list(d.loc[gene_ordered]['mml_opt_cost_in_prefix_alignment'])) 
mat = pd.DataFrame(mat).transpose()
mat.index = gene_ordered

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
t = scaler.fit_transform(mat.T)
t = pd.DataFrame(t.T)
t.index = mat.index

c = sb.heatmap(t.loc[genes])
c = sb.heatmap(np.log1p(mat.loc[genes]) )


from gsea_api.molecular_signatures_db import MolecularSignaturesDatabase
from adjustText import adjust_text
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import zscore

def visualize_gene_alignment(al_obj, vs):
    
    vega_20 = [
            '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728',
            '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2',
            '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5',
        ]

    fig = plt.figure(figsize=(4,3))
    heights = [1, 3, 1] 
    gs = plt.GridSpec(3, 1, height_ratios=heights)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0],sharex=ax1)
    ax3 = fig.add_subplot(gs[2, 0],sharex=ax1)

    #plt.subplots(3,1, gridspec_kw={'height_ratios': [1,10]})
    plt.subplot(3,1,1)
    vs.metaS.apply(lambda x: x*100/sum(x), axis=1).plot(kind='bar',stacked=True,color=vega_20, grid = False, legend=False, width=0.7, ax=ax1)
    vs.metaT.apply(lambda x: x*100/sum(x), axis=1).plot(kind='bar',stacked=True,color=vega_20, grid = False, legend=False, width=0.7,ax=ax3)
    plt.subplot(3,1,2)
    for i in range(al_obj.matched_region_DE_info.shape[0]):
                S_timebin = int(al_obj.matched_region_DE_info.iloc[i]['ref_bin'])
                T_timebin = int(al_obj.matched_region_DE_info.iloc[i]['query_bin'])
                #print(S_timebin, T_timebin)
                x_vals = [al_obj.matched_region_DE_info.iloc[i]['query_pseudotime'],al_obj.matched_region_DE_info.iloc[i]['ref_pseudotime']]#       [al_obj.matched_region_DE_info.iloc[i]['ref_pseudotime'],al_obj.matched_region_DE_info.iloc[i]['query_pseudotime']] 
                x_vals = [T_timebin+1, S_timebin+1]
                y_vals = [0,1]
                plt.plot(x_vals, y_vals, marker='.', color='black', linewidth=0.5)
  #  ax2.set_xlim([0,14])

    def set_grid_off(ax):

        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.figure.tight_layout()
        ax.grid(False)

    set_grid_off(ax1); set_grid_off(ax2); set_grid_off(ax3); 
    ax1.set_ylabel('Reference')
    ax3.set_ylabel('Query')
    fig.text(0.5, -0.05, 'Pseudotime bins with cell type composition', ha='center')
    ax1.set_title('Alignment w.r.t cell type compositions')



def plot_alignment_path_on_given_matrix(mat, paths, cmap='viridis'):

        fig,ax = plt.subplots(1,1, figsize=(7,7))
        sb.heatmap(mat, square=True,  cmap='viridis', ax=ax, cbar=True,xticklabels=True, yticklabels=True)  
        for path in paths: 
            path_x = [p[0]+0.5 for p in path]
            path_y = [p[1]+0.5 for p in path]
            ax.plot(path_y, path_x, color='white', linewidth=6)
        plt.xlabel("Reference pseudotime",fontweight='bold')
        plt.ylabel("Organoid pseudotime",fontweight='bold')
        ax.xaxis.tick_top() # x axis on top
        ax.xaxis.set_label_position('top')        
        
def get_pathway_alignment_stat(aligner, pathway_name, cluster=False):
    
    print('PATHWAY ======= ',pathway_name)
    GENE_LIST = IGS.SETS[pathway_name]
    perct_A = []
    perct_S = []
    perct_T = []
    for gene in GENE_LIST:
        series_match_percent = aligner.results_map[gene].get_series_match_percentage()
        perct_A.append(series_match_percent[0])
        perct_S.append(series_match_percent[1])
        perct_T.append(series_match_percent[2])

    print('mean matched percentage: ', round(np.mean(perct_A),2),'%' )
    print('mean matched percentage wrt ref: ',round(np.mean(perct_S),2),'%'  )
    print('mean matched percentage wrt query: ', round(np.mean(perct_T),2),'%' )
   # plt.subplots(1,1,figsize=(2,2))
    average_alignment, alignment_path =  ClusterUtils.get_cluster_average_alignments(aligner, GENE_LIST)
    mat = ClusterUtils.get_pairwise_match_count_mat(aligner,GENE_LIST )
    print('Average Alignment: ', average_alignment)
    #plot_alignment_path_on_given_matrix(paths = [alignment_path], mat=mat)


  #  plt.xlabel('Ref pseudotime')
  #  plt.ylabel('Organoid pseudotime')
   # plt.savefig('Fibro_ref_organoid_'+pathway_name+'_overall_alignment.png')
    plot_mean_trend_heatmaps(aligner, pathway_name,cluster=cluster) 

def plot_DE_genes(pathway_name):
    PATHWAY_SET = IGS.SETS[pathway]
    ax=sb.scatterplot(x['l2fc'],x['sim']*100,s=50, legend=False, hue =x['sim'] ,palette=sb.diverging_palette(15, 133, s=50, as_cmap=True),edgecolor='k',linewidth=0.3)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.ylabel('Alignment Similarity %', fontsize=12, fontweight='bold')
    plt.xlabel('L2FC mean expression', fontsize = 12, fontweight='bold')
    plt.grid(False)
    plt.tight_layout()

    TEXTS = [] 
    for label, a, b in zip(x.index, x['l2fc'],x['sim']*100):
        if(label in PATHWAY_SET):# and b<=50):
            TEXTS.append(ax.text(a, b, label, color='white', fontsize=9, fontweight='bold',bbox=dict(boxstyle='round,pad=0.1', fc='black', alpha=0.75)))
    adjust_text(TEXTS, expand_points=(2, 2),arrowprops=dict(arrowstyle="->", color='black', lw=2))
    plt.title(pathway_name,fontweight='bold', fontsize=15)
    
    
def plot_heatmaps(mat_ref,mat_query,pathway_name, cluster=False):
    
    if(cluster):
        g=sb.clustermap(mat_ref, figsize=(0.4,0.4), col_cluster=False) 
        gene_order = g.dendrogram_row.reordered_ind
        df = pd.DataFrame(g.data2d) 
        df.index = IGS.SETS[pathway_name][gene_order]
    else:
        df=mat_ref
    
    plt.subplots(1,2,figsize=FIGSIZE) #8,14/7 ******************************************************
    max_val = np.max([np.max(mat_ref),np.max(mat_query)]) 
    min_val = np.min([np.min(mat_ref),np.min(mat_query)]) 
    plt.subplot(1,2,1)
    ax=sb.heatmap(df, vmax=max_val,vmin=min_val, cbar_kws = dict(use_gridspec=False,location="top")) 
    plt.title('Ref')
    ax.yaxis.set_label_position("left")
    for tick in ax.get_yticklabels():
        tick.set_rotation(360)
    plt.subplot(1,2,2)
    if(cluster):
        mat_query = mat_query.loc[IGS.SETS[pathway_name][gene_order]] 
    ax = sb.heatmap(mat_query,vmax=max_val,  vmin=min_val,cbar_kws = dict(use_gridspec=False,location="top"), yticklabels=False) 
    plt.title('Organoid')
    plt.savefig(pathway_name+'_heatmap.pdf', bbox_inches='tight')
    plt.show()
    
    
# smoothened/interpolated mean trends + Z normalisation 
def plot_mean_trend_heatmaps(aligner, pathway_name, cluster=False):
    S_mat = []
    T_mat = []
    S_zmat = []
    T_zmat = []

    for gene in IGS.SETS[pathway_name]:

        fS = pd.DataFrame([aligner.results_map[gene].S.mean_trend, np.repeat('Ref', len(aligner.results_map[gene].S.mean_trend))]).transpose()
        fT = pd.DataFrame([aligner.results_map[gene].T.mean_trend, np.repeat('Organoid', len(aligner.results_map[gene].T.mean_trend))]).transpose()
        f = pd.concat([fS,fT])
        f[0] = np.asarray(f[0], dtype=np.float64)
        from scipy.stats import zscore
        f['z_normalised'] = zscore(f[0])
        S_mat.append(np.asarray(f[f[1]=='Ref'][0]))
        T_mat.append(np.asarray(f[f[1]=='Organoid'][0]))    
        S_zmat.append(np.asarray(f[f[1]=='Ref']['z_normalised']))
        T_zmat.append(np.asarray(f[f[1]=='Organoid']['z_normalised']))  
    S_mat = pd.DataFrame(S_mat)
    T_mat = pd.DataFrame(T_mat)
    S_zmat = pd.DataFrame(S_zmat)
    T_zmat = pd.DataFrame(T_zmat)
    
    S_mat.index = IGS.SETS[pathway_name]
    T_mat.index = IGS.SETS[pathway_name]
    S_zmat.index = IGS.SETS[pathway_name]
    T_zmat.index = IGS.SETS[pathway_name]
    
   # print('Interpolated mean trends')
   # plot_heatmaps(S_mat, T_mat, pathway_name, cluster=cluster)
    print('Z-normalised Interpolated mean trends')
    plot_heatmaps(S_zmat, T_zmat, pathway_name,cluster=cluster)
    
    
class InterestingGeneSets:
    def __init__(self):
        self.SETS = {}  # 只保留基因集字典，不加载MSigDB

    def add_new_set(self, geneset, usersetname, avail_genes=None):
        """直接添加自定义基因集，无需过滤"""
        self.SETS[usersetname] = np.asarray(geneset)
        print(f"Added gene set: {usersetname} (n={len(self.SETS[usersetname])})")


# 初始化并添加自定义基因集
IGS = InterestingGeneSets()
IGS.add_new_set(topDEgenes, 'TF_compare_Lineage1')

###绘制热图
FIGSIZE = (20, 25)  # 热图尺寸 (宽度, 高度)
get_pathway_alignment_stat(aligner, 'TF_compare_Lineage1', cluster=True)

for g in topDEgenes:
    VisualUtils.plotTimeSeries(g, aligner)
    plt.savefig('Lineage1/'+g+'_pr_vs_sec.pdf')







