# conda activate stereo
import os
import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import pandas as pd
from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel

indir = 'D:/Research/skin/02_cell2location/'  # 根据需要设置路径
# indir = '/public/home/s20223040710/skin/02_cell2location/'
# outdir = '/public/home/s20223040710/skin/05_colocalization/'  # 根据需要设置路径
outdir = 'D:/Research/skin/05_colocalization/'  # 根据需要设置路径
os.chdir(outdir)


sample_names =  ["h_sheep9_cellbin_subtype", "h_sheep12_cellbin_subtype", "h_goat10_cellbin_subtype", "h_goat29_cellbin_subtype", 
           "rg_12_cellbin_subtype", "rg_10_cellbin_subtype", "s_x_33_cellbin_subtype", "s_x2_cellbin_subtype"]

for sample in sample_names:
    adata = sc.read(indir + sample + '_cell2loc_anno.h5ad')
    
    ####01 细胞类型和细胞类型之间的共定位
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()
    
    comps = li.ut.obsm_to_adata(adata, 'q05_cell_abundance_w_sf')
    #### 01 cell cell 不同细胞类型的共定位
    misty = genericMistyData(intra=comps, extra=comps, cutoff=0.05, bandwidth=200, n_neighs=6) #不同的平台是不一样的，分析临近几个点呢 200是通用的
    misty(model=RandomForestModel, n_jobs=10, verbose = True)
    #view= juxta 是点间，主图放了点间的 可以抽取自己想看的细胞类型，在q05_cell_abundance_w_sf里面就抽取，
    #还有 intra 和  和 juxta 和 para
    view='intra'
    fig = li.pl.interactions(misty, view=view, return_fig=True,figure_size=(10,10)) + p9.scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + p9.ggtitle(f"Cell Type Colocalization  {sample}") +p9.theme(plot_title=p9.element_text(size=12))
    fig.data.to_csv(outdir + '/' + sample + '_misty_data_'+ view +'.csv', index=True)  # index=False 可以避免保存行号
    fig.save(outdir + '/' + sample + '.spatial.misty.celltype.'+ view +'.colocalization.pdf',bbox_inches = 'tight',dpi = 500)
    view='juxta'
    fig = li.pl.interactions(misty, view=view, return_fig=True,figure_size=(10,10)) + p9.scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + p9.ggtitle(f"Cell Type Colocalization  {sample}") +p9.theme(plot_title=p9.element_text(size=12))
    fig.data.to_csv(outdir + '/' + sample + '_misty_data_'+ view +'.csv', index=True)  # index=False 可以避免保存行号
    fig.save(outdir + '/' + sample + '.spatial.misty.celltype.'+ view +'.colocalization.pdf',bbox_inches = 'tight',dpi = 500)
    view='para'
    fig = li.pl.interactions(misty, view=view, return_fig=True,figure_size=(10,10)) + p9.scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + p9.ggtitle(f"Cell Type Colocalization  {sample}") +p9.theme(plot_title=p9.element_text(size=12))
    fig.data.to_csv(outdir + '/' + sample + '_misty_data_'+ view +'.csv', index=True)  # index=False 可以避免保存行号
    fig.save(outdir + '/' + sample + '.spatial.misty.celltype.'+ view +'.colocalization.pdf',bbox_inches = 'tight',dpi = 500)

