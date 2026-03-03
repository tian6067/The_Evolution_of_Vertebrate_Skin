# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:32:06 2024
@author: Yahui Zhang
"""
import stereo as st
import warnings
import scanpy as sc
warnings.filterwarnings('ignore')
# %config InlineBackend.figure_format = 'retina'
import os

# file_path = '/public/home/s20223040710/skin/result/'
file_path = 'D:/Research/skin/01_read_ref/'
dir = 'D:/Research/skin/00_saw/'
os.chdir(dir)

##### 01 处理chickenL单个文件 #####
sample = "chickenL"
data_path = f'{dir}{sample}/{sample}_r1/Y02349D313.{sample}_r1.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=True)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r1_cellbin.h5ad')


data_path = f'{dir}{sample}/{sample}_r2/Y02349D313.{sample}_r2.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=True)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r2_cellbin.h5ad')

##### 02 处理axolotl单个文件 #####
sample = "axolotl"
data_path = f'{dir}{sample}/{sample}_r1/Y02350P411.{sample}_r1.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=False)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r1_cellbin.h5ad')


data_path = f'{dir}{sample}/{sample}_r2/Y02350P411.{sample}_r2.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=False)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r2_cellbin.h5ad')


##### 03 处理goldfish单个文件 #####
sample = "goldfish"
data_path = f'{dir}{sample}/{sample}_r1/Y02352G613.{sample}_r1.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=True)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r1_cellbin.h5ad')


data_path = f'{dir}{sample}/{sample}_r2/Y02352G613.{sample}_r2.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=True)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r2_cellbin.h5ad')



##### 04 处理snake单个文件 #####
sample = "snake"
data_path = f'{dir}{sample}/{sample}_r1/Y02349MA14.{sample}_r1.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=False)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r1_cellbin.h5ad')

data_path = f'{dir}{sample}/{sample}_r2/Y02349MA14.{sample}_r2.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=False)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r2_cellbin.h5ad')

data_path = f'{dir}{sample}/{sample}_r3/Y02349MA14.{sample}_r3.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=False)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_r3_cellbin.h5ad')


##### 05 处理cy2单个文件 #####
sample = "cy2"
data_path = f'{dir}{sample}/{sample}/Y02352GD12.{sample}.label.cellbin.gef'
data_info = st.io.read_gef_info(data_path)  # 读取数据的meta信息
data_info
#*bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
data = st.io.read_gef(file_path=data_path, bin_type='cell_bins', gene_name_index=True)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
data
#查看细胞名，基因名
data.plt.spatial_scatter(cells_key=['dnbCount'], dot_size=None, palette='rainbow')
data.cells.cell_name, data.genes.gene_name
adata = st.io.stereo_to_anndata(data,flavor='scanpy')
adata.obs['sample'] = sample
adata.obs['orig.ident'] = sample
adata.write(dir + sample + '_cellbin.h5ad')