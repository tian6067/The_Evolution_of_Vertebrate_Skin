library(Seurat)
sce<-readRDS("path/to/.rds")
input.matrix<-as.data.frame(t(GetAssayData(sce,layer="count")))
write.csv(input.matrix,"C:/Users/Desktop/test.csv")

#change.py
import os,sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;

x=sc.read_csv("test.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);


#run pySCENIC
dir=/public/home/2020304010114/pyscenic/database
tfs=$dir/allTFs_hg38.txt
feather=$dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
input_loom=./sample.loom
ls $tfs  $feather  $tbl  

pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs 

pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20  \
--mask_dropouts

pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10