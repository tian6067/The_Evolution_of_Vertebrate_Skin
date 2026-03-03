#####R code#####
sce<-readRDS("seurat")
input.matrix<-as.data.frame(t(GetAssayData(sce,layer="count"))) #转置
write.csv(input.matrix,"C:/Users/周晨/Desktop/camel.csv")


#####python code #####
import loompy as lp;
import numpy as np;
import scanpy as sc;

x=sc.read_csv("camel.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
animals = ['Buffalo', 'Camel', 'Cat','Cattle','Deer', 'Dog',  'Donkey'  ]
for i in animals:
  counts_file_path = '/public/home/2020304010114/A_SKIN/AAA_data_all/cellphonedb/count/'+i+'_cellphonedb_count.txt'
meta_file_path = '/public/home/2020304010114/A_SKIN/AAA_data_all/cellphonedb/meta/'+i+'_cellphonedb_meta.txt'
out_path = '/public/home/2020304010114/A_SKIN/AAA_data_all/cellphonedb/out/'+i
cpdb_file_path = '/public/home/2020304010114/A_SKIN/AAA_data_all/cellphonedb/cellphonedb_cellchat.zip'
cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_file_path,
  meta_file_path =  meta_file_path,
  counts_file_path = counts_file_path,
  counts_data = 'hgnc_symbol',
  output_path = out_path,threads = 6)