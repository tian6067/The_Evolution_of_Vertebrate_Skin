###SAMap_pycode
import pandas as pd
import scanpy as sc
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
adata = sc.read_10x_mtx("path/to/SAMaptest/YN-A")
adata.write('path/to/h5ad/yn.h5ad')
fn1="path/to/h5ad/goat_pr.h5ad"
fn2="path/to/h5ad/sheep_pr.h5ad"
fn3="path/to/h5ad/yn_pr.h5ad"

filenames = {'go':fn1,'sh':fn2,'yn':fn3}

sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)

sam3=SAM()
sam3.load_data(fn3)

sams={'go':sam1,'sh':sam2,'yn':sam3}

sm = SAMAP(sams,f_maps = 'path/to/SAMaptest/maps',)

sm.run(pairwise=True)
samap = sm.samap

keys = {'ch':'major_cell','du':'major_cell','go':'celltype'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)

os.chdir("path/to/SAMap/") #修改所需文件的路径
MappingTable.to_csv('path/to/SAMap/chdugo.txt', sep='\t', index=True)

#SAMap通过biomart数据框将蛋白ID转为基因symbol
blast_output<-read.table("D:/SAMap/maps/shyn/yn_to_sh.txt")

first_gene<-unique(blast_output$V1)##
second_gene<-unique(blast_output$V2)

library(biomaRt)
mart<- useMart("ENSEMBL_MART_ENSEMBL")
list<-listDatasets(mart) 
#a<-listAttributes(species)

#goat:chircus_gene_ensembl
#sheep:oarambouillet_gene_ensembl
#yn:btaurus_gene_ensembl
species1 <- useMart('ensembl',dataset = "btaurus_gene_ensembl") ##
ids<-getBM(attributes=c('ensembl_peptide_id_version','external_gene_name',"ensembl_gene_id"),
           filters='ensembl_peptide_id_version',
           values=first_gene,
           mart=species1)

library(dplyr)
ids <- ids %>%
  mutate(external_gene_name = if_else(external_gene_name == "", ensembl_gene_id, external_gene_name))
ids<-ids[,1:2]

colnames(ids)<-c("V1","name1")
new<-left_join(blast_output,ids)


species2 <- useMart('ensembl',dataset = "oarambouillet_gene_ensembl") ##
ids<-getBM(attributes=c('ensembl_peptide_id_version','external_gene_name',"ensembl_gene_id"),
           filters='ensembl_peptide_id_version',
           values=second_gene,
           mart=species2)

ids <- ids %>%
  mutate(external_gene_name = if_else(external_gene_name == "", ensembl_gene_id, external_gene_name))
ids<-ids[,1:2]

colnames(ids)<-c("V2","name2")
new1<-left_join(new,ids)

output<-data.frame(new1$name1,new1$name2,new1[,3:12])
write.table(output,"D:/SAMap/newmaps/shyn/new_yn_to_sh.txt",col.names = F,row.names = F,sep = "\t")


#SAMap通过GTF文件将蛋白ID转为基因symbol
library(rtracklayer) #读取gtf
library(dplyr)

SAMap_GTF<-function(path_to_input_txt,
                    species_id1,
                    source1,
                    species_id2,
                    source2,
                    path_to_species1_gtf,
                    path_to_species2_gtf,
                    output){
  
  inputwd<-path_to_input_txt
  species1_id<-species_id1
  species2_id<-species_id2
  species1_gtf<-path_to_species1_gtf
  species2_gtf<-path_to_species2_gtf
  outwd<-output
  
  blast_output<-read.table(inputwd)
  
  if (source1=='NCBI'){
    print('The first source is NCBI')
    #对于V1
    species1_gene<-data.frame(protein_id=unique(blast_output$V1))
    
    gtf1<-as.data.frame(import(species1_gtf)) 
    gtf1<-gtf1[c("gene","protein_id")]
    gtf1 <- na.omit(gtf1)
    
    
    unique_rows <- gtf1 %>%
      distinct(protein_id, .keep_all = TRUE)
    
    merged<-left_join(species1_gene,unique_rows)
    
    colnames(merged)<-c("V1","name1")
    new<-left_join(blast_output,merged)
    
  }else if (source1 == "ENSEMBL") {
    print('The first source is ENSEMBL')
    species1_gene<-data.frame(protein_id=unique(blast_output$V1))
    
    gtf1<-as.data.frame(import(species1_gtf)) 
    gtf1<-gtf1[c('gene_id',"gene_name","protein_id",'protein_version')]
    gtf1$gene_name <- ifelse(is.na(gtf1$gene_name), gtf1$gene_id, gtf1$gene_name)
    gtf1<-gtf1[c("gene_name","protein_id",'protein_version')]
    
    gtf1 <- na.omit(gtf1)
    gtf1$merged <- paste(gtf1$protein_id, gtf1$protein_version, sep = ".")
    gtf1<-gtf1[c("gene_name",'merged')]
    names(gtf1)<-c('gene_name','protein_id')
    
    unique_rows <- gtf1 %>%
      distinct(protein_id, .keep_all = TRUE)
    
    merged<-left_join(species1_gene,unique_rows)
    
    colnames(merged)<-c("V1","name1")
    new<-left_join(blast_output,merged)
  } else {
    # 可选：处理其他情况
    warning("未知的source1值:", source1)
  }
  
  if (source2=='NCBI'){
    print('The second source is NCBI')
    species2_gene<-data.frame(protein_id=unique(blast_output$V2))
    
    gtf2<-as.data.frame(import(species2_gtf)) 
    gtf2<-gtf2[c("gene","protein_id")]
    gtf2 <- na.omit(gtf2)
    
    unique_rows <- gtf2 %>%
      distinct(protein_id, .keep_all = TRUE)
    
    merged<-left_join(species2_gene,unique_rows)
    
    colnames(merged)<-c("V2","name2")
    new1<-left_join(new,merged)
  }else if (source2 == "ENSEMBL") {
    print('The second source is ENSEMBL')
    species2_gene<-data.frame(protein_id=unique(blast_output$V2))
    
    gtf2<-as.data.frame(import(species2_gtf)) 
    gtf2<-gtf2[c('gene_id',"gene_name","protein_id",'protein_version')]
    gtf2$gene_name <- ifelse(is.na(gtf2$gene_name), gtf2$gene_id, gtf2$gene_name)
    gtf2<-gtf2[c("gene_name","protein_id",'protein_version')]
    
    gtf2 <- na.omit(gtf2)
    gtf2$merged <- paste(gtf2$protein_id, gtf2$protein_version, sep = ".")
    gtf2<-gtf2[c("gene_name",'merged')]
    names(gtf2)<-c('gene_name','protein_id')
    
    unique_rows <- gtf2 %>%
      distinct(protein_id, .keep_all = TRUE)
    
    merged<-left_join(species2_gene,unique_rows)
    
    colnames(merged)<-c("V2","name2")
    new1<-left_join(new,merged)
  } else {
    # 可选：处理其他情况
    warning("未知的source2值:", source2)
  }
  
  output<-data.frame(new1$name1,new1$name2,new1[,3:12])
  
  dir.create(paste0(outwd,species1_id,species2_id))
  
  write.table(output,file = paste0(outwd,species1_id,species2_id,"/",species1_id,"_to_",species2_id,".txt"),
              col.names = F,
              row.names = F,
              sep = "\t")
}


SAMap_GTF(path_to_input_txt="D:/BT-KYFW-2023-8028-分析结果/SAMap/cw_snake/maps/shsk/sk_to_sh.txt",
          species_id1="sk",
          source1 = 'NCBI',
          species_id2="sh",
          source2= 'ENSEMBL',
          path_to_species1_gtf="D:/BT-KYFW-2023-8028-分析结果/SAMap/cw_snake/gtf/snake.gtf",
          path_to_species2_gtf="D:/BT-KYFW-2023-8028-分析结果/SAMap/cw_snake/gtf/sheep.gtf",
          output="D:/BT-KYFW-2023-8028-分析结果/SAMap/cw_snake/newmaps/")

