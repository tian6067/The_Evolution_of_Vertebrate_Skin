library(biomaRt)
library(dplyr)
library(Seurat)
library(patchwork)  
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(openxlsx)


listEnsemblArchives() 

mart_dec2021<- useMart("ENSEMBL_MART_ENSEMBL",host = "https://dec2021.archive.ensembl.org")
mart_oct2024<- useMart("ENSEMBL_MART_ENSEMBL",host = "https://oct2024.archive.ensembl.org")

list_dec2021<-listDatasets(mart_dec2021)  #查找数据库中物种信息
list_oct2024<-listDatasets(mart_oct2024) 

human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",
                 host = "https://dec2021.archive.ensembl.org/")

species <- useMart('ensembl',dataset = "chircus_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/") 

tmp<-"D:/BT-KYFW-2023-8028-分析结果/new_rds/3、去完双细胞后rds/goat-a"  ###需要修改

sce<-readRDS(tmp) 
sce.matrix<-as.data.frame(sce@assays$RNA$counts)
rm(sce)
species_genes<-rownames(sce.matrix)
#str_subset(rownames(sce.matrix),pattern="S100A12")

homology_name <- getLDS(attributes = c("external_gene_name"),
                        filters = "external_gene_name", values = species_genes,
                        mart = species,
                        attributesL = c("external_gene_name"),
                        martL = human,
                        uniqueRows = TRUE)

homology_ID <- getLDS(attributes = c("ensembl_gene_id"),
                      filters = "ensembl_gene_id", values = species_genes,
                      mart = species,
                      attributesL = c("external_gene_name"),
                      martL = human,
                      uniqueRows = TRUE)

colnames(homology_name) <- c("species", "human")
colnames(homology_ID) <- c("species", "human")
homology<-rbind(homology_name,homology_ID)
rm(homology_ID)
rm(homology_name)

homology <- homology %>%
  group_by(species) %>%
  filter(n() == 1)

homology <- homology %>%
  group_by(human) %>%
  filter(n() == 1)

species_genes<-data.frame(species=rownames(sce.matrix))

matrix<-left_join(species_genes,homology)

matrix <- matrix %>%
  mutate(human = ifelse(is.na(human), species, human))

duplicate_elements <- unique(matrix$human
                             [duplicated(matrix$human) | 
                                 duplicated(matrix$human, fromLast = TRUE)])

homology <- homology[!homology$human %in% duplicate_elements, ]


sce.matrix$species<-rownames(sce.matrix) 

sce.matrix<-left_join(sce.matrix,homology)

sce.matrix <- sce.matrix %>%
  mutate(human = ifelse(is.na(human), species, human))

rownames(sce.matrix)<-sce.matrix$human 
sce.matrix<- subset(sce.matrix, select = -c(species,human))

sce<-CreateSeuratObject(counts = sce.matrix  , 
                        min.cells = 0,
                        min.features = 0,
                        project = basename(tmp))

saveRDS(sce,paste0("D:/BT-KYFW-2023-8028-分析结果/new_rds/4、同源转换后rds/",basename(tmp),".rds")) 
write.xlsx(homology, paste0("D:/BT-KYFW-2023-8028-分析结果/同源基因table/",basename(tmp),".xlsx"))
basename(tmp)
length(homology$species)
length(sce.matrix$human)