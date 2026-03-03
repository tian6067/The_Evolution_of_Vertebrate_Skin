library(rtracklayer) 
library(seqinr) 
library(stringr)
library(openxlsx
        
##step1## 

awk 'BEGIN{FS=OFS="\t"}{if($3=="gene"||$3=="mRNA")print $1,$3,$4,$5,$7,$9}' duck.gff |awk '{ n=split($NF, arr, ";"); $NF=""; for (i=1; i<NF; i++) printf "%s\t", $i; print arr[1] "\t" arr[2] }' FS='\t' OFS='\t' > input_file

python /public/home/2020304010114/orthofinder/gff_trans_finder.py input_file    


###step件### 

gtf <- import('D:/BT-KYFW-2023-8028-分析结果/Orthofinder/pigeon/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.gtf') ###gtf
df_gtf<-as.data.frame(gtf) 

a<-read.table("C:/Users/86137/Desktop/orthofinder_input/pigeon/longest_gene_trans.txt")
b<-as.data.frame(str_split_fixed(a$V2, "-", 2)[,2])
names(b)<-"transcript_id"

df.merged<-merge(b,df_gtf)

df2<-na.omit(unique(df.merged$protein_id))

fa <- read.fasta("C:/Users/86137/Desktop/orthofinder_input/pigeon/GCF_036013475.1_bColLiv1.pat.W.v2_protein.faa"
                 ,forceDNAtolower = F) 

fl.fa<-fa[names(fa)%in%df2]

write.fasta(sequences = fl.fa,
            names = names(fl.fa),
            nbchar = 60,
            file.out = 'C:/Users/86137/Desktop/orthofinder_input/pigeon/fl.pigeon.fa') ###


###step3###
orthofinder -f . -S blast

###step4###
orthogroup<-read.table("D:/BT-KYFW-2023-8028-分析结果/Orthofinder/R_LW/Orthogroups.tsv", header=T, sep="\t")  ###

singlecopy<-read.table("D:/BT-KYFW-2023-8028-分析结果/Orthofinder/R_LW/Orthogroups_SingleCopyOrthologues.txt") ###
names(singlecopy)<-"Orthogroup"

singlecopy_orthogroup<-merge(singlecopy,orthogroup)  

human_gtf<-import("D:/BT-KYFW-2023-8028-分析结果/Orthofinder/human/GCF_000001405.40_GRCh38.p14_genomic.gtf")

human<-as.data.frame(singlecopy_orthogroup[,"fl.human"])
names(human)<-"protein_id"

merge_human<-merge(human,human_gtf,sort=F)

gene_human<-as.data.frame(unique(merge_human$gene))
names(gene_human)<-"human_gene"

singlecopy_orthogroup<-cbind(singlecopy_orthogroup,gene_human)

#table(is.na(merge_human$gene_id)) 

other<-as.data.frame(singlecopy_orthogroup[,"fl.R_LW"])  ###
names(other)<-"protein_id"

gtf<-import('D:/BT-KYFW-2023-8028-分析结果/Orthofinder/R_LW/GCF_905171775.1_aRanTem1.1_genomic.gtf')
merge_other<-merge(other,gtf,sort=F)

gene_other<-as.data.frame(unique(merge_other$gene))

names(gene_other)<-"other_gene"
singlecopy_orthogroup<-cbind(singlecopy_orthogroup,gene_other)

singlecopy_orthogroup<-singlecopy_orthogroup[,c(4,5)]
colnames(singlecopy_orthogroup)[1]<-"human"
colnames(singlecopy_orthogroup)[2]<-"species" 

homology=singlecopy_orthogroup
#
homology <- homology %>%
  group_by(species) %>%
  filter(n() == 1)

homology <- homology %>%
  group_by(human) %>%
  filter(n() == 1)

tmp="D:/BT-KYFW-2023-8028-分析结果/new_rds/3、去完双细胞后rds/R-LW" 

sce.data<-readRDS(tmp) ###
sce.matrix<-as.data.frame(sce.data@assays$RNA$counts)

species_genes<-data.frame(species=rownames(sce.matrix))

matrix<-left_join(species_genes,homology)

matrix <- matrix %>%
  mutate(human = ifelse(is.na(human), species, human))

duplicate_elements <- unique(matrix$human
                             [duplicated(matrix$human) | 
                                 duplicated(matrix$human, fromLast = TRUE)])

homology <- homology[!homology$human %in% duplicate_elements, ]

#homology <- homology[homology$human != "RAC2", ]

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

print(basename(tmp))
print(paste0("homo gene",length(homology$human)))
print(paste0("total gene",length(rownames(sce))))
a<-data.frame(gene=rownames(sce))
write.table(a,paste0("D:/BT-KYFW-2023-8028-分析结果/各个样本基因名/",basename(tmp),".txt"))
saveRDS(sce,paste0("D:/BT-KYFW-2023-8028-分析结果/new_rds/4、同源转换后rds/",basename(tmp),".rds")) 





























repeats <- singlecopy_orthogroup %>%
  group_by(species) %>%
  filter(n() > 1) %>%
  filter(species == human)

non_repeats <- singlecopy_orthogroup %>%
  group_by(species) %>%
  filter(n() == 1)

singlecopy_orthogroup <- bind_rows(repeats, non_repeats)


sce.data<-readRDS("D:/BT-KYFW-2023-8028-分析结果/new_rds/3、去完双细胞后rds/duck-b") ###
sce.matrix<-as.data.frame(sce.data@assays$RNA$counts)
a<-data.frame(species=rownames(sce.data))
a<-left_join(a,singlecopy_orthogroup)
sum(!is.na(a$human))
df_non_na <- a[!is.na(a$human), ]
write.xlsx(df_non_na, "D:/BT-KYFW-2023-8028-分析结果/同源基因table/duck-b.xlsx")###

rm(sce.data)
sce.matrix$species<-rownames(sce.matrix) 

sce.matrix<-left_join(sce.matrix,singlecopy_orthogroup)


sce.matrix <- sce.matrix %>%
  mutate(human = ifelse(is.na(human), species, human))

rownames(sce.matrix)<-sce.matrix$human #会出现重复列名(没有的话就谢天谢地)


#如果出现重复
df<- subset(sce.matrix, select = -c(species))

same <- df %>%
  group_by(human) %>%
  filter(n() > 1) %>%
  ungroup()

unsame <- df %>%
  group_by(human) %>%
  filter(n() == 1) %>%
  ungroup()


df1 <- same %>%
  group_by(human) %>%
  summarise_all(sum)

sce.matrix<-rbind(df1,unsame)

input.matrix<- sce.matrix[,2:length(colnames(sce.matrix))]

rownames(input.matrix)<-sce.matrix$human


sce<-CreateSeuratObject(counts = input.matrix  , 
                        min.cells = 0,
                        min.features = 0,
                        project = "LV_H")

saveRDS(sce,"D:/BT-KYFW-2023-8028-分析结果/new_rds/4、同源转换后rds/LV_H") ###



