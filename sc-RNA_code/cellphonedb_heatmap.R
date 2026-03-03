###network
##单个计算 seurat名称更换
library(tidygraph)
sig_p <- read.delim("./statistical_analysis_pvalues.txt", check.names = F)
ks_cpdb_sig_interCount <- function(pval_data,#cellphonedb V5输出结果，pval.txt文件
                                   significant=0.05,#互作显著性
                                   source_celltype,#ource celltype
                                   color_set=NULL,#celltype颜色设置
                                   celltype_order=NULL,#celltype颜色设置，与color set一致
                                   showInter=F,#是否展示互作数目（节点颜色和大小是否设置为互作数目，False则代表节点设置为固定大小或celltype 细胞数）
                                   celltype.size=F,#是否展示celltype数量。
                                   scRNA=NULL,#做cellphonedb的seurat obj
                                   cellanno=NULL,#celltype所在那一列metadata的列名
                                   group.size=F)#默认大小
  
{
  
  requireNamespace("reshape2")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  requireNamespace("tidygraph")
  requireNamespace("dplyr")
  requireNamespace("igraph")
  
  all_intr <- pval_data
  col_start <- which(colnames(all_intr) == "classification")
  
  intr_pairs <- all_intr$interacting_pair
  all_intr <- t(all_intr[, -c(1:col_start)])
  colnames(all_intr) <- intr_pairs
  all_count <- reshape2::melt(all_intr)
  
  all_count$significant <- all_count$value < significant
  
  count1x <- all_count %>%
    group_by(Var1) %>%
    summarise(COUNT = sum(significant)) %>%
    as.data.frame()
  
  
  tmp <- lapply(count1x[, 1], function(x) strsplit(as.character(x), "\\|"))
  tmp <- lapply(tmp, function(x) x[[1]])
  tmp <- as.data.frame(do.call(rbind, tmp))
  colnames(tmp) <- c("SOURCE", "TARGET")
  count1x <- as.data.frame(cbind(count1x, tmp))
  all_count <- count1x[, c("SOURCE", "TARGET", "COUNT")]
  
  
  
  if (any(all_count$COUNT) > 0) {
    
    count_mat <- reshape2::acast(SOURCE ~ TARGET, data = all_count, value.var = "COUNT")
    count_mat[is.na(count_mat)] <- 0
    
  }else{
    
    stop("There are no significant results using p-value of: ", significant, call. = FALSE)
    
  }
  
  # return(count_mat)
  
  
  #plot networks
  plot_mat <- as.data.frame(count_mat)
  sourecell <- which(colnames(plot_mat)==source_celltype)
  plot_mat <- plot_mat[order(plot_mat[,sourecell], decreasing = TRUE),]
  
  
  df <- data.frame(from = rep(source_celltype,nrow(plot_mat)),
                   to = rownames(plot_mat),
                   inter_num = plot_mat[,source_celltype])
  
  #设置celltype颜色
  if(is.null(color_set)){
    
    color.use <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
                   "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416")
    
    color.use <- color.use[1:nrow(plot_mat)]
    names(color.use) <- df$to
    
  }else{
    
    color.use <- color_set
    
    if(is.null(celltype_order)){
      
      names(color.use) <- df$to
      
    }else{
      
      names(color.use) <- celltype_order
    }
    
  }
  
  
  
  #nodes
  nodes <- data.frame(name = df$to)
  nodes$inter_num <- df$inter_num
  
  
  if(group.size==F){
    
    size = rep(5, nrow(plot_mat))
    nodes$size <- size
    
  }else{
    
    metadata <- cashmere@meta.data
    celltypesize = as.data.frame(table(metadata[,cellanno]))
    colnames(celltypesize) <- c('name',"size")
    nodes = merge(nodes, celltypesize, by='name', all=F)
    
    nodes <- nodes[order(nodes[,"size"], decreasing = TRUE),]
    
  }
  
  #edge
  edges <- df[c("from","to","inter_num")]
  
  #network plot
  net <- tbl_graph(nodes = nodes, edges = edges)
  
  #plot
  p=ggraph(net,layout='igraph', algorithm = 'circle') +
    geom_edge_bend(mapping = aes(edge_width = inter_num),
                   strength = 0.2,alpha = 0.8,
                   flipped =F, edge_color = "#A9AAAA",
                   n=50, show.legend = F,
                   check_overlap =T)+
    geom_edge_loop(aes(edge_width = inter_num,
                       direction = (from - 1)*360 / length(net)),
                   colour = "#A9AAAA",
                   alpha = 0.5, show.legend = F)+
    scale_edge_width_continuous(range = c(0,5))
  
  
  if(showInter==F){
    
    
    p = p+geom_node_point(aes(size=size,colour = name), show.legend = F) +
      geom_node_point(aes(size=size), show.legend = F,
                      shape=21,colour = 'black',stroke = 1.5)+
      geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                     angle=0,hjust=0, size=3) + # 设置点的注释
      scale_size_continuous(range = c(1, 15))+
      scale_color_manual(values = color.use)+
      theme_graph()+
      theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
    
  }else{
    
    if(celltype.size==T){
      
      p = p+geom_node_point(aes(size=size,colour = inter_num)) +
        geom_node_point(aes(size=size), show.legend = F,
                        shape=21,colour = 'black',stroke = 1.5)+
        geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                       angle=0,hjust=0, fontface="bold",size=3) + # 设置点的注释
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100))+
        scale_size_continuous(range = c(1, 10))+
        theme_graph()+
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      
    }else{
      
      p = p+geom_node_point(aes(size=inter_num,colour = inter_num)) +
        geom_node_point(aes(size=inter_num), show.legend = F,
                        shape=21,colour = 'black',stroke = 1.5)+
        geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                       angle=0,hjust=0, fontface="bold",size=3) + # 设置点的注释
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100))+
        scale_size_continuous(range = c(1, 10))+
        theme_graph()+
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      
    }
    
    
  }
  
  return(p)
  
  
}


a<-list()
count<-data.frame()
df<-as.character(unique(cashmere$celltypedb))
for (i in 1:length(df)) {
  a[[i]]<-ks_cpdb_sig_interCount(pval_data = sig_p,
                                 significant=0.05,
                                 source_celltype = df[[i]],
                                 showInter=T,
                                 celltype.size=T,
                                 scRNA = cashmere,
                                 cellanno = 'celltypedb',
                                 group.size = T)
  a[[i]][["data"]]$source<-df[[i]]
  a[[i]][["data"]]$x<-NULL
  a[[i]][["data"]]$y<-NULL
  a[[i]][["data"]]$size<-NULL
  a[[i]][["data"]]$.ggraph.orig_index<-NULL
  a[[i]][["data"]]$circular<-NULL
  a[[i]][["data"]]$.ggraph.index<-NULL
  colnames(a[[i]][["data"]])<-c('TARGET','count','SOURCE')
  count<-rbind(count,a[[i]][["data"]])
}

items <- unique(c(as.character(count$TARGET), as.character(count$SOURCE)))

# Initialize an empty data frame to store the sums
sums_df <- data.frame(TARGET = character(), count = numeric(), SOURCE = character())

# Loop over each item and calculate the sums
for(item in items) {
  sum_count <- sum(count$count[count$TARGET == item | count$SOURCE == item])
  sums_df <- rbind(sums_df, data.frame(TARGET = item, count = sum_count, SOURCE = "sum"))
}

sums_df <- sums_df[order(sums_df$TARGET),]

##绘图
library(tidyr)
library(CellChat)
df.net <- spread(count, TARGET, count)# 长数据变宽
rownames(df.net) <- df.net$SOURCE
df.net <- df.net[, -1]
df.net <- as.matrix(df.net)
netVisual_circle(df.net, vertex.weight = sums_df$count,
                 weight.scale = T, label.edge= F,edge.width.max = 10,arrow.size=0.4,
                 title.name = "Number of interactions")
write.csv(count,file = 'count.csv')
write.csv(sums_df,file = 'sum.csv')

##热图
library(tidyr)
count_matrix<-spread(count, TARGET, count)
rownames(count_matrix) <- count_matrix$SOURCE
df <- count_matrix$SOURCE
count_matrix <- count_matrix[, -1]
rownames(count_matrix) <- df
count_matrix <- as.matrix(count_matrix)

count_matrix_cashmere<-count_matrix
selected_columns <- c("Dermal_papilla", "Duct_cells", "Hair_follicle_stem_cell_KRT15","Hair_follicle_stem_cell_KRT15_LGR5",
                      "Hair_follicle_stem_cell_LGR5",'Hair_Matrix','Inner_root_sheath','Interfollicle_epidermis',
                      'Outer_root_sheath','Proliferating_cells_matrix','Proliferating_cells_basal','Sebaceous_gland',
                      'Sweat_cells','Transit_amplifying_cells')
count_matrix_cashmere <- count_matrix_cashmere[selected_columns,]

custom_col_order <-cell$Var1
count_matrix_cashmere <- count_matrix_cashmere[, custom_col_order]

library(pheatmap)
pheatmap(count_matrix_cashmere, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "none", 
         cluster_cols = FALSE,
         cluster_rows = FALSE, 
         fontsize_row = 10, 
         fontsize_col = 10,
         main = "", 
         treeheight_row = 0, 
         family = 'Arial',
         color = colorRampPalette(c("dodgerblue4", 'peachpuff', 'deeppink4'))(1000),
         treeheight_col = 0,
         fontsize_number = 12,
         number_format = "%.0f", 
         legend_labels = c(0, 300),
         breaks = seq(0, 300, length.out = 1001))