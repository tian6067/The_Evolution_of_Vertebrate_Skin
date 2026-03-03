library(Biostrings)

read_gene_list <- function(file_path) {
  gene_list <- readLines(file_path)
  return(gene_list)
}

extract_gene_name <- function(header) {
  match <- regexpr("gene=([A-Za-z0-9_]+)", header)
  if (match == -1) {
    return(NULL)  
  }
  gene_name <- regmatches(header, match)
  gene_name <- sub("gene=", "", gene_name) 
  return(gene_name)
}


extract_longest_transcript <- function(fasta_file, gene_list, output_file) {
  
  sequences <- readDNAStringSet(fasta_file)

  longest_transcripts <- list()
  
  for (seq_id in names(sequences)) {
    gene_name <- extract_gene_name(seq_id)
    #gene_name <- toupper(gene_name) 
    if (!is.null(gene_name) && gene_name %in% gene_list) {
      if (!gene_name %in% names(longest_transcripts)) {
        longest_transcripts[[gene_name]] <- list()
      }
      longest_transcripts[[gene_name]][[seq_id]] <- sequences[[seq_id]]
    }
  }
  
  longest_transcripts_final <- list()
  

  for (gene_name in names(longest_transcripts)) {
    transcripts <- longest_transcripts[[gene_name]]
    

    lengths <- sapply(transcripts, length)  # 获取每个转录本的长度
    longest_id <- names(transcripts)[which.max(lengths)]  # 获取最长转录本的ID

    longest_transcripts_final[[gene_name]] <- transcripts[[longest_id]]
  }

  if (length(longest_transcripts_final) > 0) {
    writeXStringSet(DNAStringSet(longest_transcripts_final), output_file)
    cat("提取完成，结果保存在", output_file, "\n")
  } else {
    cat("没有找到匹配的基因名。\n")
  }
}

# 主程序

main <- function(x) {
  gene_list_file <- "D:/BT-KYFW-2023-8028-分析结果/re_dnds/mam/哺乳动物同源基因.txt" 
  fasta_file <- x 
  output_file <- paste0("d:/BT-KYFW-2023-8028-分析结果/re_dnds/mam/2、最长转录本cds/",basename(x)) 

  gene_list <- read_gene_list(gene_list_file)

  extract_longest_transcript(fasta_file, gene_list, output_file)
}

file<-list.files('d:/BT-KYFW-2023-8028-分析结果/re_dnds/mam/1、mouse_rat_cds/',pattern = '.fna$',full.names = T)
for (i in file) {
  main(x=i)
}


process_file <- function(file) {
  lines <- readLines(file)
  headers <- grepl("^>", lines)
  gene_names <- sub("^>", "", lines[headers])
  starts <- which(headers) + 1
  ends <- c(which(headers)[-1] - 1, length(lines))
  sequences <- sapply(seq_along(gene_names), function(i) {
    paste(lines[starts[i]:ends[i]], collapse = "")
  })
  data.frame(gene = gene_names, sequence = sequences, stringsAsFactors = FALSE)
}

file_list <- list.files('D:/BT-KYFW-2023-8028-分析结果/re_dnds/mam/2、最长转录本cds/',
                        #pattern = "\\.fa$", 
                        full.names = TRUE)

species_data <- list()
species_names <- character()
for (file in file_list) {
  df <- process_file(file)
  species_name <- sub("\\.(fa|fna)$", "", basename(file)) 
  species_data[[species_name]] <- df
  species_names <- c(species_names, species_name)
}


gene_lists <- lapply(species_data, function(x) x$gene)
common_genes <- Reduce(intersect, gene_lists)


output_dir <- "d:/BT-KYFW-2023-8028-分析结果/re_dnds/mam/3、相同基因不同物种cds/"

for (gene in common_genes) {
  output_file <- file.path(output_dir, paste0(gene, ".fa"))
  con <- file(output_file, "w")
  for (species in species_names) {
    df_species <- species_data[[species]]
    idx <- which(df_species$gene == gene)
    if (length(idx) == 0) next  # 理论上不会发生
    seq <- df_species$sequence[idx[1]]
    writeLines(paste0(">", species), con)
    writeLines(seq, con)
  }
  close(con)
}


cd /public/home/2020304010114/dnds/aves/1、common_gene
ls *.fa | sed 's/\.fa$//' > ../list
for i in $(cat ../list); do
prank -d=/public/home/2020304010114/dnds/aves/1、common_gene/${i}.fa -t=/public/home/2020304010114/dnds/aves/aves_tree.nwk -once -codon -f=paml -o=/public/home/2020304010114/dnds/aves/2、phy/${i}.out
done

cd /public/home/2020304010114/dnds/aves
for i in $(cat list)
do
sed "s/AADAT/$i/" codeml.ctl > new_codeml.ctl
codeml new_codeml.ctl
done



