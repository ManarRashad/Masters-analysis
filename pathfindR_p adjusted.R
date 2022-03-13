###########################################
## pathfindR using p adjusted for RNA,PROTEOMICS,CYTOKINES##
## 25.7.2020, DONE BY MANAR        ##

#Required packages
install.packages("pathfindR")
library(pathfindR)
library(readr)

directory = getwd()
res_files = list.files(path = directory, pattern = "res.txt", full.names = FALSE)

input_pathfindR = list()
for(i in 1:length(res_files)){
  df = read_delim(res_files[i], 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
  # adjust the input for pathfindR
  names(df)[names(df) == "hits"] <- "Gene.symbol"
  if(res_files[i] %in% c("cytokine_Class_hits.res.txt","proteome_Class_hits.res.txt")){
    names(df)[names(df) == "lfc.diff"] <- "logFC"
    names(df)[names(df) == "t.pval.adj"] <- "adj.P.Val"
  }
    
  if(res_files[i] == "DESeq_class_hits_res.txt"){
    names(df)[names(df) == "log2FoldChange"] <- "logFC"
    names(df)[names(df) == "padj"] <- "adj.P.Val"
  }
    
  df = df[,c("Gene.symbol","logFC","adj.P.Val")]
  df = as.data.frame(df)
  input_pathfindR[[res_files[i]]]=df
  
}


#change the file of input_pathfindR list only 
#PIN:Genemania for cytokines, Biogrid: for RNA/protein, geneset:KEGG, adj_method:fdr, p_threshold in 
#cytokines and protein :0.1
pathways_cytokine <- run_pathfindR(input_pathfindR$cytokine_Class_hits.res.txt,
                                   gene_sets = "KEGG", # change from default ("KEGG")
                                   pin_name_path = "STRING", # change from default ("Biogrid")
                                   output = "updated25.2/cytokine_pathway2", # change output directory
                          adj_method = "fdr", p_val_threshold = 0.1)
                          
                          
write.table(cytokines_pathways, file = "cytokines_pathways.txt", sep = "\t",
            col.names = T, quote = F, row.names = F)

RNA_clustered <- cluster_enriched_terms(pathways_RNA, plot_dend = T, plot_clusters_graph = T)
enrichment_chart(RNA_clustered, plot_by_cluster = TRUE, num_bubbles = 3)
term_gene_heatmap(result_df = pathways_RNA, genes_df = input_pathfindR$DESeq_class_hits_res.txt)
term_gene_graph(result_df = pathways_RNA, use_description = TRUE, node_size = "num_genes")
UpSet_plot(result_df = pathways_RNA, genes_df = input_pathfindR$DESeq_class_hits_res.txt)# color acc to FC value
UpSet_plot(result_df = pathways_RNA)
