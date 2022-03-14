#complexheatmap_visualisation for common pathways between omics data
library(ComplexHeatmap)
library(circlize)
library(readr)
library(RColorBrewer)

final_common_FDR <- read_delim("update_commonpathways_26.2.2021/final_common_FDR_2.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

final_common_FDR[,(3:dim(final_common_FDR)[2])] <- -10*log10(final_common_FDR[,(3:dim(final_common_FDR)[2])])
# adjust colors of pathways
coul <- c("#66C2A5","#FC8D62","#8DA0CB","#e78a9e" )
hn = HeatmapAnnotation(Pathway_type = final_common_FDR$`kegg pathway map` , annotation_name_side = "left", col = list(Pathway_type= c("Others" = "#66C2A5", "Amino acid metabolism" = "#FC8D62",
                                                                                                                                     "Carbohydrate metabolism"= "#e78a9e" , "Lipid metabolism" = "#8DA0CB")))
#Pathway_type = final_common_FDR$`kegg pathway map`
final_common_FDR[,2] = NULL
x = final_common_FDR$Pathway
final_common_FDR[,1] = NULL
rownames(final_common_FDR) = x
final_common_FDR = as.matrix(final_common_FDR)
mat = t(final_common_FDR)

#adjust colors of omics features
f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "yellow", "red"))
Heatmap(mat, rect_gp = gpar(col = "grey", lwd = 2),
        column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontface = "bold", fontsize= 9), 
        heatmap_legend_param = list(title = "-10log10p-adj", legend_height = unit(4, "cm")), 
        col = f1, width = ncol(mat)*unit(0.5, "cm"), 
        height = nrow(mat)*unit(1, "cm"), clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean", clustering_method_rows = "complete",
        clustering_method_columns = "complete", top_annotation = hn)
