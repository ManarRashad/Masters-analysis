library(circlize)
library(tidyverse)
#basic
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df

chordDiagram(mat)
circos.clear()

#ordering

par(mfrow = c(1, 2))
chordDiagram(mat, order = c("S2", "S1", "S3", "E4", "E1", "E5", "E2", "E6", "E3"))
circos.clear() # msh far2a tt7t hna aw b3d 2li gya

chordDiagram(mat, order = c("S2", "S1", "E4", "E1", "S3", "E5", "E2", "E6", "E3"))
circos.clear()

# gaps between sectors ???????????
#5 points between rows and 15 between rows and columns then 5 between columns
circos.par(gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15)) 
chordDiagram(mat)
circos.clear()

#??????????????????
circos.par(gap.after = c("S1" = 5, "S2" = 5, "S3" = 15, "E1" = 5, "E2" = 5,
                         "E3" = 5, "E4" = 5, "E5" = 5, "E6" = 15))
chordDiagram(mat)
circos.clear()

###### to gap between rows and columns ???? to gap between features and pathways
chordDiagram(mat, big.gap = 30)

######### position of gap ???? needed 
par(mfrow = c(1, 2))
circos.par(start.degree = 85, clock.wise = FALSE)
chordDiagram(mat)
circos.clear()

circos.par(start.degree = 85)
chordDiagram(mat, order = c(rev(colnames(mat)), rev(rownames(mat))))


####change the color ???
#grid.col for sector colors, col for links between sectors
grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
             E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
chordDiagram(mat, grid.col = grid.col)
chordDiagram(t(mat), grid.col = grid.col)

### changes links
chordDiagram(mat, grid.col = grid.col, transparency = 0.5)
# change according to my choice
col_fun = colorRamp2(range(mat), c("#FFEEEE", "#FF0000"), transparency = 0.5)
chordDiagram(mat, grid.col = grid.col, col = col_fun)

# to highlight the links
border_mat2 = matrix("black", nrow = 1, ncol = ncol(mat))
rownames(border_mat2) = rownames(mat)[2]
colnames(border_mat2) = colnames(mat)
chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.border = border_mat2)

# To specify color of links ????
col_mat = rand_color(length(mat), transparency = 0.5)
dim(col_mat) = dim(mat)  # to make sure it is a matrix
col_mat[mat < 12] = "#00000000"
chordDiagram(mat, grid.col = grid.col, col = col_mat) 
### ?????????? can be used
col_fun = function(x) ifelse(x < 12, "#00000000", "#FF000080")
chordDiagram(mat, grid.col = grid.col, col = col_fun)

# wide links forward and narrow links backward ???
chordDiagram(mat, grid.col = grid.col, transparency = 0)
chordDiagram(mat, grid.col = grid.col, transparency = 0, link.zindex = rank(mat))


####imp note but not used
# Row names and column names in mat can also overlap. In this case, showing direction of 
#the link is important to distinguish them

#???? direction form columns to rows
chordDiagram(mat, grid.col = grid.col, directional = 1) ## -1 or 1

# 
df = expand.grid(letters[1:3], LETTERS[1:4])
df1 = df
df1$value = sample(10, nrow(df), replace = TRUE)
df2 = df
df2$value = -sample(10, nrow(df), replace = TRUE)
df = rbind(df1, df2)
grid.col = structure(1:7, names = c(letters[1:3], LETTERS[1:4]))
chordDiagram(df, col = ifelse(df$value > 0, "red", "green"), grid.col = grid.col)

circos.clear()
library(RColorBrewer)
####3 try on my own data #amino
library(readxl)
amino_acids <- read_excel("amino_acids_circlize.xlsx")
amino_acids = as.data.frame(amino_acids)
rownames(amino_acids) = amino_acids[,1]
amino_acids[,1] = NULL
amino_acids = as.matrix(amino_acids)

circos.par(start.degree = 40,track.margin = c(0.01, 0.05),
           track.height = 0.05,
           gap.after = c(rep(3, nrow(amino_acids)-1), 20, "L-Alanine"	= 2, "Hydroxyproline"= 2,	"Creatine"=2,
                         "Pyroglutamic acid" = 2,"Glycine"=2,"Orthohydroxyphenyl acetic acid"=20,
                         "Gut.K00259"=2,"Gut.K00262"=2,	"Gut.K00265"=2,"Gut.K00266"=2,
                         "Gut.K00278"=2,"Gut.K00609"=2,	"Gut.K00610"=2,"Gut.K00764"=2,
                         "Gut.K00811"=2,"Gut.K00820"=2,"Gut.K01424"=2,"Gut.K01425"=2,
                         "Gut.K01580"=2,"Gut.K01744"=2,"Gut.K01755"=2,"Gut.K01756"=2,
                         "Gut.K01914"=2,"Gut.K01915"=2,	"Gut.K01939"=2,"Gut.K01940"=2,
                         "Gut.K01953"=2,"Gut.K01955"=2,"Gut.K01956"=2,"Gut.K09758"=2,
                         "Gut.K11358"=2,"Gut.K14260"=2,"Gut.K00145"=2,"Gut.K00147"=2,
                         "Gut.K00286"=2,"Gut.K00619"=2,"Gut.K00657"=2,"Gut.K00818"=2,
                         "Gut.K00819"=2,	"Gut.K00926"=2,"Gut.K00930"=2,"Gut.K00931"=2,
                         "Gut.K01438"=2,"Gut.K01585"=2,	"Gut.K01777"=2,	"Gut.K10536"=2,
                         "Gut.K12251"=2,	"Gut.K13747"=2,	"Gut.K00031"=2,"Gut.K00033"=2,"Gut.K00036"=2,	"Gut.K01270"=2,
                         "Gut.K01917"=2,	"Gut.K00018"=2,	"Gut.K00058"=2,"Gut.K00133"=2,"Gut.K00281"=2,	"Gut.K00302"=2,
                         "Gut.K00600"=2,	"Gut.K00831"=2,	"Gut.K00865"=2,"Gut.K00928"=2,"Gut.K00998"=2,	"Gut.K01079"=2,
                         "Gut.K01620"=2,	"Gut.K01695"=2,	"Gut.K01696"=2,"Gut.K01733"=2,"Gut.K01752"=2,	"Gut.K01754"=2,
                         "Gut.K01834"=2,	"Gut.K02203"=2,	"Gut.K06001"=2,"Gut.K06718"=2,"Gut.K12524"=2,	"Gut.K00588"=2,
                         "Gut.K00817"=2,	"Gut.K00956"=2,	"Gut.K00957"=2,"Gut.K01008"=2,
                         "Gut.K14155" = 10,"Nares.K00259"=2, "Nares.K00608"=2,"Nares.K01424"=2,
                         "Nares.K01744"=2,	"Nares.K13821"=2,	"Nares.K00318"=2,
                         "Nares.K01485"=2,"Nares.K01611"=2,	"Nares.K08687"=2,"Nares.K00031"=2,
                         "Nares.K00033"=2,"Nares.K01256"=2,"Nares.K01920"=2,"Nares.K00018"=2,
                         "Nares.K00281"=2,"Nares.K00306"=2,	"Nares.K00605"=2,"Nares.K00865"=2,
                         "Nares.K00872"=2,"Nares.K01752"=2,	"Nares.K02204"=2,"Nares.K02437"=2,
                         "Nares.K01739"=2,"Nares.K11717"=20))

           


grid.col = c("Arginine and proline metabolism" = "#FC8D62",
             "Glutathione metabolism" = "#FC8D62",
             "Glycine, serine and threonine metabolism" = "#FC8D62",
             "Phenylalanine metabolism" = "#FC8D62",
             "Selenocompound metabolism" = "#FC8D62", 
             "Alanine, aspartate and glutamate metabolism" = "#FC8D62")


col_fun = function(x) ifelse(x > 0, "red", "blue")



#2
chordDiagram(amino_acids, col = col_fun, grid.col = grid.col, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.15, 0.15),preAllocateTracks=2.5, title(main = "Amino acid metabolism pathways"))

metabolome = amino_acids[,1:6]
human = amino_acids[,1:6]
Microbiome = amino_acids[,7:105]
gut = amino_acids[,7:81]
nares = amino_acids[,82:105]


#highlight.sector(colnames(gut), track.index = 3, col = "blue", 
                 #text = "Gut", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)
#highlight.sector(colnames(nares), track.index = 3, col = "blue", 
                 #text = "Nares", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)
highlight.sector(colnames(Microbiome), track.index = 3, col = "blue", 
                 text = "Microbiome", cex = 0.8, text.col = "white", niceFacing = TRUE, font=2)
highlight.sector(colnames(metabolome), track.index = 3, col = "springgreen4", 
                 text = "Metabolome", cex = 0.55, text.col = "white", niceFacing = TRUE, font = 2)
#highlight.sector(colnames(human), track.index = 2, col = "springgreen4", 
                 #text = "Human", cex = 0.6, text.col = "white", niceFacing = TRUE, font = 2)



circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.8)
}, bg.border = NA)




circos.clear()
library(ComplexHeatmap)

lgd_lines = Legend(at = c("upregulated", "downregulated"), type = "lines", 
                   legend_gp = gpar(col = c("red","blue"), lwd = 2))
lgd_list_vertical = packLegend(lgd_lines)
lgd_list_vertical
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left","center"))


#carbohydrates
library(readxl)
carbohydrates <- read_excel("carbohydrates_circlize.xlsx")
carbohydrates = as.data.frame(carbohydrates)
rownames(carbohydrates) = carbohydrates[,1]
carbohydrates[,1] = NULL
carbohydrates = as.matrix(carbohydrates)

circos.par(start.degree = 70,track.margin = c(0.01, 0.05),
           track.height = 0.07,
           gap.after = c(rep(3, nrow(carbohydrates)-1), 20, "PFKL"	= 2, "ENO1"= 2,	"GAPDHS"=2,
                         "LDHAL6B" = 2,"ACAT2"=10,"Glycine"=2,"Pregnanediol-3-glucuronide"=20,
                         "Gut.K00001"=2,"Gut.K00845"=2,"Gut.K00850"=2,"Gut.K00873"=2,
                         "Gut.K00918"=2,"Gut.K00927"=2,"Gut.K01610"=2,"Gut.K01689"=2,
                         "Gut.K01785"=2,"Gut.K01803"=2,"Gut.K01810"=2,"Gut.K01834"=2,
                         "Gut.K01895"=2,"Gut.K01905"=2,"Gut.K04041"=2,"Gut.K06859"=2,
                         "Gut.K00018"=2,"Gut.K00024"=2,"Gut.K00042"=2,"Gut.K00600"=2,
                         "Gut.K00865"=2,"Gut.K01091"=2,"Gut.K01625"=2,"Gut.K01647"=2,
                         "Gut.K01681"=2,"Gut.K01847"=2,"Gut.K01915"=2,"Gut.K00012"=2,
                         "Gut.K00040"=2,"Gut.K00041"=2,"Gut.K00854"=2,"Gut.K01051"=2,
                         "Gut.K01195"=2,"Gut.K01685"=2,"Gut.K01686"=2,"Gut.K01786"=2,
                         "Gut.K01804"=2,"Gut.K01805"=2,"Gut.K01812"=2,"Gut.K00625"=2,
                         "Gut.K00656"=2,"Gut.K00925"=2,"Gut.K01006"=2,"Gut.K01026"=2,
                         "Gut.K01572"=2,"Gut.K01573"=2,"Gut.K01649"=2,"Gut.K01676"=2,
                         "Gut.K01960"=2,"Gut.K01961"=2,"Gut.K15024"=10,"Nares.K00016"=2,
                         "Nares.K00850"=2,"Nares.K00886"=2,"Nares.K01623"=2,"Nares.K01785"=2,
                         "Nares.K01835"=2,"Nares.K03841"=2,"Nares.K00018"=2,"Nares.K00865"=2,
                         "Nares.K01681"=2,"Nares.K02437"=2,"Nares.K05606"=2,"Nares.K00853"=2,
                         "Nares.K01786"=2,"Nares.K01804"=2,"Nares.K00029"=2,
                         "Nares.K01679"=2,"Nares.K11263"=20))


					 	


grid.col = c("Glycolysis / Gluconeogenesis" = "#e78a9e",
             "Glyoxylate and dicarboxylate metabolism" = "#e78a9e",
             "Pentose and glucuronate interconversions" = "#e78a9e",
             "Pyruvate metabolism" = "#e78a9e")




col_fun = function(x) ifelse(x > 0, "red", "blue")

#1                        
chordDiagram(carbohydrates, col = col_fun, grid.col = grid.col, 
             directional = 1, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.09, 0.04))

#2
chordDiagram(carbohydrates, col = col_fun, grid.col = grid.col, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.15, 0.15),preAllocateTracks=2.5, title(main = "Carbohydrate metabolism pathways"))

RNA = carbohydrates[,1:5]
metabolome = carbohydrates[,6:7]
#human = carbohydrates[,1:6]
Microbiome = carbohydrates[,8:76]
# gut = carbohydrates[,7:81]
# nares = carbohydrates[,82:105]


#highlight.sector(colnames(gut), track.index = 3, col = "blue", 
#text = "Gut", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)
#highlight.sector(colnames(nares), track.index = 3, col = "blue", 
#text = "Nares", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)

highlight.sector(colnames(Microbiome), track.index = 3, col = "blue", 
                 text = "Microbiome", cex = 0.8, text.col = "white", niceFacing = FALSE, font=2)
highlight.sector(colnames(metabolome), track.index = 3, col = "springgreen4", 
                 text = "Metabolome", cex = 0.5, text.col = "white", facing = "reverse.clockwise", font = 2)
highlight.sector(colnames(RNA), track.index = 3, col = "springgreen4", 
                 text = "RNA-Seq", cex = 0.8, text.col = "white", niceFacing = TRUE, font = 2)
#highlight.sector(colnames(human), track.index = 2, col = "springgreen4", 
#text = "Human", cex = 0.6, text.col = "white", niceFacing = TRUE, font = 2)



circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.8)
}, bg.border = NA)

circos.clear()

#lipid
library(readxl)
lipid <- read_excel("lipid_circlize.xlsx")
lipid = as.data.frame(lipid)
rownames(lipid) = lipid[,1]
lipid[,1] = NULL
lipid = as.matrix(lipid)

circos.par(start.degree = 70,track.margin = c(0.01, 0.05),
           track.height = 0.07,
           gap.after = c(rep(3, nrow(lipid)-1), 20, "LPCAT1"	= 2, "DGKE"= 2,	"LPCAT3"=10,
                         "Glycine"=2,"Glycerylphosphorylethanolamine"=2,"Sulfogalactosylceramide"=20,
                         "Gut.K00655"=2,"Gut.K00901"=2,"Gut.K00968"=2,
                         "Gut.K00981"=2,"Gut.K00998"=2,"Gut.K01613"=2,
                         "Gut.K06131"=2,"Gut.K01442"=2,"Gut.K01190"=2,
                         "Gut.K01201"=10,"Nares.K01058"=2,"Nares.K01095"=2,
                         "Nares.K06131"=2,"Nares.K04708"=20))




grid.col = c("Ether lipid metabolism" = "#8DA0CB",
             "Glycerophospholipid metabolism" = "#8DA0CB",
             "Primary bile acid biosynthesis" = "#8DA0CB",
             "Sphingolipid metabolism" = "#8DA0CB")



col_fun = function(x) ifelse(x > 0, "red", "blue")

#1                        
chordDiagram(lipid, col = col_fun, grid.col = grid.col, 
             directional = 1, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.09, 0.04))

#2
chordDiagram(lipid, col = col_fun, grid.col = grid.col, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.1, 0.1),preAllocateTracks=2.5, title(main = "Lipid metabolism pathways"))

RNA = lipid[,1:3]
metabolome = lipid[,4:6]
#human = lipid[,1:6]
Microbiome = lipid[,7:20]
# gut = lipid[,7:81]
# nares = lipid[,82:105]


#highlight.sector(colnames(gut), track.index = 3, col = "blue", 
#text = "Gut", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)
#highlight.sector(colnames(nares), track.index = 3, col = "blue", 
#text = "Nares", cex = 0.6, text.col = "white", niceFacing = TRUE, font=2)

highlight.sector(colnames(Microbiome), track.index = 3, col = "blue", 
                 text = "Microbiome", cex = 0.8, text.col = "white", niceFacing = FALSE, font=2)
highlight.sector(colnames(metabolome), track.index = 3, col = "springgreen4", 
                 text = "Metabolome", cex = 0.7, text.col = "white", niceFacing = TRUE, font = 2)
highlight.sector(colnames(RNA), track.index = 3, col = "springgreen4", 
                 text = "RNA-Seq", cex = 0.7, text.col = "white", niceFacing = TRUE, font = 2)
#highlight.sector(colnames(human), track.index = 2, col = "springgreen4", 
#text = "Human", cex = 0.6, text.col = "white", niceFacing = TRUE, font = 2)



circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.8)
}, bg.border = NA)

circos.clear()

#others
library(readxl)
others <- read_excel("others_pathways_circlize.xlsx")
others = as.data.frame(others)
rownames(others) = others[,1]
others[,1] = NULL
others = as.matrix(others)

circos.par(start.degree = 50,track.margin = c(0.01, 0.05),
           track.height = 0.07,
           gap.after = c(rep(3, nrow(others)-1), 20,"PTPRB"=2,"ACTB"=2,"WASF2"=2,
                         "IQGAP1"=2, "ACTG1"=2, "PTPRJ"=2,"TGFBR2"=2,"RHOA"=2,
                        "HLA-DPA1"=2, "HLA-DQA1"=2,"HLA-B"=2,"HLA-E"=2,"HLA-A"=2,
                       "HLA-DQB1"=2,"HLA-C"=2, "CLTB"=2,"HCLS1"=2,"ARPC5"=2,"PXN"=2,
                       "SOAT1"=2,"LRP1"=2,"CCL5"=2,"TNFRSF1B"=2,"PPBP"=2,"CX3CR1"=2,
                       "IL4R"=2,"IL10RA"=2,"IL17RC"=2,"CSF2RB"=2,"IL17RA"=2,"BMP6"=2,
                       "CD27"=2,"GDF9"=2,	"TLN1"=2,	"FLNA"=2,	"GRB2"=2,	"VAV1"=2,
                       "ZYX"=2,	"THBS3"=2,	"MYLK"=2,	"RAF1"=2,	"CHAD"=2,	"COL9A3"=2,
                       "S100A8"=2,	"S100A9"=2,	"ANAPC5"=2,	"UBC"=2,	"UBB"=2,	"ARF6"=2,
                       "CD14"=2,	"MTOR"=2,	"CYTH1"=2,	"NCK2"=2,	"ZAP70"=2,	"LCK"=10,
                       "VCL"=2,	"LPA"=10,	"Pantothenate"=2,	"Glycine"=2,
                       "Biliverdin"=10,	"CD40L"=2,	"IL2"=2,	"IL5"=2,	"IFNB1"=2,	"MIP1B"=20,
                       "Gut.K00077"=2,	"Gut.K00826"=2,	"Gut.K00859"=2,	"Gut.K00954"=2,
                       "Gut.K01579"=2,	"Gut.K01652"=2,	"Gut.K01653"=2,	"Gut.K01687"=2,
                       "Gut.K01918"=2,	"Gut.K03525"=2,	"Gut.K09680"=2,	"Gut.K13038"=2,
                       "Gut.K00595"=2,	"Gut.K00768"=2,	"Gut.K01195"=2,	"Gut.K01599"=2,
                       "Gut.K01719"=2,	"Gut.K02188"=2,	"Gut.K02190"=2,	"Gut.K02224"=2,
                       "Gut.K02227	"=2,"Gut.K02230"=2,	"Gut.K02231"=2,	"Gut.K02232"=2,
                       "Gut.K02233"=2,	"Gut.K03394"=2,	"Gut.K05934"=10,	"Nares.K00867"=2,
                       "Nares.K00231"=2,	"Nares.K00595"=2,	"Nares.K00768"=2,	"Nares.K13541"=2,
                       "Nares.K13542"=20))




grid.col = c("Adherens junction" = "#66C2A5",
             "Allograft rejection" = "#66C2A5",
             "Autoimmune thyroid disease" = "#66C2A5",
             "Bacterial invasion of epithelial cells" = "#66C2A5",
            "Cholesterol metabolism" = "#66C2A5",
            "Cytokine-cytokine receptor interaction" = "#66C2A5",
            "Focal adhesion" = "#66C2A5",
            "IL-17 signaling pathway" = "#66C2A5",
            "Pantothenate and CoA biosynthesis" = "#66C2A5",
            "Porphyrin and chlorophyll metabolism" = "#66C2A5",
            "Shigellosis" = "#66C2A5",
            "T cell receptor signaling pathway" = "#66C2A5")



col_fun = function(x) ifelse(x > 0, "red", "blue")

#1                        
chordDiagram(others, col = col_fun, grid.col = grid.col, 
             directional = 1, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.09, 0.04))

#2
chordDiagram(others, col = col_fun, grid.col = grid.col, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.1, 0.1),preAllocateTracks=2.5, title(main = "Other pathways"))

RNA = others[,1:55]
protein = others[,56:57]
metabolome = others[,58:60]
cytokine = others[,61:65]
#human = others[,1:6]
Microbiome = others[,66:98]
# gut = others[,7:81]
# nares = others[,82:105]




highlight.sector(colnames(Microbiome), track.index = 3, col = "blue", 
                 text = "Microbiome", cex = 0.95, text.col = "white", niceFacing = TRUE, font=2)
highlight.sector(colnames(metabolome), track.index = 3, col = "springgreen4", 
                 text = "Metabolome", cex = 0.55, text.col = "white", facing = "reverse.clockwise", font = 2)
highlight.sector(colnames(RNA), track.index = 3, col = "springgreen4", 
                 text = "RNA-Seq", cex = 0.95, text.col = "white", niceFacing = TRUE, font = 2)
highlight.sector(colnames(protein), track.index = 3, col = "springgreen4", 
                 text = "Proteome", cex = 0.65, text.col = "white", facing = "reverse.clockwise", font = 2)
highlight.sector(colnames(cytokine), track.index = 3, col = "springgreen4", 
                 text = "Cytokine", cex = 0.755, text.col = "white", niceFacing = TRUE, font = 2)




circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.8)
}, bg.border = NA)

library(ComplexHeatmap)

lgd_lines = Legend(at = c("upregulated", "downregulated"), type = "lines", 
                   legend_gp = gpar(col = c("red","blue"), lwd = 2))
lgd_list_vertical = packLegend(lgd_lines)
lgd_list_vertical
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left","bottom"))

circos.clear()
