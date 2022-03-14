#### Circlize figures to integrate between pathways and omics features

library(circlize)
library(tidyverse)
library(RColorBrewer)
#### try on my own data 
library(readxl)
amino_acids <- read_excel("amino_acids_circlize.xlsx")
amino_acids = as.data.frame(amino_acids)
rownames(amino_acids) = amino_acids[,1]
amino_acids[,1] = NULL
amino_acids = as.matrix(amino_acids)

# adjust distances between pathways and features
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

           

# colors of pathways
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

# Determine the sectors of omics data
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

# adjust the distances in the final figure
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.8)
}, bg.border = NA)
circos.clear()


# These steps repeated on other pathways categories as the circlize package can't be used in loop form
