# Quantification of a5c8 ratio from e12.5 single cell spinal motor neurons, specifically in LMC. 
# Figure 3E

#load necessary packages for analysis
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
#use_python("", required = FALSE)
# reticulate::use_python("/anaconda3/bin/python", required = FALSE)

# set working directory
# setwd("~/")

#load subset data containing LMC cells only
LMC <- readRDS(file = "LMC.rds")


# There are 1662 cells in LMC
dim(LMC@data)
# [1] 19034  1662

# Expression distribution of Hoxa5, Hoxc8, Hoxc4, Hoxc9 and Hoxa7 genes
which(row.names(LMC@data) == "Hoxa5") # [1] 6795
which(row.names(LMC@data) == "Hoxc8") # [1] 16438
which(row.names(LMC@data) == "Hoxc4") # [1] 16440
which(row.names(LMC@data) == "Hoxc9") # [1] 16437
which(row.names(LMC@data) == "Hoxa7") # [1] 6797
pdf('a5_normalized_distribution.pdf', width = 6, height = 5)
hist(as.numeric(LMC@data[6795, ]))
dev.off()
pdf('c8_normalized_distribution.pdf', width = 6, height = 5)
hist(as.numeric(LMC@data[16438, ]))
dev.off()
pdf('c4_normalized_distribution.pdf', width = 6, height = 5)
hist(as.numeric(LMC@data[16440, ]))
dev.off()
pdf('c9_normalized_distribution.pdf', width = 6, height = 5)
hist(as.numeric(LMC@data[16437, ]))
dev.off()
pdf('a7_normalized_distribution.pdf', width = 6, height = 5)
hist(as.numeric(LMC@data[6797, ]))
dev.off()


# create identity with threshold setting
# A cut-off to distinguish cells from each population could be set at 
# a local minimum between these two peak populations in the observed 
# single cell expression distribution. 
# To account for the relatively low expression levels of transcription factors, 
# we used the first local minimum closest to the “off” population. 
# An exception in this regard was made for Hoxc4, which was likely to be in excess of bimodal
# thus manual adjustment was applied
composite_ident <- data.frame(
  c4 = LMC@data["Hoxc4", ] > 0.6, 
  a7 = LMC@data["Hoxa7", ] > 0.6,
  c9 = LMC@data["Hoxc9", ] > 0.6 
)

LMC@meta.data$c4a7c9 <-  apply(composite_ident, 1, function(x) {
  possible_ids <- c("c4+", "a7+", "c9+")
  if (any(x)) {
    return(paste(possible_ids[x], collapse = " "))
  }
  return("other")
})
pdf('Manual_setting_c4a7c9.pdf', width = 8.27, height = 5.83)
DimPlot(LMC, group.by = "c4a7c9", reduction.use = "umap", pt.size = 1.3)
dev.off()


composite_ident <- data.frame(
  a5 = LMC@data["Hoxa5", ] > 0.6,
  c8 = LMC@data["Hoxc8", ] > 1.5
)

LMC@meta.data$a5c8 <-  apply(composite_ident, 1, function(x) {
  possible_ids <- c("a5+", "c8+")
  if (any(x)) {
    return(paste(possible_ids[x], collapse = " "))
  }
  return("other")
})
pdf('Manual_a5c8_forc4a7c9.pdf', width = 8.27, height = 5.83)
DimPlot(LMC, group.by = "a5c8", reduction.use = "umap", pt.size = 1.3)
dev.off()

#function for a5c8 ratio quantification in defined regions
a5c8ratio <- function(object = NULL, color = c("red", "orange", "green", "gray92")) {
  fn <- paste(deparse(substitute(object)))
  
  a5c8_num <- data.frame(table(object@meta.data$a5c8))
  a5c8_count <- mutate(a5c8_num, ratio = Freq/length(rownames(object@meta.data)))
  
  pdf(paste0("a5c8_ratio_for_", fn, ".pdf"), width = 8, height = 5)
  p <- ggplot(a5c8_count, aes(x = Var1, y = ratio, fill=Var1)) +
    geom_bar(stat = "identity") +
    ylim(0.0, 1.0) +
    scale_fill_manual(values= color) +
    ggtitle(paste0("a5c8 ratio for ", fn, " regions")) +
    labs(x = "", 
         y = "ratio (n/total)", 
         fill = "") + theme_cowplot() +
    geom_text(aes(label=round(ratio, digits = 2)), vjust=1.2, color="black", size=4.5)
  print(p)
  dev.off()
  
  write.csv(a5c8_count, file = paste0("a5c8_ratio_for", fn, ".csv"))
  return(object)
  
}

#subset c4a7c9 triple positive cells and perform quantification of a5,a5c8 and c8 cells ratio
# in this subset of cells
LMC <- SetAllIdent(LMC, id = "c4a7c9")
LMC_c4a7c9 <- SubsetData(LMC, ident.use = "c4+ a7+ c9+")
LMC_c4a7c9 <- SetAllIdent(LMC_c4a7c9, id = "a5c8")

pdf('GenePlot_c4a7c9_axes.pdf', width = 6, height = 6)
GenePlot(LMC_c4a7c9, gene1 = "Hoxa5", "Hoxc8", xlim = c(0.0, 3.0), ylim = c(0.0, 4.0)
         , col.use = c("red", "orange", "green", "gray92"))
dev.off()

# apply function a5c8ratio for quantification
a5c8ratio(LMC_c4a7c9)


#subset c4a7 double positive cells
LMC <- SetAllIdent(LMC, id = "c4a7c9")
LMC_c4a7 <- SubsetData(LMC, ident.use = "c4+ a7+")
LMC_c4a7 <- SetAllIdent(LMC_c4a7, id = "a5c8")

pdf('GenePlot_c4a7_axes.pdf', width = 6, height = 6)
GenePlot(LMC_c4a7, gene1 = "Hoxa5", "Hoxc8", xlim = c(0.0, 3.0), ylim = c(0.0, 4.0)
         , col.use = c("red", "orange", "green", "gray92"))
dev.off()


# quantify a5,a5c8 and c8 cells ratio
# for c4a7 double positive cells
a5c8ratio(LMC_c4a7)


#subset a7c9 double positive cells
LMC <- SetAllIdent(LMC, id = "c4a7c9")
LMC_a7c9 <- SubsetData(LMC, ident.use = "a7+ c9+")
LMC_a7c9 <- SetAllIdent(LMC_a7c9, id = "a5c8")

pdf('GenePlot_a7c9_axes.pdf', width = 6, height = 6)
GenePlot(LMC_a7c9, gene1 = "Hoxa5", "Hoxc8", xlim = c(0.0, 3.0), ylim = c(0.0, 4.0)
         , col.use = c("red", "orange", "green", "gray92"))
dev.off()


# quantify a5,a5c8 and c8 cells ratio
# for a7c9 double positive cells
a5c8ratio(LMC_a7c9)

#subset c4 positive cells
LMC <- SetAllIdent(LMC, id = "c4a7c9")
LMC_c4 <- SubsetData(LMC, ident.use = "c4+")
LMC_c4 <- SetAllIdent(LMC_c4, id = "a5c8")

pdf('GenePlot_c4_axes.pdf', width = 6, height = 6)
GenePlot(LMC_c4, gene1 = "Hoxa5", "Hoxc8", xlim = c(0.0, 3.0), ylim = c(0.0, 4.0)
         , col.use = c("red", "orange", "green", "gray92"))
dev.off()


# quantify a5,a5c8 and c8 cells ratio
# for c4 positive cells
a5c8ratio(LMC_c4)


#subset c9 positive cells
LMC <- SetAllIdent(LMC, id = "c4a7c9")
LMC_c9 <- SubsetData(LMC, ident.use = "c9+")
LMC_c9 <- SetAllIdent(LMC_c9, id = "a5c8")

pdf('GenePlot_c9_axes.pdf', width = 6, height = 6)
GenePlot(LMC_c9, gene1 = "Hoxa5", "Hoxc8", xlim = c(0.0, 3.0), ylim = c(0.0, 4.0)
         , col.use = c("orange", "green", "gray92"))
dev.off()


# quantify a5,a5c8 and c8 cells ratio
# for c9 positive cells
a5c8ratio(LMC_c9, color = c("orange", "green"))




