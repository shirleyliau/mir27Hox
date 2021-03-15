#load necessary packages for analysis
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
reticulate::use_python("/anaconda3/bin/python", required = TRUE)

# set working directory
#Read 10x dataset and create Seurat object
data_dir <- "./filtered_feature_bc_matrix/output"
data <- Read10X(data.dir = data_dir)
brac<- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, 
                          project = "brac.data")

# function
process <- function(object = NULL, res = 3.0,
                    sigPC = 40, cell = 500, pcs.compute = 40,
                    force.recalc = FALSE, filename = NULL) {
  if (is.null(filename)) {
    fn <- paste(deparse(substitute(object)), 
                format(Sys.time(), "%Y%m%d_%H%M%S"))
  } else {
    fn <- filename
  }
  
  obj <- object
  if(!dir.exists(paste0("./",fn))){
    dir.create(paste0("./",fn)) 
  }
  setwd(paste0("./",fn))
  obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR,
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  obj <- ScaleData(object =obj, vars.to.regress = c("nUMI", "percent.mito"),
                   genes.use = obj@var.genes, model.use = "negbinom")
  obj <- RunPCA(object = obj, do.print = TRUE, pcs.print = 1:5,
                genes.print = 5, pcs.compute = pcs.compute)
  obj <- ProjectPCA(object = obj, do.print = FALSE)
  obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = 1:sigPC, 
                      resolution = res, print.output = 0, save.SNN = TRUE,
                      force.recalc = force.recalc) 
  obj <- RunUMAP(object = obj, dims.use = 1:sigPC, min_dist = 0.5)
  
  pdf(paste0(fn,"_VariablegenePlot.pdf"), width = 8.27, height = 5.83)
  print(VariableGenePlot(object = obj))
  dev.off()
  
  pdf(paste0(fn,"_PCAPlot.pdf"), width = 8.27, height = 5.83)
  print(PCAPlot(object = obj, dim.1 = 1, dim.2 = 2, pt.size = 0.5))
  dev.off()
  
  pdf(paste0(fn,"_PCAHeatmap.pdf"), width = 35, height = 35)
  print(PCHeatmap(object = obj, pc.use = 1:pcs.compute, cells.use = cell, do.balanced = TRUE,
                  label.columns = FALSE, use.full = FALSE))
  dev.off()

  pdf(paste0(fn,"_UMAP.pdf"), width = 8.27, height = 5.83)
  DimPlot(obj, pt.size = 1.8, do.label=T, label.size = 5, reduction.use = "umap")
  dev.off()
  
  #saveRDS(obj, file = paste0(fn,".rds"))
  
  print("Scaling, feature selection and dimensionality reduction and clustering is done, 
        there should be Rdata, PCA plot and heatmap, UMAP in working folder.")
  setwd("../")
  return(obj)
}


#Data preprocessing
#quality filtering
mito.genes <- grep(pattern = "^mt-", x = rownames(x = brac@data), value = TRUE, ignore.case = T)
percent.mito <- Matrix::colSums(brac@raw.data[mito.genes, ])/Matrix::colSums(brac@raw.data)
brac <- AddMetaData(object =brac, metadata = percent.mito, col.name = "percent.mito")


# Figure S3a and S3b
pdf("Quality_filtering_Vlnplot.pdf", width = 8.27, height = 5.83)
VlnPlot(object =brac, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use = 0.1, cols.use = "dodgerblue3")
dev.off()
pdf("Quality_filtering_Scatterplot.pdf", width = 8.27, height = 5.83)
par(mfrow = c(1, 2))
GenePlot(object = brac, gene1 = "nUMI", gene2 = "percent.mito", col.use = "pink", cex.use = 0.5)
GenePlot(object = brac, gene1 = "nUMI", gene2 = "nGene", col.use = "lightblue", cex.use = 0.5)
dev.off()


# setting cutoffs for outliers (97% percentile)
UMI <- quantile(brac@meta.data$nUMI, c(.05, .10, .50, .90, .95, .97, .995, .997))
UMI
#      5%       10%       50%       90%       95%       97%     99.5%     99.7% 
# 2867.90  10606.60  29795.00  52524.40  64260.70  72198.48  99490.28 102781.63 
nGene <- quantile(brac@meta.data$nGene, c(.05, .10, .50, .90, .95, .97, .995, .997))
nGene
#         5%      10%      50%      90%      95%      97%    99.5%    99.7% 
#   1370.200 3621.800 5675.000 7152.400 7648.500 8002.000 8802.620 8957.758 
percent.mito <- quantile(brac@meta.data$percent.mito, c(.05, .10, .50, .90, .95, .97, .995, .997))
percent.mito
#           5%        10%        50%        90%        95%        97%      99.5%      99.7% 
#   0.01958109 0.02281037 0.03499278 0.05212886 0.07186809 0.28561851 0.75156154 0.79221367 
brac <- FilterCells(object = brac, subset.names = c("nGene", "nUMI", "percent.mito"),
                    low.thresholds = c(1300, -Inf, -Inf), high.thresholds = c(8000, 72200, 0.1))

saveRDS(brac, file = "brac_quality_filtering.rds")

#normalization of data
brachial <- NormalizeData(object =brac, normalization.method = "LogNormalize", 
                          scale.factor = 10000)

#feature selection (variable genes), scaling of data and dimensionality reduction
sigPC = 40
res = 0.18
brachial <- process(brachial, res = res, sigPC = sigPC, pcs.compute = 40)

# Run this only to determine significant PC
# brachial <- RunPCA(object = brachial, do.print = TRUE, pcs.print = 1:5,
#               genes.print = 5, pcs.compute = 100)
# brachial <- ProjectPCA(object = brachial, do.print = FALSE)
#Jackstraw
# brachial <- JackStraw(brachial, num.pc = 40, num.replicate = 100, prop.freq = 0.01)
# JackStrawPlot(brachial, plot.y.lim = 0.9, nCol=5, PCs = 1:40)

# #Elbow plot
# PCElbowPlot(brachial, num.pc = 100)


#visualize PC genes
pdf('brachial_VizPCA1-10.pdf', width = 14, height = 12)
VizPCA(object = brachial, pcs.use = 1:10, nCol = 2)
dev.off()
pdf('brachial_VizPCA11-20.pdf', width = 14, height = 12)
VizPCA(object = brachial, pcs.use = 11:20, nCol = 2)
dev.off()
pdf('brachial_VizPCA21-30.pdf', width = 14, height = 12)
VizPCA(object = brachial, pcs.use = 21:30, nCol = 2)
dev.off()
pdf('brachial_VizPCA31-40.pdf', width = 14, height = 12)
VizPCA(object = brachial, pcs.use = 31:40, nCol = 2)
dev.off()


#assign identity
currentids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
newids <- c("MMC", "Early neurons", "Interneurons", "LMCl", "LMCm", 
            "LMCl", "Slc17a6+", "HMC", "MMC", "MN_I", "LMCm")
brachial@ident <- plyr::mapvalues(x = brachial@ident, 
                                      from = currentids, to = newids)
myorder <- c("MMC", "LMCl",  "LMCm", "HMC", "Early neurons", "MN_I","Slc17a6+", "Interneurons")
brachial@ident <- factor(brachial@ident, levels = myorder)
brachial <- StashIdent(object = brachial, save.name = "ClusterNames_0.18")
saveRDS(brachial, file = "brachial.rds")



# Figure S3C and D - UMAP of all collected cells 
pdf("UMAP_rename_identity.pdf", width = 8.27, height = 5.83)
DimPlot(brachial, pt.size = 1.5, do.label=T, label.size = 5, reduction.use = "umap")
dev.off()

pdf("VlnPlot_MN_markers.pdf", width = 15, height = 12)
VlnPlot(brachial, c("Slc5a7", "Slc18a3", "Chat", "Mnx1", 
                    "Foxp1", "Aldh1a2", "Isl1", "Lhx1", 
                    "Lhx3", "Lhx4", "Ebf2", "Neurod2"), point.size.use = 0.05, nCol = 2)
dev.off()


# remove non Motor neurons
MN <- SubsetData(brachial, ident.remove = "Interneurons")
#rerun with similar pre-processing steps
MNsigPC = 40
res = 0.28
MN <- process(MN, res = res, sigPC = MNsigPC, force.recalc = TRUE)
currentids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
newids <- c("Early neurons", "MMC", "LMCl", "Slc17a6+", "LMCm", "HMC", 
            "LMCl", "MMC", "MN_I", "LMCl", "LMCm", "LMCm", "PGC")
MN@ident <- plyr::mapvalues(x = MN@ident, from = currentids, to = newids)
myorder <- c("MMC", "LMCl", "LMCm", "HMC", "Early neurons", "MN_I", "PGC", "Slc17a6+")
MN@ident <- factor(MN@ident, levels = myorder)
MN <- StashIdent(object = MN, save.name = "ClusterNames_MN_0.28")

#UMAP with assigned identities
pdf("UMAP_rename_MN.pdf", width = 8.27, height = 5.83)
DimPlot(MN, pt.size = 1.5, do.label=T, label.size = 5, reduction.use = "umap", 
        cols.use = c("gray92", "blue4", "dodgerblue", "gray92","gray92", "gray92", "gray92", "gray92"))
dev.off()
saveRDS(MN, file = "MN.rds")

#Motor column markers to check cluster identities assignment
pdf("VlnPlot_MN_columnmarkers.pdf", width = 15, height = 12)
VlnPlot(MN, c("Foxp1", "Aldh1a2", "Lhx1", "Isl1", "Lhx3", "Lhx4", "Ebf2", "Zeb2", "Nos1",
              "Slc5a7", "Slc17a6"), point.size.use = 0.05, nCol = 2)
dev.off()


# focusing on brachial motor neurons LMC for quantification of Hoxa5-on, Hoxc8-on, Hoxa5-on/c8-on ratio
#LMC subclustering
MN <-SetAllIdent(MN, id = "ClusterNames_MN_0.28")
LMC <- SubsetData(MN, ident.use = c("LMCm", "LMCl"))
LMCsigPC = 30
res = 0.78
LMC <- process(LMC, res = res, sigPC = LMCsigPC, force.recalc = TRUE)

# Run this only to check on PC used
# Determining significant PC (run before clustering)
# LMC <- RunPCA(object = LMC, pc.genes = LMC@var.genes, do.print = TRUE, pcs.compute = 100)
# PCAPlot(object = LMC, dim.1 = 1, dim.2 = 2)
# LMC <- ProjectPCA(object = LMC, do.print = T)

# #Jackstraw
# LMC <- JackStraw(LMC, num.pc = 40, num.replicate = 100, prop.freq = 0.01)
# JackStrawPlot(LMC, plot.y.lim = 0.9, nCol=5, PCs = 1:40)

# #Elbow
# PCElbowPlot(LMC, num.pc = 100)

myorder <- c("2", "4", "7", "0", "10", "8", "1", "5", "12", "11","9", "6", "3", "13")
LMC@ident <- factor(LMC@ident, levels = myorder)

#Figure 3C-UMAP of LMC only
pdf("UMAP_LMC.pdf", width = 8.27, height = 5.83)
DimPlot(LMC, pt.size = 1.5, do.label=T, label.size = 5, reduction.use = "umap")
dev.off()
#Figure 3D - expression of Hoxa5, Hoxc8 and motor pool marker
pdf('LMC_Hox_ETS.pdf', width = 5.83, height = 8.27)
VlnPlot(LMC,c("Hoxa5", "Hoxc8", "Etv4", "Pou3f1"), point.size.use = 0.05, nCol = 1)
dev.off()

saveRDS(LMC, file = "LMC.rds")







