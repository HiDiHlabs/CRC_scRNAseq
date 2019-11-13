# Loading Libraries for Seurat Run
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggpubr)

# Preparing the Seurat Objects
t1.data <- Read10X(data.dir = "~/Data_Colon/T1/hg19/")
t1 <- CreateSeuratObject(counts = t1.data, min.cells = 3, min.features = 200, project = "t1")
t1$sample <- "t1"
t1 <- RenameCells(t1, add.cell.id = "t1")

t2.data <- Read10X(data.dir = "~/Data_Colon/T2/hg19/")
t2 <- CreateSeuratObject(counts = t2.data, min.cells = 3, min.features = 200, project = "t2")
t2$sample <- "t2"
t2 <- RenameCells(t2, add.cell.id = "t2")

t3.data <- Read10X(data.dir = "~/Data_Colon/T3/hg19/")
t3 <- CreateSeuratObject(counts = t3.data, min.cells = 3, min.features = 200, project = "t3")
t3$sample <- "t3"
t3 <- RenameCells(t3, add.cell.id = "t3")

# Visualizing nUMI, nGenes, %Mito, and filtering
t1[["percent.mt"]] <- PercentageFeatureSet(object = t1, pattern = "^MT-")
VlnPlot(object = t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = t1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = t1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

t2[["percent.mt"]] <- PercentageFeatureSet(object = t2, pattern = "^MT-")
VlnPlot(object = t2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = t2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = t2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

t3[["percent.mt"]] <- PercentageFeatureSet(object = t3, pattern = "^MT-")
VlnPlot(object = t3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = t3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = t3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Pre-processing, normalization and finding variable genes
t1 <- subset(x = t1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
t1 <- NormalizeData(object = t1, verbose = FALSE)
t1 <- FindVariableFeatures(object = t1, selection.method = "vst", nfeatures = 2000)

t2 <- subset(x = t2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
t2 <- NormalizeData(object = t2, verbose = FALSE)
t2 <- FindVariableFeatures(object = t2, selection.method = "vst", nfeatures = 2000)

t3 <- subset(x = t3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
t3 <- NormalizeData(object = t3, verbose = FALSE)
t3 <- FindVariableFeatures(object = t3, selection.method = "vst", nfeatures = 2000)

# Perform integration
t_anchors <- FindIntegrationAnchors(object.list = list(t1, t2, t3), dims = 1:20) 
t_combined <- IntegrateData(anchorset = t_anchors, dims = 1:20)


# Selecting the assay slot from "RNA" to "integrated"
DefaultAssay(object = t_combined) <- "integrated"

# Run the standard workflow for visualization and clustering
all.genes <- rownames(x = t_combined)
t_combined <- ScaleData(object = t_combined, features = all.genes, verbose = FALSE)
t_combined <- RunPCA(object = t_combined, npcs = 30, verbose = FALSE)

DimHeatmap(object = t_combined, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(object = t_combined, dims = 16:30, cells = 500, balanced = TRUE)

# t-SNE and Clustering
t_combined <- RunUMAP(object = t_combined, reduction = "pca", dims = 1:14)
t_combined <- FindNeighbors(object = t_combined, reduction = "pca", dims = 1:14)
t_combined <- FindClusters(t_combined, resolution = 0.5)

p1 <- DimPlot(object = t_combined, reduction = "umap", group.by = "sample", pt.size = 0.5)
p2 <- DimPlot(object = t_combined, reduction = "umap", label = TRUE, pt.size = 0.5)
plot_grid(p1, p2)
DimPlot(object = t_combined, reduction = "umap", split.by = "sample", pt.size = 0.5)

# Save as an RData file
save(t_combined, file = "t_combined_cca.RData")

# Finding ALL markers
DefaultAssay(object = t_combined) <- "RNA"
t_combined_markers <- FindAllMarkers(object = t_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

### Identifying Cluster Specific Genes
DefaultAssay(object = t_combined) <- "RNA"

# Monocytes (but Macrophages are also CD68+)
FeaturePlot(object = t_combined,
            features = c("PTPRC", "CD4", "CD14", "LYZ"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

FeaturePlot(object = t_combined,
            features = c("CD68", "CD16", "CD14"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

# Activated Helper CD4+ T cells (CD3+, CD4 )
FeaturePlot(object = t_combined,
            features = c("PTPRC", "CD4", "CD3E", "CD40LG"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

# Cytotoxic CD8+ T cells (CD3+, CD8+, NKG7, GZMA)
FeaturePlot(object = t_combined,
            features = c("PTPRC", "CD3E", "CD8A", "NKG7"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

# B cells
# Plasma Cells (CD19-, SDC1+, CD27+, CD38+, TNFRSF17+, MS4A1-)
FeaturePlot(object = t_combined,
            features = c("PTPRC", "MS4A1", "CD19", "TNFRSF17"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

# Epithelial Cells
FeaturePlot(object = t_combined,
            features = c("EPCAM", "FERMT1", "VIL1", "CEACAM5"), 
            min.cutoff = "q9", 
            pt.size = 0.1)

# Fibroblasts
FeaturePlot(object = t_combined,
            features = c("VIM", "SPARC", "COL5A2", "MMP2"), 
            min.cutoff = "q9", 
            pt.size = 0.1, order = TRUE)

# All in one plot
FeaturePlot(object = t_combined,
            features = c("PTPRC", "CD4", "CD14", "LYZ", 
                         "PTPRC", "CD4", "CD3E", "CD40LG", 
                         "PTPRC", "CD3E", "CD8A", "NKG7", 
                         "PTPRC", "MS4A1", "CD19", "TNFRSF17",
                         "EPCAM", "FERMT1", "VIL1", "CEACAM5",
                         "VIM", "SPARC", "COL5A2", "MMP2"), 
            min.cutoff = "q9", 
            pt.size = 0.1,
            ncol = 4,
            order = TRUE)

###Subsetting the epithelial compartment
#EPCAM+/VIL1+/CEACAM5+/VIM-/SPARC-
subset_tumor_epithelial <- subset(t_combined, idents = '1')

# scale the data again
all.genes <- rownames(x = subset_tumor_epithelial)
subset_tumor_epithelial <- ScaleData(object = subset_tumor_epithelial, features = all.genes, verbose = FALSE)
save(subset_tumor_epithelial, file = "subset_tumor_epithelial.RData")

#PCA
subset_tumor_epithelial <- RunPCA(object = subset_tumor_epithelial, npcs = 30, verbose = FALSE)
DimHeatmap(object = subset_tumor_epithelial, dims = 1:6, balanced = TRUE)
DimHeatmap(object = subset_tumor_epithelial, dims = 7:12, balanced = TRUE)
DimHeatmap(object = subset_tumor_epithelial, dims = 13:18, balanced = TRUE)
DimHeatmap(object = subset_tumor_epithelial, dims = 19:24, balanced = TRUE)

#TSNE
DefaultAssay(object = subset_tumor_epithelial) <- "integrated"
subset_tumor_epithelial <- RunTSNE(object = subset_tumor_epithelial, reduction = "pca", dims = 1:12)
subset_tumor_epithelial <- FindNeighbors(object = subset_tumor_epithelial, reduction = "pca", dims = 1:12)
subset_tumor_epithelial <- FindClusters(subset_tumor_epithelial, resolution = 1.5)

#TSNE Plot
p1 <- DimPlot(object = subset_tumor_epithelial, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(object = subset_tumor_epithelial, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

save(subset_tumor_epithelial, file = "subset_tumor_epithelial.RData")

# Finding ALL markers
DefaultAssay(object = subset_tumor_epithelial) <- "RNA"
subset_tumor_epithelial_markers <- FindAllMarkers(object = subset_tumor_epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Finding Markers
DefaultAssay(object = subset_tumor_epithelial) <- "RNA"

FeaturePlot(object = subset_tumor_epithelial,
            features = c("SMOC2", "PROX1", "ASCL2", "RNF43"), 
            min.cutoff = "q9",
            pt.size = 1,
            reduction = "tsne")

FeaturePlot(object = subset_tumor_epithelial,
            features = c("DEFA5", "REG4", "FRZB"), 
            min.cutoff = "q9",
            pt.size = 1,
            reduction = "tsne")

FeaturePlot(object = subset_tumor_epithelial,
            features = c("MYC", "CDC20", "PA2G4"), 
            min.cutoff = "q9",
            pt.size = 1,
            reduction = "tsne")

FeaturePlot(object = subset_tumor_epithelial,
            features = c("CA9", "TFF3", "FABP1", "FXYD3"), 
            min.cutoff = "q9", 
            pt.size = 1,
            reduction = "tsne")

FeaturePlot(object = subset_tumor_epithelial,
            features = c("PROX1", "ASCL2", "SLC12A2",
                         "MYC", "CDC20", "MCM6",
                         "KRT19", "KRT20", "FABP1"), 
            min.cutoff = "q9", 
            pt.size = 1,
            reduction = "tsne",
            ncol = 3) 

# New Identities
new.cluster.ids <- c("HighMT", "Stem", "Stem", "TA", "Tdiff", "Tdiff", "Tdiff")
names(new.cluster.ids) <- levels(subset_tumor_epithelial)
subset_tumor_epithelial <- RenameIdents(subset_tumor_epithelial, new.cluster.ids)

###Subsetting the epithelial compartment WITHOUT HighMT fraction
#took cell IDs of cells that were not from HighMT
subset_tumor_epithelial_noHighMT <- subset(t_combined, cells = cell.id.epithelial.noHighMT)

# scale the data again
all.genes <- rownames(x = subset_tumor_epithelial_noHighMT)
subset_tumor_epithelial_noHighMT <- ScaleData(object = subset_tumor_epithelial_noHighMT, features = all.genes, verbose = FALSE)
save(subset_tumor_epithelial_noHighMT, file = "subset_tumor_epithelial_noHighMT.RData")

#PCA
subset_tumor_epithelial_noHighMT <- RunPCA(object = subset_tumor_epithelial_noHighMT, npcs = 30, verbose = FALSE)
DimHeatmap(object = subset_tumor_epithelial_noHighMT, dims = 1:6, balanced = TRUE)
DimHeatmap(object = subset_tumor_epithelial_noHighMT, dims = 7:12, balanced = TRUE)

#TSNE
DefaultAssay(object = subset_tumor_epithelial_noHighMT) <- "integrated"
subset_tumor_epithelial_noHighMT <- RunTSNE(object = subset_tumor_epithelial_noHighMT, reduction = "pca", dims = 1:10)
subset_tumor_epithelial_noHighMT <- FindNeighbors(object = subset_tumor_epithelial_noHighMT, reduction = "pca", dims = 1:10)
subset_tumor_epithelial_noHighMT <- FindClusters(subset_tumor_epithelial_noHighMT, resolution = 0.8)

#TSNE Plot
p1 <- DimPlot(object = subset_tumor_epithelial_noHighMT, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(object = subset_tumor_epithelial_noHighMT, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

save(subset_tumor_epithelial_noHighMT, file = "subset_tumor_epithelial_noHighMT.RData")

# Marker Genes
DefaultAssay(object = subset_tumor_epithelial_noHighMT) <- "RNA"
FeaturePlot(object = subset_tumor_epithelial_noHighMT,
            features = c("PROX1", "ASCL2", "SLC12A2",
                         "MYC", "CDC20", "MCM6",
                         "KRT19", "KRT20", "FABP1"), 
            min.cutoff = "q9", 
            pt.size = 1,
            reduction = "tsne",
            ncol = 3,
            order = TRUE) 

# New Identities
new.cluster.ids <- c("Stem", "Tdiff", "TA")
names(new.cluster.ids) <- levels(subset_tumor_epithelial_noHighMT)
subset_tumor_epithelial_noHighMT <- RenameIdents(subset_tumor_epithelial_noHighMT, new.cluster.ids)

save(subset_tumor_epithelial_noHighMT, file = "subset_tumor_epithelial_noHIghMT_new_ident.RData")
