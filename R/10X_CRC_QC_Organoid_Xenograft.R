# Loading Libraries for Seurat Run
library(Seurat)
library(dplyr)
library(cowplot)

### Organoids
o1.data <- Read10X(data.dir = "~/Data_Colon/O1/hg19/")
o1 <- CreateSeuratObject(counts = o1.data, min.cells = 3, min.features = 200, project = "o1")
o1$sample <- "o1"

o2.data <- Read10X(data.dir = "~/Data_Colon/O2/hg19/")
o2 <- CreateSeuratObject(counts = o2.data, min.cells = 3, min.features = 200, project = "o2")
o2$sample <- "o2"

o3.data <- Read10X(data.dir = "~/Data_Colon/O3/hg19/")
o3 <- CreateSeuratObject(counts = o3.data, min.cells = 3, min.features = 200, project = "o3")
o3$sample <- "o3"

# Visualizing nUMI, nGenes, %Mito, and filtering
o1[["percent.mt"]] <- PercentageFeatureSet(object = o1, pattern = "^MT-")
VlnPlot(object = o1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = o1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = o1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

o2[["percent.mt"]] <- PercentageFeatureSet(object = o2, pattern = "^MT-")
VlnPlot(object = o2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = o2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = o2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

o3[["percent.mt"]] <- PercentageFeatureSet(object = o3, pattern = "^MT-")
VlnPlot(object = o3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = o3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = o3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Pre-processing, normalization and finding variable genes
o1 <- subset(x = o1, subset = nFeature_RNA > 2000 & percent.mt < 10)
o1 <- NormalizeData(object = o1, verbose = FALSE)
o1 <- FindVariableFeatures(object = o1, selection.method = "vst", nfeatures = 2000)

o2 <- subset(x = o2, subset = nFeature_RNA > 2000 & percent.mt < 10)
o2 <- NormalizeData(object = o2, verbose = FALSE)
o2 <- FindVariableFeatures(object = o2, selection.method = "vst", nfeatures = 2000)

o3 <- subset(x = o3, subset = nFeature_RNA > 2000 & percent.mt < 10)
o3 <- NormalizeData(object = o3, verbose = FALSE)
o3 <- FindVariableFeatures(object = o3, selection.method = "vst", nfeatures = 2000)

### Xenografts
x1.data <- Read10X(data.dir = "~/Data_Colon/x1/hg19/")
x1 <- CreateSeuratObject(counts = x1.data, min.cells = 3, min.features = 200, project = "x1")
x1$sample <- "x1"

x2.data <- Read10X(data.dir = "~/Data_Colon/x2/hg19/")
x2 <- CreateSeuratObject(counts = x2.data, min.cells = 3, min.features = 200, project = "x2")
x2$sample <- "x2"

# Visualizing nUMI, nGenes, %Mito, and filtering
x1[["percent.mt"]] <- PercentageFeatureSet(object = x1, pattern = "^hg19-MT-")
VlnPlot(object = x1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = x1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = x1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

x2[["percent.mt"]] <- PercentageFeatureSet(object = x2, pattern = "^hg19-MT-")
VlnPlot(object = x2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
plot1 <- FeatureScatter(object = x2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = x2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Pre-processing, normalization and finding variable genes
x1 <- subset(x = x1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
x1 <- NormalizeData(object = x1, verbose = FALSE)
x1 <- FindVariableFeatures(object = x1, selection.method = "vst", nfeatures = 2000)

x2 <- subset(x = x2, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)
x2 <- NormalizeData(object = x2, verbose = FALSE)
x2 <- FindVariableFeatures(object = x2, selection.method = "vst", nfeatures = 2000)
