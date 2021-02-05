######################################################
########### ANALYSIS OF CRC DATA #####################
######################################################

library("stats")
library("methods")
library("utils")
library("scde") 
library("RColorBrewer") 
library("gplots") 
library("Biostrings")
library("Seurat")
library("dplyr")
library("Matrix")
library("data.table")
library("metaMA")
library("RMThreshold")

folderpath <- "D:/Teresa/Colon-final/"


######################################################
########### 1. Data Overview #########################
######################################################

# Load data for all 12 patients

printf <- function(...) cat(sprintf(...))

## CRITS
total_reads <- 100000
coding_genes <- 1000
mito_frac <- 15
non_genic_frac <- 100

ncores <- 12
oldpids <- c("90121", "90130", "92047", "92492", "92074", "92589", "91526", "100399", "91492", "92507", "92460", "92489")
pids <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12")
colors <- brewer.pal(12, "Set3")

cds <- read.csv(paste(folderpath, "resources/GRCh37_protein_coding.csv", sep=""))[,1]

print("Loading expression matrices...")
col <- c()
for (i in 1:12) {
  printf("PID: %s\n", oldpids[i])
  load(sprintf(paste(folderpath, "rawdata/expression-%s-%d_%d_%d_%d.Rdata", sep=""), oldpids[i], total_reads, coding_genes, mito_frac, non_genic_frac))
  
  outliers <- c()
  if (oldpids[i] == "92074") {
    outliers <- c("R21C10", "R3C37", "R27C61", "R32C4")
  } else if (oldpids[i] == "92492") {
    outliers <- c("R31C23", "R17C19", "R70C33", "R19C1")
  } else if (oldpids[i] == "92489") {
    outliers <- c("R24C60")
  } else if (oldpids[i] == "90121") {
    outliers <- c("R49C39")
  }
  
  outliers.idx <- unlist(lapply(strsplit(colnames(expression.matrices$count), "-"), tail, 1)) %in% outliers
  expression.matrices$count <- expression.matrices$count[,!outliers.idx]
  coding_idx <- gsub("\\..*$","",rownames(expression.matrices$count)) %in% cds
  expression.matrices$count <- expression.matrices$count[coding_idx, ]
  colnames(expression.matrices$count) <- paste0(pids[i], "-", sub(".*\\-","",colnames(expression.matrices$count)))

  if (i == 1) {
    count.matrix <- expression.matrices$count
    rownames(count.matrix) <- make.unique(as.character(gene.names)[coding_idx])
  } else {
    count.matrix <- cbind(count.matrix, expression.matrices$count)
  }
  col <- c(col, rep(colors[i], ncol(expression.matrices$count)))
}

write.csv(count.matrix, file = paste(folderpath, "rawdata/AllPatients_counts_noMeanCen_noLogTrans_incl12.csv", sep=""))


################################################################################################################
########### 2 Seurat processing of all 12 patients (with and without mean centering) #########################
################################################################################################################

# # if necessary: load count.matrix from csv file
# count.matrix <- fread(file = paste(folderpath, "rawdata/AllPatients_counts_noMeanCen_noLogTrans_incl12.csv", sep=""), header = T)    
# c <- as.data.table(count.matrix)
# all.gene.names <- t(c[,1])
# c[,1] <- NULL
# count.matrix <- as.matrix(c)
# rownames(count.matrix) <- all.gene.names
# rm(c)


### ---------------- cpm normalisation -------------------------------
tot.exp <- colSums(count.matrix)/1e6
cpm.matrix <- scale(count.matrix, center = F, scale = tot.exp)
# Calculate average expression of each gene  
AggregateExpression <- log(rowMeans(cpm.matrix)+1, base=2)
all.gene.names <- rownames(count.matrix)
names(AggregateExpression) <- all.gene.names
# Exclude genes with average expression below 3.5 --> 8222 genes
count.matrix <- count.matrix[names(AggregateExpression)[!(AggregateExpression<3.5)], ]
cpm.matrix <- cpm.matrix[names(AggregateExpression)[!(AggregateExpression<3.5)], ]
log.cpm.matrix <- log(cpm.matrix/10+1, base=2)
# Mean center 
patnames <- sub("\\-.*", "", x = colnames(log.cpm.matrix))
patnames.un <- unique(patnames)
log.cpm.matrix.mc <- NULL
for (i in 1:length(patnames.un)) {
  to.mc <- log.cpm.matrix[,patnames == patnames.un[i]]
  to.mc <- t(scale(t(to.mc), center = T, scale = F))
  log.cpm.matrix.mc <- cbind(log.cpm.matrix.mc, to.mc)
}



### ----------------- Create Seurat objects for PCA, plotting etc.  (colon.nmc = not mean centered) ----------
colon.nmc <- CreateSeuratObject(log.cpm.matrix, min.cells = 0, min.genes = 0, is.expr = -Inf, project = "Colon", names.field = 1, names.delim = "-")
colon.nmc@var.genes <- rownames(log.cpm.matrix)
colon.nmc@data <- colon.nmc@raw.data
colon.nmc@scale.data <- colon.nmc@raw.data


colon.mc <- CreateSeuratObject(log.cpm.matrix.mc, min.cells = 0, min.genes = 0, is.expr = -Inf, project = "Colon", names.field = 1, names.delim = "-")
colon.mc@var.genes <- rownames(log.cpm.matrix.mc)
colon.mc@data <- colon.mc@raw.data
colon.mc@scale.data <- colon.mc@raw.data



### ----------------- percent.mito per cell -----------------------

mito.genes <- grep(pattern = "^MT-", x = rownames(x = colon.mc@data), value = TRUE)
percent.mito <- Matrix::colSums(colon.mc@raw.data[mito.genes, ])/Matrix::colSums(colon.mc@raw.data)

# AddMetaData adds columns to object@meta.data
colon.mc <- AddMetaData(object = colon.mc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = colon.mc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


### ----------------------- PCA and tSNE ------------------------------

# Select mean-centered or not data
colon <- colon.nmc     # select colon.mc or colon.nmc



# Perform PCA and clustering
colon <- RunPCA(object = colon, pc.genes = colon@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
colon <- ProjectPCA(object = colon, do.print = FALSE)

# Look at standard deviations of PCs to draw cutoff where there is a clear elbow in the graph
filename <- paste0(folderpath,"PCElbowplot.pdf")
pdf(filename,width=4,height=4,paper='special') 
PCElbowPlot(object = colon)  
dev.off()

sigPC = 10

# Viz PCA top genes
filename <- paste0(folderpath, "VizPCATopGenes.pdf")
pdf(filename,width=4,height=8,paper='special') 
VizPCA(object = colon, pcs.use = 1:9)
dev.off()


# Perform clustering
colon <- FindClusters(object = colon, reduction.type = "pca", dims.use = 1:sigPC, 
                      resolution = 0.8, print.output = 0, save.SNN = TRUE)
# check number of clusters
length(unique(colon@ident))
# save cluster IDs 
colon <- StashIdent(object = colon, save.name = "ClusterNames_0.8")


# Perform tSNE
colon <- RunTSNE(object = colon, dims.use = 1:sigPC, do.fast = TRUE)

# Save Seurat object
save(colon, file = paste0(folderpath, "Colon_SeuratObject_after_tSNE.Rdata"))

# Save gene loadings of all genes in PCs 1-20
for (pc in 1:20){
  gene.loadings <- colon@dr$pca@gene.loadings.full[,pc]
  gene.loadings.rank <- gene.loadings
  gene.loadings <- gene.loadings[order(-gene.loadings.rank)]
  write.csv(gene.loadings, file=paste0(folderpath,sprintf("GeneLoadings_full_PC_%s.csv",pc)))
}



### ----------------------- PCA and tSNE plots ------------------------------

# colon <- load(paste0(folderpath, "12patients_MC/Colon_SeuratObject_after_tSNE.Rdata"))

# Heatmap of cells and genes ordered according to PCA scores 
# (cells.use = n 'extreme' cells on both ends of the spectrum)
filename <- paste0(folderpath, "PCHeatmap.pdf")

pdf(filename,width=8,height=8,paper='special') 
PCHeatmap(object = colon, 
          pc.use = 1:9, 
          cells.use = NULL, 
          do.balanced = T, 
          label.columns = F,
          disp.min = -5,
          disp.max = 5)
dev.off()


# PCA plot labelled by Cluster IDs and Patient IDs
plot1 <- PCAPlot(object = colon, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = F, do.label = F)
plot2 <- PCAPlot(object = colon, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "orig.ident", no.legend = F, do.label = F)
filename <- paste0(folderpath, "PCAPlot_byClusterID_byPatientID.pdf")
pdf(filename,width=16,height=7,paper='special') 
plot_grid(plot1, plot2, labels = c("Cluster ID", "Patient ID"), hjust = -1, vjust = 2)
dev.off()


# tSNE plot labelled by Cluster IDs and Patient IDs
plot1 <- TSNEPlot(object = colon, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = F, do.label = F)
plot2 <- TSNEPlot(object = colon, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "orig.ident", no.legend = F, do.label = F)
filename <- paste0(folderpath, "TSNEPlot_byClusterID_byPatientID.pdf")
pdf(filename,width=16,height=7,paper='special') 
plot_grid(plot1, plot2, labels = c("Cluster ID", "Patient ID"), hjust = -1, vjust = 2)
dev.off()


colon.nmc <- colon # colon.mc or colon.mmc 



### ----------------------- DEA for each patient ------------------------------

colon@meta.data$PatID <- as.character(sub("\\-.*","",colnames(colon@data)))

colon <- SetAllIdent(object = colon, id = "PatID")

  all_markers <- FindAllMarkers(colon, genes.use = NULL, logfc.threshold = 0.25,
                                test.use = "wilcox", min.pct = 0.1, min.diff.pct = -Inf,
                                print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf,
                                return.thresh = 0.05, do.print = TRUE, random.seed = 1,
                                min.cells = 3, latent.vars = "nUMI", assay.type = "RNA")
  
  clusterIDs <- unique(all_markers$cluster)
  
  for (c in 1:length(clusterIDs)){
    # sort DE genes by avg_logFC
    DEgenes_up <- all_markers[which(as.logical((all_markers$cluster==clusterIDs[c])*(all_markers$avg_logFC>0))),]
    DEgenes_up <- DEgenes_up[order(-DEgenes_up$avg_logFC),]
    DEgenes_down <- all_markers[which(as.logical((all_markers$cluster==clusterIDs[c])*(all_markers$avg_logFC<0))),]
    DEgenes_down <- DEgenes_down[order(DEgenes_down$avg_logFC),]
    write.csv(DEgenes_up, file=paste0(folderpath, sprintf("Cluster_%s_",clusterIDs[c]), "DEgenes_up.csv"))
    write.csv(DEgenes_down, file=paste0(folderpath, sprintf("Cluster_%s_",clusterIDs[c]), "DEgenes_down.csv"))

  }
  

    save(all_markers, file=paste0(folderpath, sprintf("AllMarkers_eachPatient.Rdata")))

    
### ----------------------- DEA for each cluster ------------------------------

  colon <- SetAllIdent(object = colon, id = "ClusterNames_0.6")
  
  
  
  all_markers <- FindAllMarkers(colon, genes.use = NULL, logfc.threshold = 0.25,
                                test.use = "wilcox", min.pct = 0.1, min.diff.pct = -Inf,
                                print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf,
                                return.thresh = 0.05, do.print = TRUE, random.seed = 1,
                                min.cells = 3, latent.vars = "nUMI", assay.type = "RNA")
  
  clusterIDs <- unique(all_markers$cluster)
  
  for (c in 0:(length(clusterIDs)-1)){
    # sort DE genes by avg_logFC
    DEgenes_up <- all_markers[which(as.logical((all_markers$cluster==c)*(all_markers$avg_logFC>0))),]
    DEgenes_up <- DEgenes_up[order(-DEgenes_up$avg_logFC),]
    DEgenes_down <- all_markers[which(as.logical((all_markers$cluster==c)*(all_markers$avg_logFC<0))),]
    DEgenes_down <- DEgenes_down[order(DEgenes_down$avg_logFC),]
    write.csv(DEgenes_up, file=paste0(folderpath, sprintf("Cluster_%s_",c), "DEgenes_up.csv"))
    write.csv(DEgenes_down, file=paste0(folderpath, sprintf("Cluster_%s_",c), "DEgenes_down.csv"))
  }
  
  save(all_markers, file=paste0(folderpath, sprintf("AllMarkers_eachCluster.Rdata")))



  ### ----------------------- hierarchical clustering based on top 10/ 20 / 50 DE genes per patient ------------------------------
  
  
Ngenes <- 10   # number of top genes to include

# collate data for heatmap
  
  load("D:/Teresa/Colon-final/DEA_eachPatient/AllMarkers_eachPatient.Rdata")    # loads all_markers
  clusterIDs <- unique(all_markers$cluster)
for (c in 1:length(clusterIDs)){
  # sort DE genes by avg_logFC
  DEgenes_up <- all_markers[which(as.logical((all_markers$cluster==clusterIDs[c])*(all_markers$avg_logFC>0))),]
  DEgenes_up <- DEgenes_up[order(-DEgenes_up$avg_logFC),]
  if (c == 1){
    top.genes <- DEgenes_up$gene[1:Ngenes]
  } else {
    top.genes <- c(top.genes, DEgenes_up$gene[1:Ngenes])
  }
}
top.genes <- unique(top.genes)

colon <- colon.nmc
hm.data <- colon@data[top.genes,]


# collate annotations for heatmap

MetaData.PatID_long <- pids
MetaData.SampleOrigin <- c("Liver metastasis", "Lung metastasis", "Liver metastasis", "Liver metastasis", 
                                  "Liver metastasis", "Primary tumor", "Liver metastasis", "Liver metastasis",
                                  "Primary tumor", "Primary tumor", "Primary tumor", "Primary tumor")
MetaData.MS <- c("MSS", "MSS", "MSS", "MSI", "MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "MSI")
# MetaData.CMS <- c("3","4","2","3","2","3","2","4","4","1","3","3") # NAs as best guess
MetaData.CMS <- c("3","NA","2","3","NA","3","NA","4","4","NA","3","3") # NAs as NAs

MetaData <- data.frame(PatID_long = MetaData.PatID_long, 
                       SampleOrigin = MetaData.SampleOrigin, 
                       MS = MetaData.MS, 
                       CMS = MetaData.CMS)

for (i in 1:dim(colon@meta.data)[1]){
  id <- which(sub("\\-.*","",MetaData.PatID_long) == colon@meta.data$orig.ident[i])
  if (i==1){
    MetaData_allCells <- MetaData[id,]
  } else{
    MetaData_allCells <- rbind(MetaData_allCells, MetaData[id,])
  }
}

rownames(MetaData_allCells) <- colon@cell.names
colon <- AddMetaData(colon, MetaData_allCells)
hm.ann <- colon@meta.data[, c("PatID_long", "SampleOrigin", "MS", "CMS")]

nmf.options(grid.patch=TRUE)

# Specify colors for meta data in heatmap
VarColors1 = brewer.pal(12, "Set3")
names(VarColors1) = unique(hm.ann[,1])
VarColors2 = brewer.pal(3, "Dark2")
names(VarColors2) = unique(hm.ann[,2])
VarColors3 = brewer.pal(2, "Dark2")
names(VarColors3) = unique(hm.ann[,3])
VarColors4 = brewer.pal(4, "Dark2")
VarColors4[[2]] <- "grey"
names(VarColors4) = unique(hm.ann[,4])
ann_colors = list(PatID_long = VarColors1,
                  SampleOrigin = VarColors2,
                  MS = VarColors3,
                  CMS = VarColors4)

# plot heatmap

# select which cells to include
#cellIDs <- c(1:10,501:510,1001:1010,1501:1510,2001:2010)
cellIDs <- 1:4661

pdf("D:/Teresa/Colon-final/heatmap.pdf",width=20,height=9,paper='special') 
aheatmap(hm.data[,cellIDs],
         color = "-RdYlBu2:100", breaks = NA, border_color = NA,
         scale = "none",
         distfun = "correlation", hclustfun = "average",
         cellwidth = 0.2, cellheight = 2,
         Rowv = TRUE, Colv = TRUE,
         annCol = hm.ann[cellIDs,],
         annColors = ann_colors,
         cexRow = 1.2, cexCol = 1.2, labCol = NA
)
dev.off()


### ------------ hierarchical clustering based on top 30 PC1+/- genes for mean-centered data -----------

Ngenes <- 30   # number of top genes to include

# collate data for heatmap

top.genes <- as.character(read.csv("D:/Teresa/Colon-final/12Patients_MC/GeneLoadings_full_PC_1.csv")$X)  # loads all_markers
top.genes <- c(top.genes[1:30], top.genes[(length(top.genes)-29):length(top.genes)])

colon <- colon.mc
hm.data <- colon@data[top.genes,]


# collate annotations for heatmap

MetaData.PatID_long <- pids
MetaData.SampleOrigin <- c("Liver metastasis", "Lung metastasis", "Liver metastasis", "Liver metastasis", 
                           "Liver metastasis", "Primary tumor", "Liver metastasis", "Liver metastasis",
                           "Primary tumor", "Primary tumor", "Primary tumor", "Primary tumor")
MetaData.MS <- c("MSS", "MSS", "MSS", "MSI", "MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "MSI")
# MetaData.CMS <- c("3","4","2","3","2","3","2","4","4","1","3","3") # NAs as best guess
MetaData.CMS <- c("3","NA","2","3","NA","3","NA","4","4","NA","3","3") # NAs as NAs

MetaData <- data.frame(PatID_long = MetaData.PatID_long, 
                       SampleOrigin = MetaData.SampleOrigin, 
                       MS = MetaData.MS, 
                       CMS = MetaData.CMS)

for (i in 1:dim(colon@meta.data)[1]){
  id <- which(sub("\\-.*","",MetaData.PatID_long) == colon@meta.data$orig.ident[i])
  if (i==1){
    MetaData_allCells <- MetaData[id,]
  } else{
    MetaData_allCells <- rbind(MetaData_allCells, MetaData[id,])
  }
}

rownames(MetaData_allCells) <- colon@cell.names
colon <- AddMetaData(colon, MetaData_allCells)
hm.ann.col <- data.frame(colon@meta.data[, "PatID_long"])
colnames(hm.ann.col) <- "PatientID"

hm.ann.rows <- data.frame(c(rep("up", 30), rep("down", 30)))
colnames(hm.ann.rows) <- "PC1"


# Specify colors for meta data in heatmap
VarColors1 = brewer.pal(12, "Set3")
names(VarColors1) = unique(as.matrix(hm.ann.col))
VarColors2 = brewer.pal(2, "Set3")
names(VarColors2) = c("up", "down")
ann_colors = list(PatientID = VarColors1,
                  PC1 = VarColors2)



nmf.options(grid.patch=TRUE)

# plot heatmap

pdf("D:/Teresa/Colon-final/heatmap_PC1_centerRow.pdf",width=20,height=9,paper='special') 
aheatmap(hm.data,
         color = "-RdYlBu2:100", breaks = seq(from = -2, to = 2, by = 0.04), border_color = NA,
         scale = "row",
         distfun = "correlation", hclustfun = "average",
         cellwidth = 0.2, cellheight = 2,
         Rowv = TRUE, Colv = TRUE,
         annCol = hm.ann.col,
         annRow = hm.ann.rows,
         annColors = ann_colors,
         cexRow = 1.2, cexCol = 1.2, labCol = NA
)
dev.off()

pdf("D:/Teresa/Colon-final/heatmap_PC1.pdf",width=20,height=9,paper='special') 
aheatmap(hm.data,
         color = "-RdYlBu2:100", breaks = NA, border_color = NA,
         scale = "none",
         distfun = "correlation", hclustfun = "average",
         cellwidth = 0.2, cellheight = 2,
         Rowv = TRUE, Colv = TRUE,
         annCol = hm.ann.col,
         annRow = hm.ann.rows,
         annColors = ann_colors,
         cexRow = 1.2, cexCol = 1.2, labCol = NA
)
dev.off()





###############################################################################################################

## Check whether or not to include HD1509

idx <- (sub("\\-.*","",colnames(colon@scale.data)) == "HD1509")

a <- colon@scale.data["LGR5",idx]
hist(a)

a <- count.matrix["LGR5",idx]
hist(a)

sum(a>0)  # is 9
a[a>0]




######################################################################################

