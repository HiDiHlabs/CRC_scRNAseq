### ----------------------- DEA for each patient separately ------------------------------

data.loc = ### INSERT DATA FOLDER HERE ###

colon_for_DEA <- colon

colon_for_DEA@meta.data$PatID <- as.character(sub("\\-.*","",colnames(colon_for_DEA@data)))

colon_for_DEA <- SetAllIdent(object = colon_for_DEA, id = "PatID")

for (i in 1:12){
  
  folderpath <- paste0(data.loc, "DEA_PatientsSeparately/",pids[i],"/")
  dir.create(folderpath, showWarnings = FALSE)
  
  # subset data from this patient
  patData <- SubsetData(colon_for_DEA, ident.use = sub("\\-.*","",pids[i]),
                        do.center = FALSE, do.scale = FALSE,
                        max.cells.per.ident = Inf, random.seed = 1, do.clean = T)
  patData <- NormalizeData(object = patData, 
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
  patData <- ScaleData(patData)
    
  # find variable genes
  patData <- FindVariableGenes(object = patData, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = -Inf, x.high.cutoff = Inf, y.cutoff = 0.5)
  
  length(patData@var.genes)
  
  ### ----------------------- PCA and tSNE ------------------------------
  
  # PCA and clustering for this patient
  patData <- RunPCA(object = patData, pc.genes = patData@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  patData <- ProjectPCA(object = patData, do.print = FALSE)
  
  # Look at standard deviations of PCs to draw cutoff where there is a clear elbow in the graph
  filename <- paste0(folderpath,"PCElbowplot.pdf")
  pdf(filename,width=4,height=4,paper='special') 
  PCElbowPlot(object = patData)  
  dev.off()
  
  sigPC = 10
  
  # Viz PCA top genes
  filename <- paste0(folderpath, "VizPCATopGenes.pdf")
  pdf(filename,width=4,height=8,paper='special') 
  VizPCA(object = patData, pcs.use = 1:9)
  dev.off()
  
  # Perform clustering
  patData <- FindClusters(object = patData, reduction.type = "pca", dims.use = 1:sigPC, 
                        resolution = 1.0, print.output = 0, save.SNN = TRUE)
  patData <- StashIdent(object = patData, save.name = "ClusterNames_1.0")
  
  # Perform tSNE
  patData <- RunTSNE(object = patData, dims.use = 1:sigPC, do.fast = TRUE)
  
  # Save Seurat object
  save(patData, file = paste0(folderpath, "patData_SeuratObject_after_tSNE.Rdata"))
  
  # Save gene loadings of all genes in PCs 1-20
  for (pc in 1:20){
    gene.loadings <- patData@dr$pca@gene.loadings.full[,pc]
    gene.loadings.rank <- gene.loadings
    gene.loadings <- gene.loadings[order(-gene.loadings.rank)]
    write.csv(gene.loadings, file=paste0(folderpath,sprintf("GeneLoadings_full_PC_%s.csv",pc)))
  }
  
  ### ----------------------- PCA and tSNE plots ------------------------------
  
  # PCA plot labelled by Cluster IDs
  filename <- paste0(folderpath, "PCAPlot_byClusterID.pdf")
  pdf(filename,width=8,height=7) 
  PCAPlot(object = patData, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "ClusterNames_1.0", no.legend = F, do.label = F)
  dev.off()
  
  # tSNE plot labelled by Cluster IDs 
  filename <- paste0(folderpath, "TSNEPlot_byClusterID.pdf")
  pdf(filename,width=8,height=7) 
  TSNEPlot(object = patData, dim.1 = 1, dim.2 = 2, do.return = TRUE, group.by = "ClusterNames_1.0", no.legend = F, do.label = F)
  dev.off()
  
  ### -------------------------- DEA for this patient ------------------------------------
  
  all_markers <- FindAllMarkers(patData, genes.use = NULL, logfc.threshold = 0.25,
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
  
}
