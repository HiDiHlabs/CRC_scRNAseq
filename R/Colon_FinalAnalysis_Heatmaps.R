########################################################################################################
########################## 4. Heatmaps of curated factors across all cells ###################################
########################################################################################################

require(MetaDE)
require(data.table)
require(RColorBrewer)
require(colorspace)
require(gplots)
require(GMD)
require(qlcMatrix)

# ------------ Load mean-centered colon data   (which we used for NNMF) ----------------------------------------

data.loc = ### INSERT DATA FOLDER HERE ###

counts <- read.csv(paste0(data.loc, "NNMF_LGR5/NMFinput_LGR5.csv"), header = TRUE, row.names = 1)

pids <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12")
selected_pids <- c("P1", "P3", "P4", "P5", "P7", "P8", "P10", "P11")

factors <- read.table(paste0(data.loc, "NNMF_LGR5/CuratedFactorList.txt"), sep="\t", header=T)
NrFactors <- dim(factors)[2]
factorIDs <- as.numeric(sub(".*\\_","",colnames(factors)))

colNames <- colnames(counts)
colPatIDs <- sub("-.*$","",colnames(counts))
colPatIDs <- as.numeric(factor(colPatIDs))

# ------------ Create matrix of cell scores:  (A) using factor scores for each cell ---------------------------------

CellScoresMatrix <- matrix(NA, nrow=NrFactors, ncol=dim(counts)[2])

for (f in 1:NrFactors){
  scores.in <- as.matrix(read.csv(sprintf(paste0(data.loc, "NNMF_LGR5/Cells_k25_Factor_%s.csv"), factorIDs[f]), row.names=1))
  for (c in 1:length(colNames)){
    CellScoresMatrix[f,c] <- scores.in[gsub("\\.","\\-",colNames[c]),1]
  }
}

rownames(CellScoresMatrix) <- c("G2/M", "G1/S", "MYC",
                                "OXPHOS_1", "OXPHOS_2", "Fatty_Acid", "Hypoxia/Glycolysis_1",  "Hypoxia/Glycolysis_2",
                                "TNFa/NFkB_1", "TNFa/NFkB_2", "Stress", 
                                "Stem", "DCS")
colnames(CellScoresMatrix) <- colNames

save(CellScoresMatrix, file=paste0(data.loc, "CellScoresMatrix_factorScores.Rdata"))

# ------------ Create matrix of cell scores:  (B) using mean expression of top N genes in factor for each cell ---------------------------------

NtopGenes <- 200

CellScoresMatrix_tG <- matrix(NA, nrow=NrFactors, ncol=dim(counts)[2])

for (f in 1:NrFactors){
  for (c in 1:length(colNames)){
    CellScoresMatrix_tG[f,c] <- mean(counts[as.character(factors[1:NtopGenes,f]),c])
  }
}

rownames(CellScoresMatrix_tG) <- c("G2/M", "G1/S", "MYC",
                                "OXPHOS_1", "OXPHOS_2", "Fatty_Acid", "Hypoxia/Glycolysis_1",  "Hypoxia/Glycolysis_2",
                                "TNFa/NFkB_1", "TNFa/NFkB_2", "Stress", 
                                "Stem", "DCS")

colnames(CellScoresMatrix_tG) <- colNames

save(CellScoresMatrix_tG, file=paste0(data.loc, "CellScoresMatrix_top200genes.Rdata"))

# ------------ Create heatmap --------------------------------------------------------

# collate data for heatmap

hm.data <- CellScoresMatrix_tG
hm.data_log <- log(hm.data+0.001)

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

MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]
# drop factor levels that are not in the selected data frame
MetaData_selectedCells <- lapply(MetaData_selectedCells, function(x) if(is.factor(x)) factor(x) else x)

hm.ann.col <- data.frame(MetaData_selectedCells)

# Specify colors for meta data in heatmap
# Specify colors for meta data in heatmap
VarColors1 = brewer.pal(8, "Set3")
names(VarColors1) = unique(hm.ann.col[,1])
VarColors2 = brewer.pal(3, "Dark2")
names(VarColors2) = unique(hm.ann.col[,2])
VarColors3 = brewer.pal(2, "Dark2")
names(VarColors3) = unique(hm.ann.col[,3])
VarColors4 = brewer.pal(4, "Dark2")
VarColors4[[3]] <- "grey"
names(VarColors4) = unique(hm.ann.col[,4])
ann_colors = list(PatID_long = VarColors1,
                  SampleOrigin = VarColors2,
                  MS = VarColors3,
                  CMS = VarColors4)

nmf.options(grid.patch=TRUE)

# plot heatmap

library(colorRamps)
cmap <- matlab.like2(101)

x_cutoff <- 1.5

pdf(paste0(data.loc, "heatmap_curatedFactors_top200g_allCells_correlation_complete_differentColorScale.pdf"),width=20,height=9,paper='special') 
aheatmap(hm.data,
         color = cmap, breaks = seq(0, x_cutoff, x_cutoff/100), border_color = NA,
         scale = "none",
         distfun = "correlation", hclustfun = "complete",
         cellwidth = 0.2, cellheight = 15,
         Rowv = TRUE, Colv = TRUE,
         annCol = hm.ann.col,
         annRow = NA,
         annColors = ann_colors,
         cexRow = 0.8, cexCol = 1.2, labCol = NA
)
dev.off()
