#######################################################################################################
##### Binary scoring of gene sets defining cell types / cell states using comparison to similarly expressed genes ###################
#######################################################################################################

library(monocle)
library(reshape2)
library(Seurat)
library(data.table)
library(RColorBrewer)
library(metaMA)
library(dplyr)
library(Matrix)
library(data.table)
library(RMThreshold)
require(RColorBrewer)


# -------------- prepare expression data -------------- 

# Read in raw counts
count.matrix <- fread(file = "D:/Teresa/Colon-final/rawdata/AllPatients_counts_noMeanCen_noLogTrans_incl12.csv", header = T)    
c <- as.data.table(count.matrix)
all.gene.names <- t(c[,1])
c[,1] <- NULL
count.matrix <- as.matrix(c)
rownames(count.matrix) <- all.gene.names
rm(c)

# Make into cpm values
tot.exp <- colSums(count.matrix)/1e6
cpm.matrix <- scale(count.matrix, center = F, scale = tot.exp)
# Calculate average expression of each gene  
AggregateExpression <- log(rowMeans(cpm.matrix)+1, base=2)
names(AggregateExpression) <- all.gene.names
# Exclude genes with average expression below 3.5 --> 8222 genes
count.matrix <- count.matrix[names(AggregateExpression)[!(AggregateExpression<3.5)], ]
cpm.matrix <- cpm.matrix[names(AggregateExpression)[!(AggregateExpression<3.5)], ]
log.cpm.matrix <- log(cpm.matrix/10+1, base=2)
patnames <- sub("\\-.*", "", x = colnames(log.cpm.matrix))
patnames.un <- unique(patnames)
log.cpm.matrix.mc <- NULL
for (i in 1:length(patnames.un)) {
  to.mc <- log.cpm.matrix[,patnames == patnames.un[i]]
  to.mc <- t(scale(t(to.mc), center = T, scale = F))
  log.cpm.matrix.mc <- cbind(log.cpm.matrix.mc, to.mc)
}

# select patients
pids <- c("HD1495-P1", "HD1509-P2", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2596-P6", "HD2779-P7", "HD2791-P8", "HD3192-P9", "HD3254-P10", "HD3371-P11", "POP1-P12")
selected_pids <- c("HD1495-P1", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2779-P7", "HD2791-P8", "HD3254-P10", "HD3371-P11")
idx <- sub("\\-.*","",colnames(log.cpm.matrix.mc)) %in% sub("\\-.*","",selected_pids) 
log.cpm.matrix.mc <- log.cpm.matrix.mc[,idx]
log.cpm.matrix <- log.cpm.matrix[,idx]


# -------------- construct control gene sets -------------- 

# Divide all genes into 25 bins of equal size based on aggregate expression levels
Expr_aggregate <- rowMeans(log.cpm.matrix)
sortIndices <- sort(Expr_aggregate, index.return=T)
bin.size <- round(length(Expr_aggregate)/25)
bin.membership <- NULL
for (i in 1:length(sortIndices$ix)){
  bin.membership[i] <- floor(sortIndices$ix[i]/length(Expr_aggregate)*25+1)
}
bin.membership[bin.membership==26] <- 25
names(bin.membership) <- names(Expr_aggregate)

# Genes in each expression bin
GenesInBin = list()
for (i in 1:25){
  GenesInBin[[i]] <- names(Expr_aggregate)[bin.membership==i]
}
# Number of genes in each expression bin
NrGenesInBin <- table(bin.membership)

# Read in gene lists defining cell states / cell types
CuratedGenes <- Definition_CellStates_GeneList   
#CuratedGenes <- Definition_CellTypes_GeneList

# Calculate average expression in control gene sets for each cell state / cell type
CtrlGeneExpression <- matrix(NA, nrow = ncol(log.cpm.matrix.mc), ncol = length(CuratedGenes))   # row: cell, col: factor
# For each cell
for (c in 1:ncol(log.cpm.matrix)){
  # For each factor
  for (f in 1:length(CuratedGenes)){
    # For each gene in gene set, randomly select 100 genes from the same expression bin
    # and add their expression ...
    bin.compare.TotalExpr <- 0
    for (g in 1:length(CuratedGenes[[f]])){
      bin.id <- bin.membership[as.character(CuratedGenes[[f]][g])]
      if (!is.na(bin.id)){     # some genes don't exist in expression matrix
        bin.compare.genes.idx <- sample(1:NrGenesInBin[bin.id],100)
        bin.compare.genes <- GenesInBin[[bin.id]][bin.compare.genes.idx]
        bin.compare.TotalExpr <- bin.compare.TotalExpr + sum(log.cpm.matrix.mc[as.vector(bin.compare.genes),c])
      }
    }    
    # ... and divide by total number to obtain average
    CtrlGeneExpression[c,f] <- bin.compare.TotalExpr/(length(CuratedGenes[[f]])*100)
    
  }
}

save(CtrlGeneExpression, file="D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM_new/DefinedCellStates_CtrlGeneExpression.Rdata")
#save(CtrlGeneExpression, file="D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM_new/DefinedCellTypes_CtrlGeneExpression.Rdata")



# -------------- calculate gene set scores for each cell -------------- 

GeneSetScores <- matrix(NA, nrow = ncol(log.cpm.matrix.mc), ncol = length(CuratedGenes))   # row: cell, col: factor

# For each cell
for (c in 1:ncol(log.cpm.matrix.mc)){
  # For each factor
  for (f in 1:length(CuratedGenes)){
    CuratedGenes_inMatrix <- as.character(CuratedGenes[[f]])
    CuratedGenes_inMatrix <- CuratedGenes_inMatrix[CuratedGenes_inMatrix %in% rownames(log.cpm.matrix.mc)]
    GeneSetScores[c,f] <- mean(log.cpm.matrix.mc[CuratedGenes_inMatrix,c])-CtrlGeneExpression[c,f]
  }
}
rownames(GeneSetScores) <- colnames(log.cpm.matrix.mc)
#colnames(GeneSetScores) <- c(paste0("Factor_", seq(1,dim(GeneSetScores)[2],1)))
colnames(GeneSetScores) <- c("CellCycle",	"OXPHOS",	"Hypoxia",	"Glycolysis",	"FattyAcid")    # cell states
#colnames(GeneSetScores) <- c("Stem",	"DCS",	"TA",	"T-diff")                                 # cell types


save(GeneSetScores, file="D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM_new/DefinedCellStates_GeneSetScores.Rdata")
#save(GeneSetScores, file="D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM_new/DefinedCellTypes_GeneSetScores.Rdata")



# -------------- set cutoff and determine significant gene sets for each cell -------------- 

# set cutoff for gene set scores by sd
sigma <- 1    # 1 -> 10406 sig. sets in 3795 cells, 2 -> 2863 sig. sets in 1871 cells
cutoffGeneSetScore <- mean(GeneSetScores) + sigma*sd(GeneSetScores)

# determine significant gene sets for each cell
significantGeneSets <- (GeneSetScores > cutoffGeneSetScore)

# number of significant gene sets per cell
hist(rowSums(significantGeneSets))
sum(rowSums(significantGeneSets)>0)

# number of cells in which gene set is enriched, per patient
enrichedCells <- matrix(NA, nrow = length(selected_pids), ncol = dim(GeneSetScores)[2])
pid.rowlabels <- sub("\\-.*","",rownames(significantGeneSets))
for (p in 1:nrow(enrichedCells)){
  for (f in 1:ncol(enrichedCells)) {
    enrichedCells[p,f] <- sum(significantGeneSets[(pid.rowlabels==sub("\\-.*","",selected_pids[p])),f])
  }
}

mypalette<-brewer.pal(dim(GeneSetScores)[2],"Dark2")
par(mfrow=c(2,ceiling(length(selected_pids)/2)))
par(las=2)  # set axis label orientation to 'always perpendicular to the axis')
for (p in 1:8){
  barplot(enrichedCells[p,]/sum(pid.rowlabels==sub("\\-.*","",selected_pids[p])),
          main=selected_pids[p],
          col=mypalette,
          ylim=c(0,0.4),
          names.arg=colnames(GeneSetScores),
          cex.names=0.7
  )
}



###################################### DELETE BELOW HERE ####################################

# ----------- determine overlap between significantly enriched gene sets across cells ---------- 

hmcol<-colorRampPalette(brewer.pal(11,"RdBu"))(101)

NrPatients <- length(selected_pids)
NrFactors <- dim(GeneSetScores)[2]

OverlapArray <- array(NA, dim=c(NrPatients,NrFactors,NrFactors))
for (p in 1:NrPatients){
  cells.to.inspect <- which(pid.rowlabels==sub("\\-.*","",selected_pids[p]))
  for (i in 1:NrFactors){
    for (j in 1:NrFactors){
      # calculate proportion of cells that express both gene sets i and j  (denominator: all cells that express either i or j)
      OverlapArray[p,i,j] <- sum(significantGeneSets[cells.to.inspect,i]*significantGeneSets[cells.to.inspect,j])/
        (sum(significantGeneSets[cells.to.inspect,i])+sum(significantGeneSets[cells.to.inspect,j]))
    }
  }
  pdf(sprintf("D:/Teresa/Colon-final/FactorScoring/Heatmap_GeneSetOverlap_Pat_%s.pdf", selected_pids[p]),width=6,height=6,paper='special') 
  heatmap(OverlapArray[p,,], 
          scale="none",
          main=selected_pids[p],
          symm=T,
          col=rev(hmcol)
  )
  dev.off()
}


# ----------- binary heatmap with positive factors per cell ---------- 

hmcol<-colorRampPalette(brewer.pal(11,"RdBu"))(101)

NrPatients <- length(selected_pids)
NrFactors <- dim(GeneSetScores)[2]

OverlapArray <- array(NA, dim=c(NrPatients,NrFactors,NrFactors))
for (p in 1:NrPatients){
  cells.to.inspect <- which(pid.rowlabels==sub("\\-.*","",selected_pids[p]))
  for (i in 1:NrFactors){
    for (j in 1:NrFactors){
      # calculate proportion of cells that express both gene sets i and j  (denominator: all cells that express either i or j)
      OverlapArray[p,i,j] <- sum(significantGeneSets[cells.to.inspect,i]*significantGeneSets[cells.to.inspect,j])/
        (sum(significantGeneSets[cells.to.inspect,i])+sum(significantGeneSets[cells.to.inspect,j]))
    }
  }
  pdf(sprintf("D:/Teresa/Colon-final/FactorScoring/Heatmap_GeneSetOverlap_Pat_%s.pdf", selected_pids[p]),width=6,height=6,paper='special') 
  heatmap(significantGeneSets[p,,], 
          scale="none",
          main=selected_pids[p],
          symm=T,
          col=rev(hmcol)
  )
  dev.off()
}


# ------------ Create heatmap --------------------------------------------------------

# collate data for heatmap

hm.data <- significantGeneSets
hm.data_log <- log(hm.data+0.001)

# load annotations for heatmap
load("D:/Teresa/Colon-final/MetaData_allCells_forHeatmaps.Rdata")

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
pdf("D:/Teresa/Colon-final/FactorScoring/heatmap_factors_binary.pdf",width=20,height=9,paper='special') 
aheatmap(hm.data,
         color = "-Spectral:100", breaks = NA, border_color = NA,
         scale = "none",
         distfun = "correlation", hclustfun = "complete",
         cellwidth = 0.2, cellheight = 15,
         Rowv = F, Colv = F,
         annCol = hm.ann.col,
         annColors = ann_colors,         
         annRow = NA,
         cexRow = 0.8, cexCol = 1.2, labCol = NA
)
dev.off()





pdf("D:/Teresa/Colon-final/FactorScoring/heatmap_factors_binary.pdf",width=5,height=9,paper='special') 
aheatmap(matrix(as.numeric(hm.data), nrow = nrow(hm.data),  ncol = ncol(hm.data)),
         color = c("white", "black"),
         breaks  = NA
)
dev.off()


binaryMat <- t(matrix(as.numeric(hm.data), nrow = nrow(hm.data),  ncol = ncol(hm.data)))
pdf("D:/Teresa/Colon-final/FactorScoring/heatmap_factors_binary.pdf",width=20,height=9,paper='special') 
NMF::aheatmap(binaryMat, 
              scale="none", 
              revC=TRUE, 
              main="Enriched factors per cell", 
              distfun = "euclidean", hclustfun = "complete",
              annCol = hm.ann.col,
              annColors = ann_colors,
              color = c("white", "black")
)
dev.off()





############################################################################################
################### Binary scoring of factors per cell using AUCell ###################
############################################################################################


# ------------------------- Read in required files --------------------

folderout <- "D:/Teresa/Colon-final/AUCell/"

factors <- read.table("D:/Teresa/Colon/180226_NNMF/CuratedFactorList_extended.txt", sep="\t", header=F)
load("D:/Teresa/Colon-final/12patients_MC/Colon_SeuratObject_after_tSNE.Rdata")



# ------------------------- AUCell ------------------------------------

library(AUCell)
library(GSEABase)

exprMatrix <- colon@data[unique(as.character(as.matrix(factors))),]

# Build gene-expression rankings for each cell  
# (ranked from highest to lowest value; genes with same expression value are shuffled)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
save(cells_rankings, file=paste0(folderout,"cells_rankings.RData"))

# Create list of gene sets
geneSets <- list()
for (i in 1:dim(factors)[2]){
  geneSets[[i]] <- GeneSet(unique(as.character(factors[,i])), setName=sprintf("geneSet_Factor_%s", i))
}

# Calculate enrichment for the gene signatures (AUC)
# ~ fraction of genes, within the top X genes in the ranking, that are included in the signature
# where 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used:
# by default, 5% of the total number of genes; common values may range from 1 to 20%
for (i in 1:dim(factors)[2]){
  cells_AUC[[i]] <- AUCell_calcAUC(geneSets[[i]], cells_rankings, aucMaxRank = ceiling(0.2*nrow(cells_rankings))) # 20% used here
}
save(cells_AUC, file=paste0(folderout,"cells_AUC_Factor_%s.RData"))

#  Determine the cells with the given gene signatures or active gene sets
set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC[[1]], plotHist=TRUE, assign=TRUE) 





cells_AUC <- matrix(NA, dim(exprMatrix)[1], dim(exprMatrix)[2])
cells_AUC_thresholds <- matrix(NA, 4,5)
cells_AUC_thresholds_nCells <- matrix(NA, 4, 14)
cells_AUC_thresholds_selected <- matrix(NA, 1, 14)
rownames(cells_AUC_thresholds) <- c("outlierOfGlobal", "Global_k1", "L_k2", "R_k3")







for (i in c(1,4,7,5,9:14,16:19)){
  if (i==1) {
    geneSets <- GeneSet(unique(as.character(as.matrix(factors[,c(1,2,3)]))), setName="geneSet") 
    j = 1} else if (i==4){
      geneSets <- GeneSet(unique(as.character(as.matrix(factors[,c(4,6,15)]))), setName="geneSet")
      j = 2} else if (i==7){
        geneSets <- GeneSet(unique(as.character(as.matrix(factors[,c(7,8)]))), setName="geneSet")
        j = 3} else {
          geneSets <- GeneSet(as.character(as.matrix(factors[,i])), setName="geneSet")
          if (i==5){
            j=4
          } else if (i < 15){
            j=i-4 } else {j=i-5 }
        }
  
  A <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = ceiling(0.2*nrow(cells_rankings)))
  cells_AUC[j,] <- a@assays$data$AUC
  
  # detect treshold automatically
  cells_assignment <- AUCell_exploreThresholds(A, plotHist=F, assign=TRUE, thrP=0.01) 
  cells_assignment$geneSet$aucThr$thresholds
  
  
  # for (c in 1:4){
  #   cells_AUC_thresholds[rownames(cells_AUC_thresholds)[c],j] <- 
  #               cells_assignment$geneSet$aucThr$thresholds[rownames(cells_AUC_thresholds)[c],"threshold"]
  #   cells_AUC_thresholds_nCells[rownames(cells_AUC_thresholds)[c],j] <- 
  #     cells_assignment$geneSet$aucThr$thresholds[rownames(cells_AUC_thresholds)[c],"nCells"]
  # }
  cells_AUC_thresholds_selected[1,j] <- cells_assignment$geneSet$aucThr$thresholds["Global_k1","threshold"]
  
  
  
}




## ----- Binarise to_plot_m matrix ------

to_plot_m_binary <- matrix(NA, ncol = dim(to_plot_m)[2], nrow = dim(to_plot_m)[1])
for (i in 1:dim(to_plot_m)[1]){
  to_plot_m_binary[i,] <- as.numeric(to_plot_m[i,]>cells_AUC_thresholds_selected[i])
}


rownames(to_plot_m_binary) <- rownames(to_plot_m) 
colnames(to_plot_m_binary) <- colnames(to_plot_m)


pdf("D:/Teresa/Colon/180226_NNMF/Heatmap_19Factors_merged_binary.pdf",width=20,height=10,paper='special') 
heatmap.3(as.matrix(to_plot_m_binary), F,
          Rowv = F,  
          Colv = T,
          dendrogram = "col",
          scale="none",
          trace="none",
          density.info="none",
          main="Mean factor expression",
          ColIndividualColors=col2[colPatIDs],
          margin.for.labRow = 2,
          margin.for.labCol = 0.0
)
dev.off()




### --------- Fraction of cells positive for each factor, for each patient -------------------


FactorFractions <- matrix(NA, nrow = 10, ncol = 14)
for (i in 1:14){
  for (p in 1:10){
    FactorFractions[p,i] <- sum(to_plot_m_binary[i,colPatIDs==p])/sum(colPatIDs==p) 
  }
}

barplot(to_plot_m_binary[], main=rownames(to_plot_m_binary)[i])

to_plot_m_binary



# by factor
par(mfrow=c(2,7))
for (i in 1:14){
  barplot(FactorFractions[,i], col = brewer.pal(10,"Set3"), main=rownames(to_plot_m_binary)[i], ylim=c(0.0,1.0))
}

# by patient
par(mfrow=c(2,5))
for (p in 1:10){
  barplot(FactorFractions[p,], col = brewer.pal(10,"Set3"), main=pids[p], ylim=c(0.0,1.0))
}

i=1
barplot(FactorFractions[,i], col = brewer.pal(10,"Set3"), main=rownames(to_plot_m_binary)[i], ylim=c(0.0,1.0), legend.text = pids)
barplot(FactorFractions[i,], col = brewer.pal(10,"Set3"), main=pids[p], ylim=c(0.0,1.0), legend.text = rownames(to_plot_m_binary))







