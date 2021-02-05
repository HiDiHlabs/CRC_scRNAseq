#######################################################################################################
##### Binary scoring of factors per cell using comparison to similarly expressed genes ###################
##### using merged factors ###################
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

data.loc = ### INSERT DATA FOLDER HERE ###

# -------------- prepare expression data -------------- 

# Read in raw counts
count.matrix <- fread(file = paste0(data.loc, "rawdata/AllPatients_counts_noMeanCen_noLogTrans_incl12.csv"), header = T)    
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
pids <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12")
selected_pids <- c("P1", "P3", "P4", "P5", "P7", "P8", "P10", "P11")
idx <- sub("\\-.*","",colnames(log.cpm.matrix.mc)) %in% sub("\\-.*","",selected_pids) 
log.cpm.matrix.mc <- log.cpm.matrix.mc[,idx]

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

# Load merged factor gene lists - use merger of top genes for merged factors
load(paste0(data.loc, "MergedFactors_GeneLists.Rdata"))

# Calculate average expression in control gene sets for each factor and cell
CtrlGeneExpression <- matrix(NA, nrow = ncol(log.cpm.matrix.mc), ncol = length(MergedFactors_GeneLists))   # row: cell, col: factor
# For each cell
for (c in 1:ncol(log.cpm.matrix.mc)){
  # For each factor
  for (f in 1:length(MergedFactors_GeneLists)){
    # For each gene in gene set, randomly select 100 genes from the same expression bin
    # and add their expression ...
    bin.compare.TotalExpr <- 0
    for (g in 1:length(MergedFactors_GeneLists[[f]])){
      bin.id <- bin.membership[as.character(MergedFactors_GeneLists[[f]][g])]
      bin.compare.genes.idx <- sample(1:NrGenesInBin[bin.id],100)
      bin.compare.genes <- GenesInBin[[bin.id]][bin.compare.genes.idx]
      bin.compare.TotalExpr <- bin.compare.TotalExpr + sum(log.cpm.matrix.mc[as.vector(bin.compare.genes),c])
    }    
    # ... and divide by total number to obtain average
    CtrlGeneExpression[c,f] <- bin.compare.TotalExpr/(length(MergedFactors_GeneLists[[f]])*100)
    
  }
}

save(CtrlGeneExpression, file=paste0(data.loc, "FactorScoring/CtrlGeneExpression_mergedFactors.Rdata"))

# -------------- calculate gene set scores for each cell -------------- 

GeneSetScores <- matrix(NA, nrow = ncol(log.cpm.matrix.mc), ncol = length(MergedFactors))   # row: cell, col: factor

# For each cell
for (c in 1:ncol(log.cpm.matrix.mc)){
  # For each factor
  for (f in 1:length(MergedFactors)){
    GeneSetScores[c,f] <- mean(log.cpm.matrix.mc[as.character(MergedFactors_GeneLists[[f]]),c])-CtrlGeneExpression[c,f]
  }
}
rownames(GeneSetScores) <- colnames(log.cpm.matrix.mc)
colnames(GeneSetScores) <- c("Immune_response", "Hypoxia/Glycolysis", "Cell_cycle", "OXPHOS",
                                "MYC", "Stem", "DCS", "Fatty_acid")

save(GeneSetScores, file=paste0(data.loc, "FactorScoring/GeneSetScores_topGenesFromMergedFactors.Rdata"))

# -------------- set cutoff and determine significant gene sets for each cell -------------- 

sigma <- 1   
cutoffGeneSetScores <- rep(0,dim(GeneSetScores)[2])
for (c in 1:length(cutoffGeneSetScores)){
  cutoffGeneSetScores[c] <- mean(GeneSetScores[,c]) + sigma*sd(GeneSetScores[,c])
}

# determine significant gene sets for each cell
significantGeneSets <- matrix(NA,dim(GeneSetScores)[1],dim(GeneSetScores)[2])
for (c in 1:length(cutoffGeneSetScores)){
  significantGeneSets[,c] <- (GeneSetScores[,c] > cutoffGeneSetScores[c])
}
rownames(significantGeneSets) <- rownames(GeneSetScores)

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
par(mfrow=c(2,ceiling(length(selected_pids)/2)),
    mar=c(5.1,4.1,4.1,2.1),
    las=2   # set axis label orientation to 'always perpendicular to the axis')
    )
for (p in 1:8){
  barplot(enrichedCells[p,]/sum(pid.rowlabels==sub("\\-.*","",selected_pids[p])),
          main=selected_pids[p],
          col=mypalette,
          ylim=c(0,0.4),
          names.arg=colnames(GeneSetScores),
          cex.names=0.7
  )
}

# ----------- determine overlap between significantly enriched gene sets across cells ---------- 

hmcol<-colorRampPalette(brewer.pal(11,"RdBu"))(101)

FactorNames <- c("Immune_response", "Hypoxia/Glycolysis", "Cell_cycle", "OXPHOS",
                 "MYC", "Stem", "DCS", "Fatty_acid")

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
  pdf(sprintf(paste0(data.loc, "FactorScoring/Heatmap_GeneSetOverlap_Pat_%s_mergedFactors.pdf"), selected_pids[p]),width=6,height=6,paper='special') 
  heatmap(OverlapArray[p,,], 
          scale="none",
          main=selected_pids[p],
          symm=T,
          col=rev(hmcol),
          labRow=FactorNames,
          labCol=FactorNames,
          cexRow=1,
          cexCol=1,
          margins=c(10,10)
  )
  dev.off()
}

# ----------- binary heatmap with positive factors per cell ---------- 

# collate data for heatmap
hm.data <- significantGeneSets

# load annotations for heatmap
load(paste0(data.loc, "MetaData_allCells_forHeatmaps.Rdata"))

MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]
# drop factor levels that are not in the selected data frame
MetaData_selectedCells <- lapply(MetaData_selectedCells, function(x) if(is.factor(x)) factor(x) else x)

hm.ann.col <- data.frame(MetaData_selectedCells)

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

binaryMat <- t(matrix(as.numeric(hm.data), nrow = nrow(hm.data),  ncol = ncol(hm.data)))
pdf(paste0(data.loc, "FactorScoring/heatmap_mergedfactors_binary.pdf"),width=20,height=9,paper='special') 
NMF::aheatmap(binaryMat, 
              scale="none", 
              revC=TRUE, 
              main="Enriched factors per cell", 
              distfun = "euclidean", hclustfun = "complete",
              annCol = hm.ann.col,
              annColors = ann_colors,
              color = c("white", "black"),
              cexRow = 0.5,
              labRow = FactorNames,
              labCol = rep("",ncol(binaryMat))
)
dev.off()

# remove cells that do not score 'ON' for any factor and plot heatmap
binaryMat_rmOFF <- binaryMat[,!(colSums(binaryMat)==0)]
pdf(paste0(data.loc, "FactorScoring/heatmap_mergedfactors_binary_rmOFF.pdf"),width=20,height=9,paper='special') 
NMF::aheatmap(binaryMat_rmOFF, 
              scale="none", 
              revC=TRUE, 
              main="Enriched factors per cell", 
              distfun = "euclidean", hclustfun = "complete",
              annCol = hm.ann.col[!(colSums(binaryMat)==0),],
              annColors = ann_colors,
              color = c("white", "black"),
              cexRow = 0.5,
              labRow = FactorNames,
              labCol = rep("",ncol(binaryMat_rmOFF))
)
dev.off()
# how many cells were excluded?
sum()
