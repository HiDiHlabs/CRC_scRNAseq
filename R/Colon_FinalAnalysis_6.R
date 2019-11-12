#######################################################################################################
##### Heatmap of gene set scores ###################
##### using merged factors ###################
#######################################################################################################

library(NMF)

# load Seurat object (12 patients, mean-centered)
load("D:/Teresa/Colon-final/12patients_MC/Colon_SeuratObject_after_tSNE.Rdata")
selected_pids <- c("HD1495-P1", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2779-P7", "HD2791-P8", "HD3254-P10", "HD3371-P11")
FactorNames <- c("Immune_response", "Hypoxia/Glycolysis", "Cell_cycle", "OXPHOS",
                 "MYC", "Stem", "DCS", "Fatty_acid")

# load GeneSetScores
load("D:/Teresa/Colon-final/FactorScoring/GeneSetScores_topGenesFromMergedFactors.Rdata")

# collate data for heatmap
hm.data <- GeneSetScores

# load annotations for heatmap
load("D:/Teresa/Colon-final/MetaData_allCells_forHeatmaps.Rdata")

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

hmcol = colorRampPalette(brewer.pal(11,"RdBu"))(101)
  
nmf.options(grid.patch=TRUE)

Mat <- t(GeneSetScores)
pdf("D:/Teresa/Colon-final/FactorScoring/heatmap_mergedfactors_GeneSetScores.pdf",width=20,height=9,paper='special') 
NMF::aheatmap(Mat, 
              scale="none", 
              revC=TRUE, 
              main="Enriched factors per cell", 
              distfun = "euclidean", hclustfun = "complete",
              annCol = hm.ann.col,
              annColors = ann_colors,
              color = rev(hmcol),
              cexRow = 0.5,
              labRow = FactorNames,
              labCol = rep("",ncol(Mat))
)
dev.off()



# ---------------- for LGR5-high cells only ----------------

selected_LGR5expr <- colon@scale.data["LGR5",rownames(GeneSetScores)]

# set cutoff at 1
lgr5high.idx <- selected_LGR5expr>1

# collate data for heatmap
hm.data <- GeneSetScores

# load annotations for heatmap
load("D:/Teresa/Colon-final/MetaData_allCells_forHeatmaps.Rdata")

MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]
# drop factor levels that are not in the selected data frame
MetaData_selectedCells <- lapply(MetaData_selectedCells, function(x) if(is.factor(x)) factor(x) else x)

hm.ann.col <- data.frame(MetaData_selectedCells)
hm.ann.col <- hm.ann.col[lgr5high.idx,]

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

hmcol = colorRampPalette(brewer.pal(11,"RdBu"))(101)

nmf.options(grid.patch=TRUE)

Mat <- t(GeneSetScores[lgr5high.idx,])
pdf("D:/Teresa/Colon-final/FactorScoring/heatmap_mergedfactors_GeneSetScores_LGR5highCells.pdf",width=20,height=9,paper='special') 
NMF::aheatmap(Mat, 
              scale="none", 
              revC=TRUE, 
              main="Enriched factors per cell", 
              distfun = "euclidean", hclustfun = "complete",
              annCol = hm.ann.col,
              annColors = ann_colors,
              color = rev(hmcol),
              cexRow = 0.5,
              labRow = FactorNames,
              labCol = rep("",ncol(Mat))
)
dev.off()




# -----------Violin plot of OXPHOS expression in DCS and Stem cells ---------------

require(ggplot2)
require(reshape2)
require(RColorBrewer)
library(data.table)

# binarise GeneSetScores to identify DCS and stem cells:
# set cutoff for gene set scores by sd
sigma <- 1    # 1 -> 10406 sig. sets in 3795 cells, 2 -> 2863 sig. sets in 1871 cells
cutoffGeneSetScores <- rep(0,dim(GeneSetScores)[2])
for (c in 1:length(cutoffGeneSetScores)){
  cutoffGeneSetScores[c] <- mean(GeneSetScores[,c]) + sigma*sd(GeneSetScores[,c])
}

# determine significant gene sets for each cell
significantGeneSets <- matrix(NA,dim(GeneSetScores)[1],dim(GeneSetScores)[2])
for (c in 1:length(cutoffGeneSetScores)){
  significantGeneSets[,c] <- (GeneSetScores[,c] > cutoffGeneSetScores[c])
}

# determine which cells are DCS (column 7) and STEM (column 6)
isDCS <- significantGeneSets[,7]
isSTEM <- significantGeneSets[,6]
isBOTH <- as.logical(isDCS*isSTEM)

sum(isDCS)
sum(isSTEM)
sum(isBOTH)
# --> sigma = 1: 162 DCS, 371 STEM, 49 classified as both
# --> sigma = 0.5: 781 DCS, 878 STEM, 365 classified as both

# determine if cells from all 8 patients are included
PatID_long <- as.vector(MetaData_selectedCells$PatID_long)
table(PatID_long[isDCS])
table(PatID_long[isSTEM])
table(PatID_long[isBOTH])

# record cell type for each cell: DCS/STEM/both/neither
CellType <- matrix(NA,dim(hm.data)[1],1)
for (c in 1:length(CellType)){
  if (isBOTH[c] == 1){
    CellType[c] <- "both"
  } else if (isDCS[c] == 1){
    CellType[c] <- "DCS"
  } else if (isSTEM[c] == 1){
    CellType[c] <- "STEM"
  } else {CellType[c]="neither"}
}

# plot violin plot for each signature score for DCS/STEM/both/neither cells

df <- data.frame(
  CellName = rownames(hm.data),
  Patient = MetaData_selectedCells$PatID_long,
  CellType = t(CellType),
  row.names = 1
)

df2 <- hm.data
colnames(df2) <- c("Immune_response"  ,  "HypoxGlycol" ,"Cell_cycle"      ,   "OXPHOS"  ,           "MYC"           ,     "Stem"    ,          
                   "DCS"          ,      "Fatty_acid"  )
df2 <- as.data.frame(df2)
df <- cbind(df,df2)


# select factor

i="OXPHOS"

pdf(sprintf("D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM/Factor_%s_inDCSvsSTEM.pdf", i),width=6,height=6,paper='special') 

print(ggplot(data = df, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T) +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()


# -----------Violin plots of factor expression in DCS, Stem and cycling cells (based on factor scores) ---------------

# determine which cells are cycling (column 7) and STEM (column 6)
isCycling <- significantGeneSets[,3]
isSTEMandCycling <- as.logical(isSTEM*isCycling)
isDCSandCycling <- as.logical(isDCS*isCycling)
isSTEMandDCSandCycling <- as.logical(isDCS*isSTEMandCycling)

sum(isCycling)
sum(isSTEMandCycling)
sum(isDCSandCycling)
sum(isSTEMandDCSandCycling)
# --> 1140 cycling, of which 166 STEM and 105 DCS

# determine if cells from all 8 patients are included
table(PatID_long[isCycling])
table(PatID_long[isSTEMandCycling])
table(PatID_long[isDCSandCycling])
table(PatID_long[isSTEMandDCSandCycling])

# record cell type for each cell: DCS/STEM/Cycling/neither
CellType <- matrix(NA,dim(hm.data)[1],1)
for (c in 1:length(CellType)){
  if (isSTEMandDCSandCycling[c] == 1){
    CellType[c] <- "Cycling_STEM_DCS"
  } else if (isDCSandCycling[c] == 1){
    CellType[c] <- "Cycling_DCS"
  } else if (isSTEMandCycling[c] == 1){
    CellType[c] <- "Cycling_STEM"
  } else if (isCycling[c] == 1){
    CellType[c] <- "Cycling"
  } else if (isDCS[c] == 1){
    CellType[c] <- "DCS"
  } else if (isSTEM[c] == 1){
    CellType[c] <- "STEM"
  } else {CellType[c]="neither"}
}

# plot violin plot for each signature score for DCS/STEM/both/neither cells

df <- data.frame(
  CellName = rownames(hm.data),
  Patient = PatID_long,
  CellType = t(CellType),
  row.names = 1
)

df2 <- hm.data
colnames(df2) <- c("Immune_response"  ,  "HypoxGlycol" ,"Cell_cycle"      ,   "OXPHOS"  ,           "MYC"           ,     "Stem"    ,          
                   "DCS"          ,      "Fatty_acid"  )
df2 <- as.data.frame(df2)
df <- cbind(df,df2)


# select factor

i="OXPHOS"

pdf(sprintf("D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM/Factor_%s_inDCSvsSTEMvsCycling.pdf", i),width=12,height=6,paper='special') 

print(ggplot(data = df, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T) +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()




# -----------Violin plots of factor expression in DCS, Stem and cycling cells (individually, based on factor scores) ---------------

# determine which cells are cycling (column 7) and STEM (column 6)
isCycling <- significantGeneSets[,3]
isSTEMandCycling <- as.logical(isSTEM*isCycling)
isDCSandCycling <- as.logical(isDCS*isCycling)
isSTEMandDCSandCycling <- as.logical(isDCS*isSTEMandCycling)

sum(isCycling)
sum(isSTEMandCycling)
sum(isDCSandCycling)
sum(isSTEMandDCSandCycling)
# --> 1140 cycling, of which 166 STEM and 105 DCS

# determine if cells from all 8 patients are included
table(PatID_long[isCycling])
table(PatID_long[isSTEMandCycling])
table(PatID_long[isDCSandCycling])
table(PatID_long[isSTEMandDCSandCycling])

# record cell type for each cell: DCS/STEM/Cycling/neither
CellType <- matrix(NA,dim(hm.data)[1],1)
for (c in 1:length(CellType)){
  if (isSTEMandDCSandCycling[c] == 1){
    CellType[c] <- "Cycling_STEM_DCS"
  } else if (isDCSandCycling[c] == 1){
    CellType[c] <- "Cycling_DCS"
  } else if (isSTEMandCycling[c] == 1){
    CellType[c] <- "Cycling_STEM"
  } else if (isBOTH[c] == 1){
    CellType[c] <- "STEM_DCS"
  } else if (isCycling[c] == 1){
    CellType[c] <- "Cycling"
  } else if (isDCS[c] == 1){
    CellType[c] <- "DCS"
  } else if (isSTEM[c] == 1){
    CellType[c] <- "STEM"
  } else {CellType[c]="neither"}
}

# plot violin plot for each signature score for DCS/STEM/both/neither cells

df <- data.frame(
  CellName = rownames(hm.data),
  Patient = PatID_long,
  CellType = CellType,
  row.names = 1
)

df2 <- hm.data
colnames(df2) <- c("Immune_response"  ,  "HypoxGlycol" ,"Cell_cycle"      ,   "OXPHOS"  ,           "MYC"           ,     "Stem"    ,          
                   "DCS"          ,      "Fatty_acid"  )
df2 <- as.data.frame(df2)
df <- cbind(df,df2)


# select factor

i="HypoxGlycol"

pdf(sprintf("D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM/Factor_%s_inDCSvsSTEMvsCycling.pdf", i),width=12,height=6,paper='special') 

print(ggplot(data = df, aes(x=CellType, y=HypoxGlycol, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()



# -----------Violin plots of factor expression in DCS, Stem and cycling cells (individually, based on top 200 genes from 8 merged factors) ---------------

# Load merged factor gene lists - use merger of top genes for merged factors
load("D:/Teresa/Colon-final/MergedFactors_GeneLists.Rdata")

df <- data.frame(
  CellName = rownames(hm.data),
  Patient = PatID_long,
  CellType = t(CellType),
  row.names = 1
)


topGenes.data <- matrix(NA, nrow = dim(hm.data)[1], ncol =length(MergedFactors_GeneLists))
for (f in 1:length(MergedFactors_GeneLists)){
  for (c in 1:dim(hm.data)[1]){
    topGenes.data[c,f] <- mean(log.cpm.matrix.mc[MergedFactors_GeneLists[[f]],c])
  }
}
df2 <- topGenes.data
colnames(df2) <- c("Immune_response"  ,  "HypoxGlycol" ,"Cell_cycle"      ,   "OXPHOS"  ,           "MYC"           ,     "Stem"    ,          
                   "DCS"          ,      "Fatty_acid"  )
df2 <- as.data.frame(df2)
df <- cbind(df,df2)


# select factor

i="HypoxGlycol"

pdf(sprintf("D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM/topGenes_Factor_%s_inDCSvsSTEMvsCycling.pdf", i),width=12,height=6,paper='special') 

print(ggplot(data = df, aes(x=CellType, y=HypoxGlycol, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()


# -----------Violin plots of factor expression in DCS, Stem and cycling cells (individually, based on top 200 genes from 8 merged factors) ---------------

# Load unmerged factor gene lists
CuratedGenes <- read.table("D:/Teresa/Colon-final/NNMF_LGR5/CuratedFactorList.txt",header = T)[1:200,]

df <- data.frame(
  CellName = rownames(hm.data),
  Patient = PatID_long,
  CellType = CellType,
  row.names = 1
)


topGenes.data <- matrix(NA, nrow = dim(hm.data)[1], ncol = dim(CuratedGenes)[2])
for (f in 1:dim(CuratedGenes)[2]){
  for (c in 1:dim(hm.data)[1]){
    topGenes.data[c,f] <- mean(log.cpm.matrix.mc[as.character(CuratedGenes[,f]),c])
  }
}
df2 <- topGenes.data
colnames(df2) <- colnames(CuratedGenes)
df2 <- as.data.frame(df2)
df <- cbind(df,df2)


# select factor

i="Factor_23"

pdf(sprintf("D:/Teresa/Colon-final/OXPHOS_in_DCS_vs_STEM/unmergedTopGenes_%s_inDCSvsSTEMvsCycling.pdf", i),width=12,height=6,paper='special') 

print(ggplot(data = df, aes(x=CellType, y=Factor_23, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()




# ----------- t-test of OXPHOS and HypoxGlycol in DCS vs. STEM---------------

t.test(,)
t.test(df$OXPHOS[isDCS],df$OXPHOS[isSTEM])

dcs <- data.frame(HypoxGlycol = df$HypoxGlycol[isDCS],
                  OXPHOS = df$OXPHOS[isDCS])

stem <- data.frame(HypoxGlycol = df$HypoxGlycol[isSTEM],
                  OXPHOS = df$OXPHOS[isSTEM])

dcs$celltype <- 'DCS'
stem$celltype <- 'STEM'

toplot <- rbind(dcs, stem)

ggplot(toplot, aes(HypoxGlycol, fill = celltype)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
ggplot(toplot, aes(OXPHOS, fill = celltype)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')




# ----------- Does LGR5 correlate with OXPHOS? ---------------

# number of cells in the different cell types identified
cutoff_low <- 0
cutoff_high <- 2.5

# select same number of cells with highest expression of cell type markers
isLGR5high <- (log.cpm.matrix.mc["LGR5",]>cutoff_low)
isDEFA5high <- (log.cpm.matrix.mc["DEFA5",]>cutoff_low)
isTOP2Ahigh <- (log.cpm.matrix.mc["TOP2A",]>cutoff_high)
isMKI67high <- (log.cpm.matrix.mc["MKI67",]>cutoff_high)
isKRT20high <- (log.cpm.matrix.mc["KRT20",]>cutoff_high)
isTFF3high <- (log.cpm.matrix.mc["TFF3",]>cutoff_high)
isPRDX3high <- (log.cpm.matrix.mc["PRDX6",]>cutoff_high)

lgr5high <- data.frame(HypoxGlycol = df$HypoxGlycol[isLGR5high],
                  OXPHOS = df$OXPHOS[isLGR5high])

defa5high <- data.frame(HypoxGlycol = df$HypoxGlycol[isDEFA5high],
                   OXPHOS = df$OXPHOS[isDEFA5high])

top2ahigh <- data.frame(HypoxGlycol = df$HypoxGlycol[isTOP2Ahigh],
                    OXPHOS = df$OXPHOS[isTOP2Ahigh])

mki67high <- data.frame(HypoxGlycol = df$HypoxGlycol[isMKI67high],
                    OXPHOS = df$OXPHOS[isMKI67high])

lgr5high$celltype <- 'lgr5high'
defa5high$celltype <- 'defa5high'
top2ahigh$celltype <- 'top2ahigh'
mki67high$celltype <- 'mki67high'


toplot <- rbind(lgr5high, defa5high, top2ahigh, mki67high)
toplot <- rbind(lgr5high, defa5high)

ggplot(toplot, aes(HypoxGlycol, fill = celltype)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
ggplot(toplot, aes(OXPHOS, fill = celltype)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')

t.test(lgr5high$HypoxGlycol,defa5high$HypoxGlycol)
t.test(lgr5high$OXPHOS,defa5high$OXPHOS)

# as control: randomly sampled cells

rndctrl1 <- sample(1:dim(log.cpm.matrix.mc)[1], 600)
rndctrl2 <- sample(1:dim(log.cpm.matrix.mc)[1], 600)

ctrl1 <- data.frame(HypoxGlycol = df$HypoxGlycol[rndctrl1],
                    OXPHOS = df$OXPHOS[rndctrl1],
                    celltype = 'ctrl1')

ctrl2 <- data.frame(HypoxGlycol = df$HypoxGlycol[rndctrl2],
                    OXPHOS = df$OXPHOS[rndctrl2],
                    celltype = 'ctrl2')

t.test(ctrl1$HypoxGlycol,ctrl2$HypoxGlycol)
t.test(ctrl1$OXPHOS,ctrl2$OXPHOS)



