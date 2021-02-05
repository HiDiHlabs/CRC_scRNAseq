###################################################################################
##### Comparison of metabolism expression in DCS and stem cells ###################
###################################################################################

require(ggplot2)
require(reshape2)
require(RColorBrewer)
library(data.table)
library(NMF)

data.loc = ### INSERT DATA FOLDER HERE ###

# load Seurat object (12 patients, mean-centered)
load(paste0(data.loc, "12patients_MC/Colon_SeuratObject_after_tSNE.Rdata"))
selected_pids <- c("P1", "P3", "P4", "P5", "P7", "P8", "P10", "P11")
FactorNames <- c("Immune_response", "Hypoxia/Glycolysis", "Cell_cycle", "OXPHOS",
                 "MYC", "Stem", "DCS", "Fatty_acid")
log.cpm.matrix.mc <- colon@raw.data[,sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids)]

load(paste0(data.loc, "FactorScoring/GeneSetScores_topGenesFromMergedFactors.Rdata"))

# load annotations
load(paste0(data.loc, "MetaData_allCells_forHeatmaps.Rdata"))
MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]

# -------------- 1a: Cell type identification based on factor scores -------------- 

# binarise GeneSetScores to identify DCS and stem cells:
# set cutoff for gene set scores by sd
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

# determine which cells are DCS (column 7) and STEM (column 6) 
# and Cell_cycle (column 3) and OXPHOS (column 4) and HypoxGlycol (column 2)
isDCS <- significantGeneSets[,7]
isSTEM <- significantGeneSets[,6]
isCycling <- significantGeneSets[,3]
isOxphos <- significantGeneSets[,4]
isHypox <- significantGeneSets[,2]
isNeither <- (!isDCS)*(!isSTEM)*(!isCycling)*(!isOxphos)*(!isHypox)

# -------------- 1b: Cell type identification based on marker genes -------------- 

cutoff_low <- 0
cutoff_high <- 2.5

# select cells with highest expression of cell type markers
isLGR5high <- (log.cpm.matrix.mc["LGR5",]>cutoff_low)
isDEFA5high <- (log.cpm.matrix.mc["DEFA5",]>cutoff_low)
isTOP2Ahigh <- (log.cpm.matrix.mc["TOP2A",]>cutoff_high)
isMKI67high <- (log.cpm.matrix.mc["MKI67",]>cutoff_high)
isKRT20high <- (log.cpm.matrix.mc["KRT20",]>cutoff_high)
isTFF3high <- (log.cpm.matrix.mc["TFF3",]>cutoff_high)
isPRDX3high <- (log.cpm.matrix.mc["PRDX6",]>cutoff_high)

isDCS <- isDEFA5high
isSTEM <- isLGR5high
isCycling <- isTOP2Ahigh
isOxphos <- isPRDX3high
isHypox <- isKRT20high
isNeither <- (!isDCS)*(!isSTEM)*(!isCycling)*(!isOxphos)*(!isHypox)

# -------------- 2: combine all identified cells into one data frame -------------- 

sum(isDCS)
sum(isSTEM)
sum(isCycling)
sum(isOxphos)
sum(isHypox)
sum(isNeither)

# determine if cells from all 8 patients are included
PatID_long <- as.vector(MetaData_selectedCells$PatID_long)
table(PatID_long[isDCS])
table(PatID_long[isSTEM])
table(PatID_long[isCycling])
table(PatID_long[isOxphos])
table(PatID_long[isHypox])

# assemble data frame with all cells (N.B. cells that are classified as >1+ will appear >once)
cellsDCS <- data.frame(
  CellName = rownames(GeneSetScores)[isDCS],
  Patient = PatID_long[isDCS],
  CellType = 'DCS'
)
cellsSTEM <- data.frame(
  CellName = rownames(GeneSetScores)[isSTEM],
  Patient = PatID_long[isSTEM],
  CellType = 'Stem'
)
cellsCycling <- data.frame(
  CellName = rownames(GeneSetScores)[isCycling],
  Patient = PatID_long[isCycling],
  CellType = 'Cycling'
)
cellsOxphos <- data.frame(
  CellName = rownames(GeneSetScores)[isOxphos],
  Patient = PatID_long[isOxphos],
  CellType = 'Oxphos'
)
cellsHypox <- data.frame(
  CellName = rownames(GeneSetScores)[isHypox],
  Patient = PatID_long[isHypox],
  CellType = 'Hypox'
)
df1 <- rbind(cellsDCS,cellsSTEM,cellsCycling,cellsOxphos,cellsHypox)
rownames(df1) <- 1:dim(df1)[1]


# -------------- 3a: Add expression data from factor scores --------------

df2 <- GeneSetScores[as.character(df1$CellName),c(7,6,3,4,2)]
colnames(df2) <- c("DCS"      , "STEM",   "Cycling",   "OXPHOS"  ,   "HypoxGlycol"  )

df2 <- as.data.frame(df2)
rownames(df2) <- 1:dim(df1)[1]
df <- cbind(df1,df2)



# -------------- 3b: Add expression data from top 200 genes --------------

# load merged factor gene lists 
load(paste0(data.loc, "MergedFactors_GeneLists.Rdata"))

topGenes.data <- matrix(NA, nrow = ncol(log.cpm.matrix.mc), ncol = length(MergedFactors_GeneLists))
for (f in 1:length(MergedFactors_GeneLists)){
  for (c in 1:ncol(log.cpm.matrix.mc)){
    topGenes.data[c,f] <- mean(log.cpm.matrix.mc[MergedFactors_GeneLists[[f]],c])
  }
}
rownames(topGenes.data) <- colnames(log.cpm.matrix.mc)


df2 <- topGenes.data[as.character(df1$CellName),c(7,6,3,4,2)]
colnames(df2) <- c("DCS"      , "STEM",   "Cycling",   "OXPHOS"  ,   "HypoxGlycol"  )

df2 <- as.data.frame(df2)
rownames(df2) <- 1:dim(df1)[1]
df <- cbind(df1,df2)



# ---------------------------- 4: Violin plots and t-tests ----------------------------

# violin plots of OXPHOS and HypoxGlycol in the populations

writefolder <- "MarkerGenes_MeanExpression"

pdf(sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/Factor_OXPHOS_inDCSvsSTEM.pdf"),writefolder),width=7,height=5,paper='special') 
print(ggplot(data = df, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

pdf(sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/Factor_HypoxGlycol_inDCSvsSTEM.pdf"),writefolder),width=7,height=5,paper='special') 
print(ggplot(data = df, aes(x=CellType, y=HypoxGlycol, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()


# violin plots of OXPHOS and HypoxGlycol excluding quenched populations 

pdf(sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/Factor_OXPHOS_inDCSvsSTEM_exclQuenched.pdf"),writefolder),width=6,height=5,paper='special') 
df.exclQuenched <- df[!df$CellType=='Oxphos',]
print(ggplot(data = df.exclQuenched, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

pdf(sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/Factor_HypoxGlycol_inDCSvsSTEM_exclQuenched.pdf"),writefolder),width=6,height=5,paper='special') 
df.exclQuenched <- df[!df$CellType=='Hypox',]
print(ggplot(data = df.exclQuenched, aes(x=CellType, y=HypoxGlycol, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

# pair-wise t-tests of OXPHOS and HypoxGlycol expression between the populations
popNames <- c('DCS',
              'Stem',
              'Cycling',
              'Hypox',
              'Oxphos'
)

pvalues_OXPHOS <- rep(NA, choose(length(popNames),2))
pvalues_Hypox <- rep(NA, choose(length(popNames),2))
i=1
for (a in 2:length(popNames)){
  for (b in 1:(a-1)){
    pvalues_OXPHOS[i] <- t.test(df$OXPHOS[df$CellType==popNames[a]],df$OXPHOS[df$CellType==popNames[b]])[[3]]
    pvalues_Hypox[i] <- t.test(df$HypoxGlycol[df$CellType==popNames[a]],df$HypoxGlycol[df$CellType==popNames[b]])[[3]]
    names(pvalues_OXPHOS)[i] <- paste(popNames[a],popNames[b],sep='_')
    names(pvalues_Hypox)[i] <- paste(popNames[a],popNames[b],sep='_')
    i = i+1
  }
}
write.csv(pvalues_OXPHOS, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues_OXPHOS.csv"),writefolder))
write.csv(pvalues_Hypox, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues_Hypox.csv"),writefolder))


# control: random sets of cells
pvalues_ctrl_OXPHOS <- rep(NA, choose(length(popNames),2))
pvalues_ctrl_Hypox <- rep(NA, choose(length(popNames),2))
i=1
for (a in 2:length(popNames)){
  for (b in 1:(a-1)){
    ctrl1 <- sample(1:dim(df)[1],sum(df$CellType==popNames[a]))
    ctrl2 <- sample(1:dim(df)[1],sum(df$CellType==popNames[b]))
    pvalues_ctrl_OXPHOS[i] <- t.test(df$OXPHOS[ctrl1],df$OXPHOS[ctrl2])[[3]]
    pvalues_ctrl_Hypox[i] <- t.test(df$HypoxGlycol[ctrl1],df$HypoxGlycol[ctrl2])[[3]]
    names(pvalues_ctrl_OXPHOS)[i] <- paste(popNames[a],popNames[b],sep='_')
    names(pvalues_ctrl_Hypox)[i] <- paste(popNames[a],popNames[b],sep='_')
    i = i+1
  }
}
write.csv(pvalues_ctrl_OXPHOS, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues_ctrl_OXPHOS.csv"),writefolder))
write.csv(pvalues_ctrl_Hypox, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues_ctrl_Hypox.csv"),writefolder))

# ---- check overlap between cells selected based on factor scores vs. marker genes ----

isDCS.Factor <- isDCS
isSTEM.Factor <- isSTEM
isCycling.Factor <- isCycling
isOxphos.Factor <- isOxphos
isHypox.Factor <- isHypox
isNeither.Factor <- isNeither

isDCS.MarkerGenes <- isDCS
isSTEM.MarkerGenes <- isSTEM
isCycling.MarkerGenes <- isCycling
isOxphos.MarkerGenes <- isOxphos
isHypox.MarkerGenes <- isHypox
isNeither.MarkerGenes <- isNeither

2*sum(isDCS.MarkerGenes*isDCS.Factor)/(sum(isDCS.MarkerGenes)+sum(isDCS.Factor))
2*sum(isSTEM.MarkerGenes*isSTEM.Factor)/(sum(isSTEM.MarkerGenes)+sum(isSTEM.Factor))
2*sum(isCycling.MarkerGenes*isCycling.Factor)/(sum(isCycling.MarkerGenes)+sum(isCycling.Factor))
2*sum(isOxphos.MarkerGenes*isOxphos.Factor)/(sum(isOxphos.MarkerGenes)+sum(isOxphos.Factor))
2*sum(isHypox.MarkerGenes*isHypox.Factor)/(sum(isHypox.MarkerGenes)+sum(isHypox.Factor))
2*sum(isNeither.MarkerGenes*isNeither.Factor)/(sum(isNeither.MarkerGenes)+sum(isNeither.Factor))

sum(isDCS.MarkerGenes*isDCS.Factor)/min(sum(isDCS.MarkerGenes),sum(isDCS.Factor))
sum(isSTEM.MarkerGenes*isSTEM.Factor)/min(sum(isSTEM.MarkerGenes),sum(isSTEM.Factor))
sum(isCycling.MarkerGenes*isCycling.Factor)/min(sum(isCycling.MarkerGenes),sum(isCycling.Factor))
sum(isOxphos.MarkerGenes*isOxphos.Factor)/min(sum(isOxphos.MarkerGenes),sum(isOxphos.Factor))
sum(isHypox.MarkerGenes*isHypox.Factor)/min(sum(isHypox.MarkerGenes),sum(isHypox.Factor))
sum(isNeither.MarkerGenes*isNeither.Factor)/min(sum(isNeither.MarkerGenes),sum(isNeither.Factor))
