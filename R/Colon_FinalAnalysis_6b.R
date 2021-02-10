###################################################################################
##### Comparison of metabolism expression in DCS and stem cells ###################
###################################################################################

# -----------Violin plot of OXPHOS expression in DCS and Stem cells ---------------

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

# load GeneSetScores
load(paste0(data.loc, "FactorScoring/GeneSetScores.Rdata"))

# load annotations
load(paste0(data.loc, "MetaData_allCells_forHeatmaps.Rdata"))
MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]

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
# and Cell_cycle (column 3) and OXPHOS (column 4) and HypoxGlycol (column 2)
isDCS <- significantGeneSets[,7]
isSTEM <- significantGeneSets[,6]
isCycling <- significantGeneSets[,3]
isOxphos <- significantGeneSets[,4]
isHypox <- significantGeneSets[,2]
isNeither <- (!isDCS)*(!isSTEM)*(!isCycling)*(!isOxphos)*(!isHypox)

sum(isDCS)
sum(isSTEM)
sum(isCycling)
sum(isOxphos)
sum(isHypox)
sum(isNeither)
# --> sigma = 1: 377 DCS, 460 STEM, 664 Cycling, 478 Oxphos, 495 Hypox, 1454 neither

# determine if cells from all 8 patients are included
PatID_long <- as.vector(MetaData_selectedCells$PatID_long)
table(PatID_long[isDCS])
table(PatID_long[isSTEM])
table(PatID_long[isCyling])
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
cellsHypox <- data.frame(
  CellName = rownames(GeneSetScores)[isHypox],
  Patient = PatID_long[isHypox],
  CellType = 'Hypox'
)
cellsOxphos <- data.frame(
  CellName = rownames(GeneSetScores)[isOxphos],
  Patient = PatID_long[isOxphos],
  CellType = 'OxPhos'
)
df <- rbind(cellsDCS,cellsSTEM,cellsCycling,cellsHypox,cellsOxphos)
rownames(df) <- 1:dim(df)[1]
# load GeneSetScores
load(paste0(data.loc, "FactorScoring/GeneSetScores_topGenesFromMergedFactors.Rdata"))
df2 <- GeneSetScores[as.character(df$CellName),]
colnames(df2) <- c("Immune_response"  ,  "HypoxGlycol" ,"Cell_cycle"      ,   "OXPHOS"  ,           "MYC"           ,     "Stem"    ,          
                   "DCS"          ,      "Fatty_acid"  )
df2 <- as.data.frame(df2)
rownames(df2) <- 1:dim(df)[1]
df <- cbind(df,df2)

# violin plots of OXPHOS and HypoxGlycol in the populations

pdf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/Factor_OXPHOS_inDCSvsSTEM.pdf"),width=7,height=5,paper='special') 
print(ggplot(data = df, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

pdf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/Factor_HypoxGlycol_inDCSvsSTEM.pdf"),width=7,height=5,paper='special') 
print(ggplot(data = df, aes(x=CellType, y=HypoxGlycol, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

# violin plots of OXPHOS and HypoxGlycol excluding quenched populations 

pdf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/Factor_OXPHOS_inDCSvsSTEM_exclQuenched.pdf"),width=6,height=5,paper='special') 
df.exclQuenched <- df[!df$CellType=='OxPhos',]
print(ggplot(data = df.exclQuenched, aes(x=CellType, y=OXPHOS, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

pdf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/Factor_HypoxGlycol_inDCSvsSTEM_exclQuenched.pdf"),width=6,height=5,paper='special') 
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
popIdx <- cbind(isDCS,
            isSTEM,
            isCycling,
            isOxphos,
            isHypox
)
popNames <- c('DCS',
             'Stem',
             'Cycling',
             'Hypox',
             'OxPhos'
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
write.csv(pvalues_OXPHOS, file = paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/pvalues_OXPHOS.csv"))
write.csv(pvalues_Hypox, file = paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/pvalues_Hypox.csv"))

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
write.csv(pvalues_ctrl_OXPHOS, file = paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/pvalues_ctrl_OXPHOS.csv"))
write.csv(pvalues_ctrl_Hypox, file = paste0(data.loc, "OXPHOS_in_DCS_vs_STEM/Factor_Factor/pvalues_ctrl_Hypox.csv"))
