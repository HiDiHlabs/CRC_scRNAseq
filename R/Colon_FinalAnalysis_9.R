###################################################################################
##### Define cell types vs. cell states ###################
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

# load GeneSetScores list object
load(paste0(data.loc, "FactorScoring/GeneSetScores_topGenesFromMergedFactors.Rdata"))

# load gene lists defining cell types and cell states
Definition_CellTypes_in <- read.table(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/Definition_CellTypes.csv"), header = T, sep = ";")
Definition_CellStates_in <- read.table(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/Definition_CellStates.csv"), header = T, sep = ";")    
# turn into list objects
Definition_CellTypes_GeneList <- list()
for (c in 1:dim(Definition_CellTypes_in)[2]){
  Definition_CellTypes_GeneList[[c]] <- as.character(Definition_CellTypes_in[,c])
}
save(Definition_CellTypes_GeneList, file=paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/Definition_CellTypes_GeneList.Rdata"))
Definition_CellStates_GeneList <- list()
for (c in 1:dim(Definition_CellStates_in)[2]){
  Definition_CellStates_GeneList[[c]] <- as.character(Definition_CellStates_in[,c])[1:(min(c((which(Definition_CellStates_in[[c]]=="")-1),length(Definition_CellStates_in[[c]]))))]
}
save(Definition_CellStates_GeneList, file=paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/Definition_CellStates_GeneList.Rdata"))
  
# load annotations
load(paste0(data.loc, "MetaData_allCells_forHeatmaps.Rdata"))
MetaData_selectedCells <- MetaData_allCells[sub("\\-.*","",colon@cell.names) %in% sub("\\-.*","",selected_pids),]

#########  run Colon_FinalAnalysis_9b to obtain GeneSetScores matrices for cell types & cell states ###########  

# -------------- 1: Cell type identification based on cell type genes -------------- 

load(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/DefinedCellTypes_GeneSetScores.Rdata"))

# binarise GeneSetScores to identify cell types:
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
isStem <- significantGeneSets[,1]
isDCS <- significantGeneSets[,2]
isTA <- significantGeneSets[,3]
isTdiff <- significantGeneSets[,4]
isNeither <- (!isStem)*(!isDCS)*(!isTA)*(!isTdiff)

# -------------- 2: combine all identified cells into one data frame -------------- 

sum(isStem)
sum(isDCS)
sum(isTA)
sum(isTdiff)
sum(isNeither)

# determine if cells from all 8 patients are included
PatID_long <- as.vector(MetaData_selectedCells$PatID_long)
table(PatID_long[isStem])
table(PatID_long[isDCS])
table(PatID_long[isTA])
table(PatID_long[isTdiff])

# assemble data frame with all cells (N.B. cells that are classified as >1+ will appear >once)
cellsStem <- data.frame(
  CellName = rownames(GeneSetScores)[isStem],
  Patient = PatID_long[isStem],
  CellType = 'Stem'
)
cellsDCS <- data.frame(
  CellName = rownames(GeneSetScores)[isDCS],
  Patient = PatID_long[isDCS],
  CellType = 'DCS'
)
cellsTA <- data.frame(
  CellName = rownames(GeneSetScores)[isTA],
  Patient = PatID_long[isTA],
  CellType = 'TA'
)
cellsTdiff <- data.frame(
  CellName = rownames(GeneSetScores)[isTdiff],
  Patient = PatID_long[isTdiff],
  CellType = 'Tdiff'
)
df1 <- rbind(cellsStem,cellsDCS,cellsTA,cellsTdiff)
rownames(df1) <- 1:dim(df1)[1]


# -------------- 3: Add scores for cell states --------------

load(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/DefinedCellStates_GeneSetScores.Rdata"))

df2 <- GeneSetScores[as.character(df1$CellName),]
colnames(df2) <- c("CellCycle" , "OXPHOS"  ,   "Hypoxia"   , "Glycolysis" ,"FattyAcid"  )

df2 <- as.data.frame(df2)
rownames(df2) <- 1:dim(df1)[1]
df <- cbind(df1,df2)

# ---------------------------- 4: Violin plots and t-tests ----------------------------

# violin plots of cell states

writefolder <- "CellTypes_CellStates/FattyAcid"
#  "CellCycle"  "OXPHOS"     "Hypoxia"    "Glycolysis" "FattyAcid" 

pdf(sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/CellState_in_CellTypes.pdf"),writefolder),width=7,height=5,paper='special') 
print(ggplot(data = df, aes(x=CellType, y=FattyAcid, fill=CellType)) + 
        geom_violin(trim=T, scale='count') +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))
dev.off()

# pair-wise t-tests of cell state expression between the populations
popNames <- c('Stem',
              'DCS',
              'TA',
              'Tdiff'
)

pvalues <- rep(NA, choose(length(popNames),2))
i=1
for (a in 2:length(popNames)){
  for (b in 1:(a-1)){
    pvalues[i] <- t.test(df$FattyAcid[df$CellType==popNames[a]],df$FattyAcid[df$CellType==popNames[b]])[[3]]
    names(pvalues)[i] <- paste(popNames[a],popNames[b],sep='_')
    i = i+1
  }
}
write.csv(pvalues, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues.csv"),writefolder))


# control: random sets of cells
pvalues_ctrl <- rep(NA, choose(length(popNames),2))
i=1
for (a in 2:length(popNames)){
  for (b in 1:(a-1)){
    ctrl1 <- sample(1:dim(df)[1],sum(df$CellType==popNames[a]))
    ctrl2 <- sample(1:dim(df)[1],sum(df$CellType==popNames[b]))
    pvalues_ctrl[i] <- t.test(df$FattyAcid[ctrl1],df$FattyAcid[ctrl2])[[3]]
    names(pvalues_ctrl)[i] <- paste(popNames[a],popNames[b],sep='_')
    i = i+1
  }
}
write.csv(pvalues_ctrl, file = sprintf(paste0(data.loc, "OXPHOS_in_DCS_vs_STEM_new/%s/pvalues_ctrl.csv"),writefolder))
