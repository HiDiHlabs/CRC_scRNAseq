########################################################################################################
########################## 3.2 NMF of 8 Liver metastases (mean-centered) ###################################
########################################################################################################

library(gplots)
library(Seurat)
library(MetaDE)
library(colorspace)
library(GMD)

selected_pids <- c("P1", "P3", "P4", "P5", "P7", "P8")
idx <- colon@meta.data$orig.ident %in% sub("\\-.*","",selected_pids)
colon.v <- SubsetData(colon, cells.use = colon@cell.names[idx])
colon.v <- FindVariableGenes(colon.v, mean.function = ExpMean,
                             dispersion.function = LogVMR, do.plot = TRUE, set.var.genes = TRUE,
                             x.low.cutoff = 0, x.high.cutoff = Inf, 
                             y.cutoff = 0.5,
                             y.high.cutoff = Inf)
length(colon.v@var.genes)

data <- colon.v@data
data.pos <- data
data.pos[data.pos < 0] <- 0
# res <- nmf(data.pos, 25) -- very slow in R!

write.csv(data.pos, file = paste0(data.loc, "NMFinput_Liver.csv"))

## ----------------------------------------------------
## ---------- run NNMF in MATLAB -----------------------
## -----------------------------------------------------

## Select which k was used
ks <- c("k15","k20","k25")
kn <- c(15,20,25)

k=3

nnmf.genes <- as.matrix(read.csv(paste0(data.loc, "NNMF_Liver/nnmf_Liver_genes_",ks[k],"_defaultParams.csv"), header=F))
nnmf.cells <- read.csv(paste0(data.loc, "NNMF_Liver/nnmf_Liver_cells_",ks[k],"_defaultParams.csv"), header=F)

n <- as.matrix(nnmf.genes)

clustfun <- function(x){
  h <- hclust(x, method = "ward.D2")
  return(h)
}

hmcol<-diverge_hsv(101)
heatmap.2(log(n+1),
          Rowv = T,  
          Colv = T,
          dendrogram = "both",
          scale="none",
          trace="none",
          density.info="none",
          col=hmcol, #PurpleAndYellow(),
          # RowSideColors=col1[gr.row], 
          # ColSideColors=col2[gr.col],
          # colsep = NULL,
          # rowsep = row.cuts,
          # sepcolor="white",
          # sepwidth=c(0.1,0.1),
          # key=T,
          cexRow = 0.1,
          #labCol = T, #rep("",dim(data_to_plot)[2]),
          cexCol= 1,
          margins=c(5,5),
          hclustfun = clustfun
)

## ---------- Save lists of genes & cells sorted by representation in NNMF factors -----------

gene.list <- rownames(colon.v@data)
cell.list <- colnames(colon.v@data)

top200 <- NULL
for (i in 1:kn[k]){
  # genes
  i_order <- sort(nnmf.genes[,i], index.return = T)
  i_save <- rev(i_order$x)
  names(i_save) <- gene.list[rev(i_order$ix)]
  write.csv(i_save, file=sprintf(paste0(data.loc, "NNMF_Liver/Genes_",ks[k],"_Factor_%s.csv"), i))
  
  # store top 200 genes per factor in separate matrix
  if (i==1){
    top200 <- i_save[1:200]
  } else {
    top200 <- c(top200, i_save[1:200])
  }
  
  # cells
  i_order <- sort(nnmf.cells[,i], index.return = T)
  i_save <- rev(i_order$x)
  names(i_save) <- cell.list[rev(i_order$ix)]
  write.csv(i_save, file=sprintf(paste0(data.loc, "NNMF_Liver/Cells_",ks[k],"_Factor_%s.csv"), i)) 
}

out <- data.frame(factor = rep(1:kn[k],each=200),
                  gene = names(top200),
                  score = top200)
write.csv(out, file=paste0(data.loc, "NNMF_Liver/AllFactors_top200genes.csv"))

## ---------- Violin plots to determine if any factors are patient-specific 

require(ggplot2)
require(reshape2)
require(RColorBrewer)
library(data.table)

for (i in 1:kn[k]){
  scores.in <- as.data.table(read.csv(sprintf(paste0(data.loc, "NNMF_Liver/Cells_",ks[k],"_Factor_%s.csv"), i)))
  colnames(scores.in) <- c("patient","score")
  scores.in$patient <- sub("-.*$","",scores.in$patient)
  
  
  pdf(sprintf(paste0(data.loc, "NNMF_Liver/CellScores_perFactor_%s.pdf"), i),width=18,height=6,paper='special') 
  
  print(ggplot(data = scores.in, aes(x=patient, y=score, fill=patient)) + 
          geom_violin(trim=T) +
          geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
          geom_boxplot(width=0.1, size=1.2) +
          scale_fill_brewer(palette="Set3") +
          scale_color_brewer(palette="Set3") +
          theme(legend.position="none"))
  
  dev.off()
}

## --------- Cell score distributions for each factor and patient --> overlap -------------------

require(reshape2)
require(RColorBrewer)
library(data.table)

pids_short <- sub("-.*$","",selected_pids)

for (i in 1:kn[k]){
  
  scores.in <- as.data.table(read.csv(sprintf(paste0(data.loc, "NNMF_Liver/Cells_",ks[k],"_Factor_%s.csv"), i)))
  colnames(scores.in) <- c("patient","score")
  scores.in$patient <- sub("-.*$","",scores.in$patient)
  
  # Plot histogram  
  pdf(sprintf(paste0(data.loc, "NNMF_Liver/",ks[k],"/CellScores_perPatient_Hist_Factor_%s.pdf"), i),width=8,height=6,paper='special') 
  
  print(ggplot(scores.in,aes(x=score, fill=patient)) + 
          geom_histogram(data=subset(scores.in,patient == pids_short[1]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[2]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[3]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[4]), col="grey", alpha = 0.3, bins = 25, position="dodge") +  
          geom_histogram(data=subset(scores.in,patient == pids_short[5]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[6]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[7]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          geom_histogram(data=subset(scores.in,patient == pids_short[8]), col="grey", alpha = 0.3, bins = 25, position="dodge") +
          scale_fill_brewer(palette="Set3"))  
  
  dev.off()
  
  # Calculate overlap
  Overlaps_RinC <- matrix(NA, ncol=8, nrow=8)
  for (p in 1:8){
    for (q in 1:8){
      p_scores <- scores.in$score[scores.in$patient==pids_short[p]]
      q_scores <- scores.in$score[scores.in$patient==pids_short[q]]
      q_min <- min(q_scores)
      q_max <- max(q_scores)
      Overlaps_RinC[p,q] <- sum(p_scores>=q_min & p_scores<=q_max)/length(p_scores)   # p in q
    }
  }
  
  # Plot histogram of overlaps
  pdf(sprintf(paste0(data.loc, "NNMF_Liver/",ks[k],"/CellScores_perPatient_Hist_Overlap_Factor_%s.pdf"), i),width=4,height=4,paper='special') 
  #hcol <- brewer.pal(10,"Set3")
  hist(Overlaps_RinC)
  dev.off() 
}
