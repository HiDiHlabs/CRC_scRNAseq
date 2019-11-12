########################################################################################################
########################## 3.1 NMF of 8 LGR5+ patients (mean-centered) ###################################
########################################################################################################

library(gplots)
library(Seurat)
library(MetaDE)
library(colorspace)
library(GMD)

selected_pids <- c("HD1495-P1", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2779-P7", "HD2791-P8", "HD3254-P10", "HD3371-P11")
load("D:/Teresa/Colon-final/12patients_MC/Colon_SeuratObject_after_tSNE.Rdata")
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



write.csv(data.pos, file = "D:/Teresa/Colon-final/NNMF_LGR5/NMFinput_LGR5.csv")


## -----------------------------------------------------
## ---------- run NNMF in MATLAB -----------------------
## -----------------------------------------------------


## Select which k was used
ks <- c("k20","k25","k30")
kn <- c(20,25,30)
k=2

## Plot

nnmf.genes <- as.matrix(read.csv(paste0("D:/Teresa/Colon-final/NNMF_LGR5/",ks[k],"/nnmf_LGR5_genes_",ks[k],"_defaultParams.csv"), header=F))
nnmf.cells <- read.csv(paste0("D:/Teresa/Colon-final/NNMF_LGR5/",ks[k],"/nnmf_LGR5_cells_",ks[k],"_defaultParams.csv"), header=F)

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
  write.csv(i_save, file=sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/Genes_",ks[k],"_Factor_%s.csv"), i))
  
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
  write.csv(i_save, file=sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/Cells_",ks[k],"_Factor_%s.csv"), i)) 

}

out <- data.frame(factor = rep(1:kn[k],each=200),
           gene = names(top200),
           score = top200)
write.csv(out, file="D:/Teresa/Colon-final/NNMF_LGR5/AllFactors_top200genes.csv")

for (i in 1:kn[k]){
  if (i ==1){
    out2 <- data.frame(out$gene[out$factor==i])
  } else {
    out2 <- cbind(out2, data.frame(out$gene[out$factor==i]))
  }
}
colnames(out2) <- paste0("Factor_", as.character(1:kn[k]))
write.csv(out2, file="D:/Teresa/Colon-final/NNMF_LGR5/AllFactors_top200genes_NamesOnly.csv")

  


## ---------- Violin plots to determine if any factors are patient-specific 


require(ggplot2)
require(reshape2)
require(RColorBrewer)
library(data.table)

for (i in 1:kn[k]){
  scores.in <- as.data.table(read.csv(sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/Cells_",ks[k],"_Factor_%s.csv"), i)))
  colnames(scores.in) <- c("patient","score")
  scores.in$patient <- sub("-.*$","",scores.in$patient)
  
  
  pdf(sprintf("D:/Teresa/Colon-final/NNMF_LGR5/CellScores_perFactor_%s.pdf", i),width=18,height=6,paper='special') 
  
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
  
Overlaps_List <- list()

for (i in 1:kn[k]){
  
  scores.in <- as.data.table(read.csv(sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/Cells_",ks[k],"_Factor_%s.csv"), i)))
  colnames(scores.in) <- c("patient","score")
  scores.in$patient <- sub("-.*$","",scores.in$patient)
  
  # Plot histogram  
  pdf(sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/",ks[k],"/CellScores_perPatient_Hist_Factor_%s.pdf"), i),width=8,height=6,paper='special') 
  
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
  pdf(sprintf(paste0("D:/Teresa/Colon-final/NNMF_LGR5/",ks[k],"/CellScores_perPatient_Hist_Overlap_Factor_%s.pdf"), i),width=4,height=4,paper='special') 
  #hcol <- brewer.pal(10,"Set3")
  hist(Overlaps_RinC)
  
  dev.off() 
  
  Overlaps_List[[i]] <- Overlaps_RinC
  
}


# exclude patient-specific factors based on overlap with other patients:
# exclude if less than 50% overlap between more than 5 patient pairings
CutoffPc = 0.5

length(Overlaps_RinC)
IncludeFactors <- list()
for (o in 1:25){
  c <- colSums(Overlaps_List[[o]]<CutoffPc)
  r <- rowSums(Overlaps_List[[o]]<CutoffPc)
if (sum(c) + sum(r) > 5){
  IncludeFactors[[o]] <- F
} else {
  IncludeFactors[[o]] <- T
}
}
