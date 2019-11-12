# ----------- determine fraction of cells expressing more than one gene set ---------- 

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
      # calculate proportion of cells that express both gene sets i and j  (denominator: all cells from this patient)
      OverlapArray[p,i,j] <- sum(significantGeneSets[cells.to.inspect,i]*significantGeneSets[cells.to.inspect,j])/
        length(cells.to.inspect)
    }
  }
  pdf(sprintf("D:/Teresa/Colon-final/FactorScoring/Heatmap_GeneSetOverlap_Pat_%s_mergedFactors.pdf", selected_pids[p]),width=6.5,height=6,paper='special') 
  aheatmap(OverlapArray[p,,], 
           Rowv=T,
           Colv="Rowv",
           scale="none",
           main=selected_pids[p],
           col=cmap,
           labRow=FactorNames,
           labCol=FactorNames,
           fontsize=12,
           cexRow=0.9,
           cexCol=0.9,
           legend=T,
           txt=matrix(as.character(round(100*OverlapArray[p,,], digits = 0)),ncol=8,nrow=8)
  )
  dev.off()
}

