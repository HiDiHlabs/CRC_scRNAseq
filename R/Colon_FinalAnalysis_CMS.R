##############################################################################
############### CMS analysis of colon data for each patient ##################
##############################################################################

### ---------------- for each cell individually ------------------------------

library(Biobase)
library(CMScaller)

pids <- c("HD1495-P1", "HD1509-P2", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2596-P6", "HD2779-P7", "HD2791-P8", "HD3192-P9", "HD3254-P10", "HD3371-P11", "POP1-P12")


for (i in 1:12){
  
  folderpath <- paste0("D:/Teresa/Colon-final/DEA_PatientsSeparately/",pids[i],"/")
  load(paste0(folderpath, "patData_SeuratObject_after_tSNE.Rdata"))
  
  # get RNA-seq counts from data
  counts <- as.matrix(patData@data)
  counts <- t(scale(t(counts), center = TRUE, scale = TRUE))
  
  # CMS prediction 
  res <- CMScaller(emat=counts, rowNames = "symbol", FDR=0.25, doPlot=T)
  res[which(!is.na(res$prediction)),]
  table(res$prediction)
  
  r <- c( table(res$prediction), (dim(patData@data)[2]-sum(table(res$prediction))))
  names(r) <- c(names(r)[1:4], "NA")
  
  if (i==1){
    CMS_predictions <- r
  } else{
    CMS_predictions <- rbind(CMS_predictions, r)
  }
  
}

rownames(CMS_predictions) <- pids
write.csv(CMS_predictions, file="D:/Teresa/Colon-final/CMS_predictions.csv")



### ---------------- pooled counts for each patient ------------------------------

library(Biobase)
library(CMScaller)

pids <- c("HD1495-P1", "HD1509-P2", "HD1664-P3", "HD1883-P4", "HD1960-P5", "HD2596-P6", "HD2779-P7", "HD2791-P8", "HD3192-P9", "HD3254-P10", "HD3371-P11", "POP1-P12")

folderpath <- "D:/Teresa/Colon-final/12patients_noMC/"
load(paste0(folderpath, "Colon_SeuratObject_after_tSNE.Rdata"))

colon@meta.data$PatID <- as.character(sub("\\-.*","",colnames(colon@data)))
PatIDs <- unique(colon@meta.data$PatID)

for (i in 1:length(PatIDs)){
  
  idx <- colon@meta.data$PatID == PatIDs[i]
  counts.add <- rowSums(colon@data[,idx])
  
  if (i==1){
    counts <- counts.add
  } else {
    counts <- cbind(counts, counts.add)
  }
}

colnames(counts) <- PatIDs

# CMS prediction 
res <- CMScaller(emat=counts, rowNames = "symbol", FDR=0.25, doPlot=T)
res[which(!is.na(res$prediction)),]
table(res$prediction)


write.csv(res, file="D:/Teresa/Colon-final/CMS_predictions_pooledCounts.csv")
