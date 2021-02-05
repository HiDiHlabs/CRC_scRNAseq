######################################################
########### Reads / Genes per Cell #####################
######################################################

folderpath <- ### INSERT DATA FOLDER HERE ###
count.matrix <- read.csv(paste(folderpath, "rawdata/AllPatients_counts_noMeanCen_noLogTrans_incl12.csv", sep=""),
                         row.names = "X")
  
# ------------- reads per cell -------------------------------------

cellSums <- colSums(count.matrix)

cellSums_dt <- as.data.table(cellSums)
colnames(cellSums_dt) <- "reads"
rownames(cellSums_dt) <- colnames(count.matrix)
cellSums_dt$patient <- sub("\\..*$","",rownames(cellSums_dt))
cellSums_dt$log_reads <- log(cellSums_dt$reads)

pdf(sprintf(paste0(data.loc, 'detectedReads.pdf'), i),width=18,height=6,paper='special') 

print(ggplot(data = cellSums_dt, aes(x=patient, y=reads, fill=patient)) + 
        geom_violin(trim=T) +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()

pdf(sprintf(paste0(data.loc, 'detectedReads_log.pdf', i),width=18,height=6,paper='special') 

print(ggplot(data = cellSums_dt, aes(x=patient, y=log_reads, fill=patient)) + 
        geom_violin(trim=T) +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()

# average across patients
mean(cellSums_dt$reads)
# average for each patient
mr <- NULL
for (i in 1:length(unique(cellSums_dt$patient))){
  if (i==1){
    mr <- mean(cellSums_dt$reads[cellSums_dt$patient == unique(cellSums_dt$patient)[i]])
  } else{
    mr <- c(mr, mean(cellSums_dt$reads[cellSums_dt$patient == unique(cellSums_dt$patient)[i]]))
  }
}
names(mr) <- unique(cellSums_dt$patient)
write.csv(mr, file = paste0(data.loc, 'detectedReads_meanPerPatient.csv'))

# --------------- genes per cell - min. 1 reads --------------------

gene.matrix <- count.matrix
gene.matrix[count.matrix < 1] <- 0
gene.matrix[count.matrix >= 1] <- 1

cellSums <- colSums(gene.matrix)

cellSums_dt <- as.data.table(cellSums)
colnames(cellSums_dt) <- "reads"
rownames(cellSums_dt) <- colnames(count.matrix)
cellSums_dt$patient <- sub("\\..*$","",rownames(cellSums_dt))

pdf(sprintf(paste0(data.loc, "detectedGenes.pdf"), i),width=18,height=6,paper='special') 

print(ggplot(data = cellSums_dt, aes(x=patient, y=reads, fill=patient)) + 
        geom_violin(trim=T) +
        geom_jitter(shape=1, position=position_jitter(0.2), size=0.1, color="gray60") +
        geom_boxplot(width=0.1, size=1.2) +
        scale_fill_brewer(palette="Set3") +
        scale_color_brewer(palette="Set3") +
        theme(legend.position="none"))

dev.off()

# average across patients
mean(cellSums_dt$reads)
# average for each patient
mr <- NULL
for (i in 1:length(unique(cellSums_dt$patient))){
  if (i==1){
    mr <- mean(cellSums_dt$reads[cellSums_dt$patient == unique(cellSums_dt$patient)[i]])
  } else{
    mr <- c(mr, mean(cellSums_dt$reads[cellSums_dt$patient == unique(cellSums_dt$patient)[i]]))
  }
}
names(mr) <- unique(cellSums_dt$patient)
write.csv(mr, file = paste0(data.loc, "detectedGenes_meanPerPatient.csv")
