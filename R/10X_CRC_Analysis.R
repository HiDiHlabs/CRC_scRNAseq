# Performed on the following samples
DefaultAssay(object = o1) <- "RNA"
DefaultAssay(object = x1) <- "RNA"
DefaultAssay(object = subset_tumor_epithelial_noHighMT) <- "RNA"

### AddModuleScore
# load genelists
genelist_NNMF <- read_delim("genelist_NNMF.csv", 
                            ";", quote = "'", escape_double = FALSE, 
                            trim_ws = TRUE)

list_oxphos8 <- list (c(genelist_NNMF$OXPHOS_8))
list_glycolysis2 <- list (c(genelist_NNMF$`Hypoxia/Glycolysis_2`))

# Calculating the geneset expression
# repeat for "nnmf_glycolysis21"
# repeat for different samples
subset_tumor_epithelial_noHighMT <- AddModuleScore(object = subset_tumor_epithelial_noHighMT, 
                                                   features = list_oxphos8, ctrl = 5, 
                                                   name = 'nnmf_oxphos8')

# Plotting (NNMF)
VlnPlot(object = subset_tumor_epithelial_noHighMT, features = "nnmf_oxphos81", pt.size = 0) +
  stat_summary(fun.y = median, geom='point', size = 3, colour = "black") 

# Same violin but with t-test
nnmf_oxphos81_df <- as.data.frame(subset_tumor_epithelial_noHighMT@meta.data$nnmf_oxphos81)
clusters_df <- as.data.frame(subset_tumor_epithelial_noHighMT@active.ident)
df_final <- cbind(nnmf_oxphos81_df, clusters_df$`subset_tumor_epithelial_noHighMT@active.ident`)
colnames(df_final) <- c("score", "cluster")

compare_means(score ~ cluster, data = df_final, method = "t.test")

my_comparisons <- list( c("TA", "Stem"), c("TA", "Tdiff"), c("Stem", "Tdiff") )

ggplot(df_final, aes(x = cluster, y = score, fill = cluster)) + 
  geom_violin(trim = F)  +
  stat_summary(fun.y = median, geom='point', size = 3, colour = "black") +
  stat_compare_means(comparisons = my_comparisons)

### Plotting Marker Genes
FeaturePlot(object = subset_tumor_epithelial_noHighMT,
            features = c("PROX1", "ASCL2",
                         "CDC20", "MCM6",
                         "FABP1", "KRT20"),
            cols = c("lightgrey", "red"), 
            pt.size = 1,
            reduction = "tsne",
            order = TRUE,
            ncol = 2, label = FALSE)

FeaturePlot(object = o1,
            features = c("PROX1", "SMOC2",
                         "CDC20", "MYC",
                         "CA9", "TFF3"),
            cols = c("lightgrey", "red"), 
            pt.size = 0.25,
            reduction = "tsne",
            order = TRUE,
            ncol = 2, label = FALSE)

FeaturePlot(object = x1,
            features = c("ASCL2", "RNF43",
                         "CDC20", "PA2G4",
                         "CA9", "ALDOA",
                         "REG4", "FRZB"),
            cols = c("lightgrey", "red"), 
            pt.size = 0.25,
            reduction = "tsne",
            order = TRUE,
            ncol = 2, label = FALSE)
