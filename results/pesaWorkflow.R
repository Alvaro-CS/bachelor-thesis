seurat3000 <- readRDS("PESA_S3000.rds")


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat3000[["percent.mt"]] <- PercentageFeatureSet(seurat3000, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(seurat3000@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(seurat3000, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seurat3000, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat3000, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Here we remove unwanted cells from the dataset
#seurat3000 <- subset(seurat3000, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#NORMALIZING THE DATA
#As default, we employ a global-scaling norm. method "LogNormalize"
seurat3000 <- NormalizeData(seurat3000, normalization.method = "LogNormalize", scale.factor = 10000)

#FEATURE SELECTION
seurat3000 <- FindVariableFeatures(seurat3000, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat3000), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat3000)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2 

#SCALING THE DATA
#Pre-processing step before dimensional reduction techniques like PCA
all.genes <- rownames(seurat3000)
seurat3000 <- ScaleData(seurat3000, features = all.genes)

#PERFORM LINEAR DIMENSIONAL REDUCTION
seurat3000 <- RunPCA(seurat3000, features = VariableFeatures(object = seurat3000))

# Examine and visualize PCA results a few different ways
print(seurat3000[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat3000, dims = 1:2, reduction = "pca")
DimPlot(seurat3000, reduction = "pca")
DimHeatmap(seurat3000, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat3000, dims = 1:15, cells = 500, balanced = TRUE)

#DETERMINE THE DIMENIONALY OF THE DATASET
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat3000 <- JackStraw(seurat3000, num.replicate = 100)
seurat3000 <- ScoreJackStraw(seurat3000, dims = 1:20)

JackStrawPlot(seurat3000, dims = 1:15)
ElbowPlot(seurat3000)


#CLUSTER THE CELLS
seurat3000 <- FindNeighbors(seurat3000, dims = 1:10)
seurat3000 <- FindClusters(seurat3000, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat3000), 5)

#RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/TSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seurat3000 <- RunUMAP(seurat3000, dims = 1:10)
seurat3000 <- RunTSNE(seurat3000, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat3000, reduction = "umap",label='FALSE')

saveRDS(seurat3000, file = "output/seurat3000_tutorial.rds")

#FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(seurat3000, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat3000, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat3000.markers <- FindAllMarkers(seurat3000, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat3000.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster1.markers <- FindMarkers(seurat3000, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Tools for visualizing marker expression
VlnPlot(seurat3000, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(seurat3000, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(seurat3000, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20? markers (or all markers if less than 20) for each cluster.
top10 <- seurat3000.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat3000, features = top10$gene) 

#ASSIGNING CELL TYPE IDENTITY TO CLUSTERS
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat3000)
seurat3000 <- RenameIdents(seurat3000, new.cluster.ids)
DimPlot(seurat3000, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(seurat3000, file = "output/seurat30003k_final.rds")
