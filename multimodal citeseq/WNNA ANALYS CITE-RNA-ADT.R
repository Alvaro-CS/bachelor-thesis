library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

#We first perform pre-processing and dimensional reduction on both assays independently.
#We use standard normalization, but you can also use SCTransform or any alternative method.

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

#For each cell, we calculate its closest neighbors in the dataset based on a weighted combination of RNA and protein similarities.
#The cell-specific modality weights and multimodal neighbors are calculated in a single function, which takes ~2 minutes to run on this dataset.
#We specify the dimensionality of each modality (similar to specifying the number of PCs to include in scRNA-seq clustering),
#but you can vary these settings to see that small changes have minimal effect on the overall results.


# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight

bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)


#We can now use these results for downstream analysis, such as visualization and clustering.
#For example, we can create a UMAP visualization of the data based on a weighted combination of RNA and protein data.
#We can also perform graph-based clustering and visualize these results on the UMAP, alongside a set of cell annotations.
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2


#We can also compute UMAP visualization based on only the RNA and protein data and compare.
#We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states 
#(the ADT panel contains markers for differentiated cells), while the converse is true of T cell states (where the ADT analysis outperforms RNA).

bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

p3 <- DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

#We can visualize the expression of canonical marker genes and proteins on the multimodal UMAP, which can assist in verifying the provided annotations:

p5 <- FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(bm, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p5 / p6


#Finally, we can visualize the modality weights that were learned for each cell. Each of the populations with the highest RNA weights represent progenitor cells,
#while the populations with the highest protein weights represent T cells. This is in line with our biological expectations,
#as the antibody panel does not contain markers that can distinguish between different progenitor populations.
VlnPlot(bm, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend()
saveRDS(bm, file = "bm_final.rds")
