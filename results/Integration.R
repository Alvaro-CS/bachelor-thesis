
library(dplyr)
library(Seurat)
library(SeuratData)
library(SingleR)
pesa <- readRDS("PESA_S3000.rds")

#We want to remove neutrophils
DimPlot(pesa, reduction = "tSNE.100")
DimPlot(pesa, reduction = "umap")
# find markers for every cluster compared to all remaining cells, report only the positive ones
pesa.markers <- FindAllMarkers(pesa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pesa.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect).
#cluster1.markers <- FindMarkers(pesa, ident.1 = "C0", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
FeaturePlot(pesa, features = c("BCL6","MME","AQP9","CSF2RB","FGR","FPR1","DNAJB1","OSM","CD63","CD68","CD14","CD33", "CD44"))
DimPlot(pesa, reduction = "umap")

#Cell annotation
library(celldex) #Provides reference to several reference datasets from bulk or microarray data mostly
hpca.se <- HumanPrimaryCellAtlasData() #Specific reference
hpca.se #matrix of log-expression values with sample-level labels.

pred.pesa <- SingleR(test = pesa@assays$RNA@counts, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)

#Stats of predictions
table(pred.pesa$labels)
table(Label=pred.pesa$labels,Lost=is.na(pred.pesa$pruned.labels))
plotScoreHeatmap(pred.pesa)

pesa$celltypes <- pred.pesa$labels
##########
UMAPPlot(pesa)
p1=DimPlot(pesa, reduction = "umap",group.by = 'Cell_Type_Experimental') #automatic
p2=DimPlot(pesa, reduction = "umap",group.by = 'celltypes') #from singleR
p1+p2
##########
#Remove the neutrophiles
pesaNN<-subset(pesa, subset=celltypes =="Neutrophils", invert = TRUE)
p1=DimPlot(pesa, reduction = "umap",group.by = 'celltypes') #old
p2=DimPlot(pesaNN, reduction = "umap",group.by = 'celltypes') #No neutrophils
p1+p2


InstallData("bmciteD")
citeD <- LoadData(ds = "bmciteD")
DefaultAssay(citeD) <- 'RNA'

#for visualizing citeD
citeD <- ScaleData(citeD, verbose = FALSE)
citeD <- RunPCA(citeD, npcs = 30, verbose = FALSE)
citeD <- RunUMAP(citeD, reduction = "pca", dims = 1:30)
citeD <- FindNeighbors(citeD, reduction = "pca", dims = 1:30)
citeD <- FindClusters(citeD, resolution = 0.5)








######           TRANSFER         #######
DefaultAssay(pesaNN) <- 'RNA'
blood.anchors <- FindTransferAnchors(reference = citeD, query = pesaNN)
predictions <- TransferData(anchorset = blood.anchors, refdata = citeD$celltype.l2)

pesaNN <- AddMetaData(pesaNN, metadata = predictions)

##CREATE BLOOD LIST
blood.list <-list(pesaNN,citeD)
# normalize and identify variable features for each dataset independently
blood.list <- lapply (X=blood.list, FUN= function(x){
  x <- NormalizeData(x)
  x<- FindVariableFeatures(x,selection.method = "vst",nfeatures=2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = blood.list)

#We then identify anchors using the FindIntegrationAnchors() function, which takes 
#a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
#TAKES TIME FINDING ALL PAIRWISE ANCHORS. Found 12791 anchors, retained 5580 anchors, 9m 32 s
blood.anchors <- FindIntegrationAnchors(object.list = blood.list, anchor.features = features)

# this command creates an 'integrated' data assay
blood.combined <- IntegrateData(anchorset = blood.anchors)

########## INTEGRATE ANALYSIS
DefaultAssay(blood.combined) <- "integrated"
blood.combined <- ScaleData(blood.combined, verbose = FALSE)
blood.combined <- RunPCA(blood.combined, npcs = 30, verbose = FALSE)
blood.combined <- RunUMAP(blood.combined, reduction = "pca", dims = 1:30)
blood.combined <- FindNeighbors(blood.combined, reduction = "pca", dims = 1:30)
blood.combined <- FindClusters(blood.combined, resolution = 0.5)

blood.combined$celltype.final<-blood.combined$celltype.l2
nacells <- is.na(blood.combined$celltype.final)
blood.combined$celltype.final[nacells] <- blood.combined$predicted.id[which(nacells)]


######PLOTS#####

pciteD<-DimPlot(citeD, reduction = "umap", group.by = "celltype.l2")
pf <- DimPlot(blood.combined, reduction = "umap", group.by = "celltype.final")

pciteD+pf

#Observe contribution from each experiment
df <- Embeddings(blood.combined,reduction = "umap")
df <- cbind(df,blood.combined@meta.data)
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=celltype.final)) +
  geom_point() +
  facet_wrap("orig.ident")



#########
#Now, we export cell names with their celular types and counts into a table for simulating simple bulks from count matrix.
cell.metadata<-blood.combined@assays$RNA@counts
cell.metadata<-matrix(cell.metadata)
saveRDS(cell.metadata, file = "count_matrix.rds")
write.table(cell.metadata, file= "count_matrix.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = T)


cell.types<-blood.combined$celltype.final
saveRDS(cell.types, file = "cell_types.rds")
write.table(cell.types, file= "cell_types.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = F)

#########
