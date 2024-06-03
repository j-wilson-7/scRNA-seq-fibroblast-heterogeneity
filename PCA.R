# Clustering, UMAP and DGE of publicly available IPF scRNA-seq dataset (fibroblasts only) from Tsukui et al.

#Install packages
BiocManager::install("scDblFinder")

#Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
#library(SingleR)
library(ggplot2)
library(scater)
library(celldex)
library(DoubletFinder)
library(scDblFinder)

## Prepare the data
#Load the filtered object
data <- readRDS("/home/jowilson/cosbi/tsukui/seurat_GSE132771_f.rds")

#Create metadata containing the sample names
sample_names <- substr(colnames(data),1,15)
data$sample_name <- sample_names

#Create metadata containing disease or non-disease info
diseases <- substr(colnames(data),12, 14)
data$disease <- diseases

## LogNormalise the data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection) - genes that differ lots between cells
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with top 10 labelled
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
pdf(file = 'variable-features.pdf', width = 12)
print(plot2)
dev.off()



## Perform PCA
#Scale the data to ensure highly expressed genes don't dominate the analysis, the default is to do ScaleData only on those top 2000 genes we identified earlier:
data <- ScaleData(data)

#The previously defined 2000 top variable features will be used as input:
data <- RunPCA(data, features = VariableFeatures(object = data))

#Look at the top genes for each PC
PC1PC2 <- VizDimLoadings(data, dims = 1:2, reduction = "pca")
PC1PC2
pdf(file = 'PC1PC2.pdf', width = 12)
print(PC1PC2)
dev.off()

#Look at a PCA plot 
pca_plot <- DimPlot(data, group.by = "sample_name", reduction = "pca") + 
  ggtitle("PCA plot")
pca_plot
pdf(file = 'PCAplot.pdf', width = 12)
print(pca_plot)
dev.off()

#Look at normal vs disease PCA plot
pca_NMLvIPF <- DimPlot(data, group.by = "disease", reduction = "pca") + 
  ggtitle("PCA plot")
pdf(file = 'PCA-NMLvIPF.pdf', width = 12)
print(pca_NMLvIPF)
dev.off()



## Determining number of PCs
#Calculate percent variation associated with each PC
  pct = data[["pca"]]@stdev/sum(data[["pca"]]@stdev)*100

#calculate cumulative percents for each PC
cumu = cumsum(pct)

#identify PC threshold where cumulative percentage > 90, but less than 5 percent of standard deviation
co1 = which(cumu > 90 & pct < 5)[1]
co1

#determine the percent difference between one PC and next
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1

#last point where difference in pct > 0.1
co2
pcs=min(co1,co2)

#Other ways to check number of clusters: JackStrawplot and Elbowplot
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:17)
ElbowPlot(data)
#################


## Clustering KNN K-nearest neighbour for pre doublet detection
#Number of PCs used in paper is 19 and resolution was 0.3
data <- FindNeighbors(data, dims = 1:19)
#The resolution is usually between 0.4-1.2 for datasets with around 3K cells
data <- FindClusters(data, resolution = seq(0.1,0.4, by = 0.1)) #Test a few different resolutions

#For post doublet detection use 21 PCs
data <- FindNeighbors(data, dims = 1:21)
data <- FindClusters(data, resolution = seq(0.1,0.4, by = 0.1))


## UMAP
#data <- RunUMAP(data, dims = 1:17) #For the pre-doublet detected dataset
data <- RunUMAP(data, dims = 1:21) #For the post-doublet detected dataset
DimPlot(data, reduction = "umap", group.by = "sample_name")
DimPlot(data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
#Save UMAP object so you can get it back later
saveRDS(data, file = "Amy-UMAP-object.rds")

Idents(data) <- "RNA_snn_res.0.1" #Changes it to use resolution 1


## Manual annotation:
# find markers for every cluster compared to all remaining cells, report only the positive ones
res1.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Report the top 2 for each cluster
res1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Save as an excel to look at the top markers for each cluster
write.csv(res1.markers, file = "markers_results.csv")
write.csv(res1.markers, file = "markers_results_postdoublet.csv")

# alternatively can find markers for a specific cluster relative to all the other clusters
#c5.markers = FindMarkers(data, only.pos = TRUE, ident.1=c(5), min.pct = 0.25, logfc.threshold = 0.25)
#c4.markers = FindMarkers(data, only.pos = TRUE, ident.1=c(4), min.pct = 0.25, logfc.threshold = 0.25)

#VinPlot() and FeaturePlot to look at gene expression
VlnPlot(data, features = c("COL1A1", "PDGFRB"))

#Mesenchymal cells:
FeaturePlot(data, features = c("PDGFRB"), label=TRUE)
FeaturePlot(data, features = c("PDGFRA"), label=TRUE)
FeaturePlot(data, features = c("MYH11"), label=TRUE)
FeaturePlot(data, features = c("LUM"), label=TRUE)

#Epithelial cells: (markers taken from Habermann et al.)
FeaturePlot(data, features = c("EPCAM"), label=TRUE)

#Immune cells:
FeaturePlot(data, features = c("PTPRC"), label=TRUE)
FeaturePlot(data, features = c("CD3E"), label=TRUE)
FeaturePlot(data, features = c("KLRB1"), label=TRUE)
FeaturePlot(data, features = c("MARCO"), label=TRUE)

#Endothelial cells:
FeaturePlot(data, features = c("PECAM1"), label=TRUE)




## Take only cluster 5 (fibroblasts)
cells_to_recluster <- which(data$RNA_snn_res.0.1 == "5")
cluster.5 <- data[, cells_to_recluster]

#Re-run seurat pipeline on it
cluster.5 <- NormalizeData(cluster.5)
cluster.5 <- ScaleData(cluster.5)
cluster.5 <- FindVariableFeatures(cluster.5, selection.method = "vst", nfeatures = 2000)
cluster.5 <- RunPCA(cluster.5, features = VariableFeatures(object = cluster.5))

#Define the dims
#Calculate percent variation associated with each PC
pct = cluster.5[["pca"]]@stdev/sum(cluster.5[["pca"]]@stdev)*100
#calculate cumulative percents for each PC
cumu = cumsum(pct)
#identify PC threshold where cumulative percentage > 90, but less than 5 percent of standard deviation
co1 = which(cumu > 90 & pct < 5)[1]
#determine the percent difference between one PC and next
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1
#last point where difference in pct > 0.1
pcs=min(co1,co2)
pcs
ElbowPlot(cluster.5)

#Plot new UMAP
cluster.5 <- RunUMAP(cluster.5, dims = 1:13)
cluster.5 <- FindNeighbors(cluster.5, dims = 1:13)
cluster.5 <- FindClusters(cluster.5, resolution = seq(0.1,0.4, by = 0.1)) #Test a few different resolutions
DimPlot(cluster.5, reduction = "umap", group.by = "RNA_snn_res.0.2", label = TRUE)
Idents(cluster.5) <- "RNA_snn_res.0.2" #Changes it to use resolution 0.2

#For post doublet detection:
DimPlot(cluster.5, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
Idents(cluster.5) <- "RNA_snn_res.0.1" #Changes it to use resolution 0.3


#Manual annotation
fibrocluster.markers <- FindAllMarkers(cluster.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save as an excel to look at the top markers for each cluster
write.csv(fibrocluster.markers, file = "fibrocluster_markers_postdoublet_results.csv")

#Check key genes
FeaturePlot(cluster.5, features = c("COL1A1"), label=TRUE)
FeaturePlot(cluster.5, features = c("ACTA2"), label=TRUE)
FeaturePlot(cluster.5, features = c("CTHRC1"), label=TRUE)
FeaturePlot(cluster.5, features = c("HAS1"), label=TRUE)
FeaturePlot(cluster.5, features = c("PLIN2"), label=TRUE)
FeaturePlot(cluster.5, features = c("LUM"), label=TRUE)
FeaturePlot(cluster.5, features = c("TWIST1"), label=TRUE)
FeaturePlot(cluster.5, features = c("MYLK"), label=TRUE)

#UMAP by disease
fibro_NMLvIPF <- DimPlot(cluster.5, group.by = "disease", label = TRUE, reduction = "UMAP")
fibro_NMLvIPF

#Check key genes for the post doublet corrected (there are two UMAPS, UMAP for pre doublet correction and umap for post)
FeaturePlot(cluster.5, features = c("COL1A1"), reduction = "umap")
FeaturePlot(cluster.5, features = c("ACTA2"), reduction = "umap")
FeaturePlot(cluster.5, features = c("CTHRC1"), reduction = "umap")
FeaturePlot(cluster.5, features = c("HAS1"), reduction = "umap")
FeaturePlot(cluster.5, features = c("PLIN2"), reduction = "umap")
FeaturePlot(cluster.5, features = c("LUM"), reduction = "umap")
FeaturePlot(cluster.5, features = c("TWIST1"), reduction = "umap")
FeaturePlot(cluster.5, features = c("MYLK"), reduction = "umap")





## Remove doublets from the original large data object with all clusters

#Make single cell experiment obj
sce <- as.SingleCellExperiment(data)

print(sce)
print(colnames(colData(sce)))

set.seed(123)

sce <- scDblFinder(sce, samples="sample_name")
table(sce$scDblFinder.class)

#Take only the singlets
singlets <- sce$scDblFinder.class == "singlet"
sce_singlets <- sce[, singlets]

#Make it back into a Seurat object and re-run everything above
data <- as.Seurat(sce_singlets)



## Remove doublets from the subset fibroblast populations (cluster 5 only)

#Make single cell experiment obj
sce5 <- as.SingleCellExperiment(cluster.5)
sce5 <- scDblFinder(sce5, samples="sample_name")
table(sce5$scDblFinder.class)

#Take only the singlets
singlets5 <- sce5$scDblFinder.class == "Singlet"
sce5_singlets <- sce5[, singlets5]
print(sce5_singlets)

#Make it back into a Seurat object and re-run everything above
cluster.5.sing <- as.Seurat(sce5_singlets)





## Redoing the analysis with integration. (Note this is for v4 Seurat, integration was updated in v5)

#Do the doublet removal first (from above)

#Then, split the seurat object based on IPF vs disease
split.data <- SplitObject(data, split.by = "disease")

#Normalise and identify variable features for each dataset independently
split.data <- lapply(X=split.data, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Select the features that are repeatedly variable across datasets (will be used as anchors during the integration)
features <- SelectIntegrationFeatures(object.list = split.data)
anchors <- FindIntegrationAnchors(object.list = split.data, anchor.features = features) #This step takes a while
int.data <- IntegrateData(anchorset = anchors)

#Then run the workflow as usual
int.data <- ScaleData(int.data)
int.data <- RunPCA(int.data, features=VariableFeatures(object = int.data))

pct = int.data[["pca"]]@stdev/sum(int.data[["pca"]]@stdev)*100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1
pcs=min(co1,co2)
pcs #number of PCs to use

int.data <- FindNeighbors(int.data, dims = 1:20)
int.data <- FindClusters(int.data, resolution = seq(0.1,0.4, by = 0.1)) #test a few different resolutions
int.data <- RunUMAP(int.data, dims = 1:20)
DimPlot(int.data, reduction = "umap", group.by = "sample_name")
DimPlot(int.data, reduction = "umap", group.by = "integrated_snn_res.0.1", label = TRUE) #check which resolution is best for your umap
Idents(int.data) <- "integrated_snn_res.0.1" #Changes it to use resolution 1
int.markers <- FindAllMarkers(int.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#UMAP by disease
DimPlot(int.data, group.by = "disease", reduction = "umap")

#Save the doublet removed, integrated object
saveRDS(int.data, file = "GSE132771_filtered_doublet_integrated.rds")

#Look at gene expression to define fibroblast clusters
VlnPlot(int.data, features = c("COL1A1", "PDGFRB"))
FeaturePlot(int.data, features = c("LUM"))

#Take the fibroblast cluster
cells_to_recluster <- which(int.data$integrated_snn_res.0.1 == "6")
cluster.6 <- data[, cells_to_recluster]

#Re-run seurat pipeline on it
cluster.6 <- NormalizeData(cluster.6)
cluster.6 <- ScaleData(cluster.6)
cluster.6 <- FindVariableFeatures(cluster.6, selection.method = "vst", nfeatures = 2000)
cluster.6 <- RunPCA(cluster.6, features = VariableFeatures(object = cluster.6))

#Define the dims
#Calculate percent variation associated with each PC
pct = cluster.6[["pca"]]@stdev/sum(cluster.6[["pca"]]@stdev)*100
#calculate cumulative percents for each PC
cumu = cumsum(pct)
#identify PC threshold where cumulative percentage > 90, but less than 5 percent of standard deviation
co1 = which(cumu > 90 & pct < 5)[1]
#determine the percent difference between one PC and next
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1
#last point where difference in pct > 0.1
pcs=min(co1,co2)
pcs

#Plot new UMAP
cluster.6 <- RunUMAP(cluster.6, dims = 1:13)
cluster.6 <- FindNeighbors(cluster.6, dims = 1:13)
cluster.6 <- FindClusters(cluster.6, resolution = seq(0.1,0.4, by = 0.1)) #Test a few different resolutions
DimPlot(cluster.6, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
Idents(cluster.6) <- "RNA_snn_res.0.1"

cluster.6 <- FindClusters(cluster.6, resolution = 0.05) #Try res 0.5
DimPlot(cluster.6, reduction = "umap", group.by = "RNA_snn_res.0.05", label = TRUE)

#Manual annotation
fibrocluster.markers <- FindAllMarkers(cluster.6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save as an excel to look at the top markers for each cluster
write.csv(fibrocluster.markers, file = "fibrocluster_markers_integrated_cluster6_results.csv")

#Look at gene expression to define the different fibroblasts
FeaturePlot(cluster.6, features = c("COL1A1"))
FeaturePlot(cluster.6, features = c("ACTA2"))
FeaturePlot(cluster.6, features = c("CTHRC1"))
FeaturePlot(cluster.6, features = c("HAS1"))
FeaturePlot(cluster.6, features = c("PLIN2"),)
FeaturePlot(cluster.6, features = c("LUM"))
FeaturePlot(cluster.6, features = c("TWIST1"))
FeaturePlot(cluster.6, features = c("MYLK"))

FeaturePlot(cluster.6, features = c("TEAD4"))

#UMAP by disease
DimPlot(cluster.6, group.by = "disease", reduction = "umap")
fibrocluster.IPFmarkers <- FindMarkers(cluster.6, group.by = "disease", ident.1 = "IPF", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(fibrocluster.IPFmarkers, file = "fibrocluster_IPF_markers_integrated_cluster6_results.csv")


###### TAKING A LOOK AT THE HABERMANN DATA (GSE135893 fibroblasts only)

#How many clusters?
pct = GSE135893_fib[["pca"]]@stdev/sum(GSE135893_fib[["pca"]]@stdev)*100
#calculate cumulative percents for each PC
cumu = cumsum(pct)
#identify PC threshold where cumulative percentage > 90, but less than 5 percent of standard deviation
co1 = which(cumu > 90 & pct < 5)[1]
#determine the percent difference between one PC and next
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1
#last point where difference in pct > 0.1
pcs=min(co1,co2)
pcs

#Plot new UMAP fibroblasts only
GSE135893_fib <- RunUMAP(GSE135893_fib, dims = 1:16)
GSE135893_fib <- FindNeighbors(GSE135893_fib, dims = 1:16)
GSE135893_fib <- FindClusters(GSE135893_fib, resolution = seq(0.1,0.4, by = 0.1)) #Test a few different resolutions
DimPlot(GSE135893_fib, reduction = "umap", group.by = "SCT_snn_res.0.1", label = TRUE)
Idents(GSE135893_fib) <- "SCT_snn_res.0.1"

#Gene expression
habermann.fib.markers <- FindAllMarkers(GSE135893_fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save as an excel to look at the top markers for each cluster
write.csv(habermann.fib.markers, file = "fibrocluster_markers_habermann_results.csv")

#Look at gene expression to define the different fibroblasts
DimPlot(GSE135893_fib, group.by = "Status", reduction = "umap")
FeaturePlot(GSE135893_fib, features = c("COL1A1"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("ACTA2"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("CTHRC1"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("HAS1"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("PLIN2"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("LUM"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("TWIST1"), label=TRUE)
FeaturePlot(GSE135893_fib, features = c("MYLK"), label=TRUE)

FeaturePlot(cluster.6, features = c("ACTA2"))
FeaturePlot(cluster.6, features = c("CTHRC1"))
FeaturePlot(cluster.6, features = c("HAS1"))
FeaturePlot(cluster.6, features = c("PLIN2"))
FeaturePlot(GSE135893_fib, features = c("LUM"), label=TRUE)
FeaturePlot(cluster.6, features = c("TWIST1"))
FeaturePlot(cluster.6, features = c("MYLK"))



#Analysing the Adams dataset (GSE136831 fibroblasts only)
#How many clusters? (12)
pct = GSE136831_fib[["pca"]]@stdev/sum(GSE136831_fib[["pca"]]@stdev)*100
#calculate cumulative percents for each PC
cumu = cumsum(pct)
#identify PC threshold where cumulative percentage > 90, but less than 5 percent of standard deviation
co1 = which(cumu > 90 & pct < 5)[1]
#determine the percent difference between one PC and next
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])>0.1), decreasing = T)[1]+1
#last point where difference in pct > 0.1
pcs=min(co1,co2)
pcs

#Plot new UMAP fibroblasts only
GSE136831_fib <- RunUMAP(GSE136831_fib, dims = 1:12)
GSE136831_fib <- FindNeighbors(GSE136831_fib, dims = 1:12)
GSE136831_fib <- FindClusters(GSE136831_fib, resolution = seq(0.1,0.4, by = 0.1)) #Test a few different resolutions
DimPlot(GSE136831_fib, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
Idents(GSE136831_fib) <- "RNA_snn_res.0.1"

#Gene expression
adams.fib.markers <- FindAllMarkers(GSE136831_fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save as an excel to look at the top markers for each cluster
write.csv(adams.fib.markers, file = "fibrocluster_markers_adams_results.csv")

#Look at gene expression to define the different fibroblasts
DimPlot(GSE136831_fib, group.by = "Disease_Identity", reduction = "umap")
FeaturePlot(GSE136831_fib, features = c("COL1A1"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("ACTA2"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("CTHRC1"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("HAS1"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("PLIN2"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("LUM"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("TWIST1"), label=TRUE)
FeaturePlot(GSE136831_fib, features = c("MYLK"), label=TRUE)










