### Run scissor package

library(Scissor)
library(org.Hs.eg.db)

####################### ADAMS scRNA-seq reference dataset with TGFbeta1 timecourse bulk data ########################

#The package requires a single cell expression matrix. First lets use Adams with the normalised counts:
setwd("~/Documents/Chambers/COSBI")
adams <- LoadSeuratRds("GSE136831_fib.rds")
adams_sc_matrix <- as.matrix(adams@assays$RNA$data) #Take the normalised counts
adams_sc_dataset <- Seurat_preprocessing(adams_sc_matrix, verbose = F) #Pre-process in the way they want, add cell types later in metadata

#The package requires a bulk expression matrix. It seems this should be normalised counts.
#If you want to use the Plasmax timecourse:
setwd("~/Documents/Chambers/RNA-seq Plasmax timecourse")
TGF_bulk_dataset <- read.csv("TGFtimecourse-normalisedcounts-JW.csv")
setwd("~/Documents/Chambers/COSBI")
TGF_bulk_dataset <- as.matrix(TGF_bulk_dataset) #Make it into a matrix and add symbols
rownames(TGF_bulk_dataset) <- TGF_bulk_dataset[,1]
TGF_bulk_dataset <- TGF_bulk_dataset[,-1]
anno2 <- mapIds(org.Hs.eg.db,keys=rownames(TGF_bulk_dataset),
               column=c("SYMBOL"),
               keytype="ENTREZID")
rownames(TGF_bulk_dataset) <- anno2
colnames(TGF_bulk_dataset) <- paste(1:ncol(TGF_bulk_dataset))


#The package requires a phenotype list (conditions in the bulk)
#All conditions for reference:
phenotype <- matrix(
  c(rep(c(0, 1), each = 3, times = 4)),
  ncol = 1,
)
colnames(phenotype) <- c("Phenotype")

#Run scissor using a binomial as we have a binary variable (0 = Control, 1 = TGF)
tag <- c('Control', 'TGF-β1')
infos2 <- Scissor(TGF_bulk_dataset, adams_sc_dataset, phenotype, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Adams.RData")

#Visualise on UMAP. Red are scissor + (TGF-B1 associated) and blue are scissor - (control associated)
Scissor_select <- rep(0, ncol(adams_sc_dataset))
names(Scissor_select) <- colnames(adams_sc_dataset)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
adams_sc_dataset <- AddMetaData(adams_sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(adams_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

#Plot gene of interest (the below are the same irrespective of which bulk dataset you give)
FeaturePlot(adams_sc_dataset, features = c("COL1A1"))
FeaturePlot(adams_sc_dataset, features = c("ACTA2"))
FeaturePlot(adams_sc_dataset, features = c("CTHRC1"))
FeaturePlot(adams_sc_dataset, features = c("HAS1"))
FeaturePlot(adams_sc_dataset, features = c("PLIN2"))

#Add disease identity to scissor object
adams_sc_dataset@meta.data$Disease_Identity <- adams@meta.data$Disease_Identity
DimPlot(adams_sc_dataset, group.by = "Disease_Identity", reduction = "umap")

#Add cell labels to scissor object
adams_sc_dataset@meta.data$c_cell_type <- adams@meta.data$c_cell_type
DimPlot(adams_sc_dataset, group.by = "c_cell_type", reduction = "umap")

#Test reliability of your scissor run. n = number of permutations
#numbers <- length(infos2$Scissor_pos) + length(infos2$Scissor_neg)
#result <- reliability.test(X, Y, network, alpha = 0.2, family = "binomial", cell_num = numbers, n = 10, nfold = 10)

#Use FindMarkers function in Seurat to look at the differentially expressed genes between Scissor - and scissor + cells. Metadata is called "scissor"
Idents(adams_sc_dataset) <- adams_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
tgf.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0.25)
ctrl.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(2), min.pct = 0.25, logfc.threshold = 0.25)

write.csv(tgf.markers, file = "tgf_scissorpos_markers.csv")
write.csv(ctrl.markers, file = "ctrl_scissorneg_markers.csv")



#Edit the files if you only want to use one time-point
#24h is the last six samples
TGF_bulk_dataset_24h <- TGF_bulk_dataset[,19:24]
phenotype_24h <- matrix(
  c(rep(c(0, 1), each = 3, times = 1)),
  ncol = 1,
)
colnames(phenotype_24h) <- c("Phenotype")
tag <- c('Control', 'TGF-β1')
infos2 <- Scissor(TGF_bulk_dataset_24h, adams_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Adams_24h.RData")

#1.5h is the first six samples (same phenotype and tag)
TGF_bulk_dataset_1.5h <- TGF_bulk_dataset[,1:6]
infos2 <- Scissor(TGF_bulk_dataset_1.5h, adams_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Adams_1.5h.RData")

#6h is the second six samples
TGF_bulk_dataset_6h <- TGF_bulk_dataset[,7:12]
infos2 <- Scissor(TGF_bulk_dataset_6h, adams_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Adams_6h.RData")

#12h is the third six samples
TGF_bulk_dataset_12h <- TGF_bulk_dataset[,13:18]
infos2 <- Scissor(TGF_bulk_dataset_12h, adams_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Adams_12h.RData")



####################### HABERMANN scRNA-seq reference dataset with TGFbeta1 timecourse bulk data ########################
#If we want to use Habermann:
habermann <- LoadSeuratRds("GSE135893_fib.rds")
habermann_sc_matrix <- as.matrix(habermann@assays$RNA$data)
habermann_sc_dataset <- Seurat_preprocessing(habermann_sc_matrix, verbose = F)

#All conditions for reference:
phenotype <- matrix(
  c(rep(c(0, 1), each = 3, times = 4)),
  ncol = 1,
)
colnames(phenotype) <- c("Phenotype")

#Run scissor on Habermann
tag <- c('Control', 'TGF-β1')
infos2 <- Scissor(TGF_bulk_dataset, habermann_sc_dataset, phenotype, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Habermann.RData")

#Visualise on UMAP. Red are scissor + (TGF-B1 associated) and blue are scissor - (control associated)
Scissor_select <- rep(0, ncol(habermann_sc_dataset))
names(Scissor_select) <- colnames(habermann_sc_dataset)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
habermann_sc_dataset <- AddMetaData(habermann_sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(habermann_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

#Plot gene of interest (the below are the same irrespective of which bulk dataset you give)
FeaturePlot(habermann_sc_dataset, features = c("COL1A1"))
FeaturePlot(habermann_sc_dataset, features = c("ACTA2"))
FeaturePlot(habermann_sc_dataset, features = c("CTHRC1"))
FeaturePlot(habermann_sc_dataset, features = c("PLIN2"))
FeaturePlot(habermann_sc_dataset, features = c("HAS1"))

#Add disease identity to scissor object
habermann_sc_dataset@meta.data$Status <- habermann@meta.data$Status
DimPlot(habermann_sc_dataset, group.by = "Status", reduction = "umap")

#Add cell labels to scissor object
habermann_sc_dataset@meta.data$celltype <- habermann@meta.data$celltype
DimPlot(habermann_sc_dataset, group.by = "celltype", reduction = "umap")

#Choose the bulk time-point above as we did for Adams, then run these:
#24h
infos2 <- Scissor(TGF_bulk_dataset_24h, habermann_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Habermann_24h.RData")

#1.5h
infos2 <- Scissor(TGF_bulk_dataset_1.5h, habermann_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Habermann_1.5h.RData")

#6h
infos2 <- Scissor(TGF_bulk_dataset_6h, habermann_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Habermann_6h.RData")

#12h
infos2 <- Scissor(TGF_bulk_dataset_12h, habermann_sc_dataset, phenotype_24h, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_Habermann_12h.RData")



############### Bulk RNA-Seq from CRISPR-Cas9 gene edited pHLFs mapped to Habermann or Adams scRNA-seq #################

#If you want to use the CRISPR dataset, grab the single cell matrix of interest from above then:
setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022")
TGF_bulk_dataset <- read.csv("scissor_normalisedcounts_CRISPR.csv")
setwd("~/Documents/Chambers/COSBI")
TGF_bulk_dataset <- as.matrix(TGF_bulk_dataset)

rownames(TGF_bulk_dataset) <- TGF_bulk_dataset[,1]
TGF_bulk_dataset <- TGF_bulk_dataset[,-1]
anno2 <- mapIds(org.Hs.eg.db,keys=rownames(TGF_bulk_dataset),
                column=c("SYMBOL"),
                keytype="ENTREZID")
rownames(TGF_bulk_dataset) <- anno2
colnames(TGF_bulk_dataset) <- paste(1:ncol(TGF_bulk_dataset))

#Remove the 614 NAs in symbol column
#has_na <- is.na(TGF_bulk_dataset[, 2])
#TGF_bulk_dataset <- TGF_bulk_dataset[!has_na, ]


#Control cells only (first six samples)
TGF_bulk_dataset_crCTRL <- TGF_bulk_dataset[,1:6]
phenotype_crCTRL <- matrix(
  c(rep(c(0, 1), each = 3, times = 1)),
  ncol = 1,
)
colnames(phenotype_crCTRL) <- c("Phenotype")
tag <- c('Control', 'TGF-β1')

#Adams alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crCTRL, adams_sc_dataset, phenotype_crCTRL, tag = tag, 
                  alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crCTRL_Adams.RData")
#Adams alpha 0.2
infos4 <- Scissor(TGF_bulk_dataset_crCTRL, adams_sc_dataset, phenotype_crCTRL, tag = tag, 
                  alpha = 0.2, 
                  family = "binomial", Save_file = "Scissor_crCTRL_Adams_alpha0.2.RData")

#Habermann alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crCTRL, habermann_sc_dataset, phenotype_crCTRL, tag = tag, 
                  alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crCTRL_Habermann.RData")

#RAPTOR/mTORC1-depleted cells only (second six samples)
TGF_bulk_dataset_crRAPTOR <- TGF_bulk_dataset[,7:12]
phenotype_crRAPTOR <- matrix(
  c(rep(c(0, 1), each = 3, times = 1)),
  ncol = 1,
)
colnames(phenotype_crRAPTOR) <- c("Phenotype")

#Adams alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crRAPTOR, adams_sc_dataset, phenotype_crRAPTOR, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crRAPTOR_Adams.RData")
#Adams alpha 0.2
infos4 <- Scissor(TGF_bulk_dataset_crRAPTOR, adams_sc_dataset, phenotype_crRAPTOR, tag = tag, 
                  alpha = 0.2, 
                  family = "binomial", Save_file = "Scissor_crRAPTOR_Adams_alpha0.2.RData")

#Haberman alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crRAPTOR, habermann_sc_dataset, phenotype_crRAPTOR, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crRAPTOR_Habermann.RData")

#Adams alpha 0.01
infos4 <- Scissor(TGF_bulk_dataset_crRAPTOR, adams_sc_dataset, phenotype_crRAPTOR, tag = tag, 
                  alpha = 0.0005, 
                  family = "binomial", Save_file = "Scissor_crRAPTOR_Adams_alpha0.01.RData")

#RICTOR/mTORC2-depleted cells only (final six samples)
TGF_bulk_dataset_crRICTOR <- TGF_bulk_dataset[,13:18]
phenotype_crRICTOR <- matrix(
  c(rep(c(0, 1), each = 3, times = 1)),
  ncol = 1,
)
colnames(phenotype_crRICTOR) <- c("Phenotype")

#Adams alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crRICTOR, adams_sc_dataset, phenotype_crRICTOR, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crRICTOR_Adams.RData")
#Adams alpha 0.2
infos4 <- Scissor(TGF_bulk_dataset_crRICTOR, adams_sc_dataset, phenotype_crRICTOR, tag = tag, 
                  alpha = 0.2, 
                  family = "binomial", Save_file = "Scissor_crRICTOR_Adams_alpha0.2.RData")

#Habermann alpha 0.5
infos3 <- Scissor(TGF_bulk_dataset_crRICTOR, habermann_sc_dataset, phenotype_crRICTOR, tag = tag, alpha = 0.5, 
                  family = "binomial", Save_file = "Scissor_crRICTOR_Habermann.RData")


#Visualisation when using Adams as scRNA-seq reference:
#Visualise on UMAP. Red are scissor + (TGF-B1 associated) and blue are scissor - (control associated)
Scissor_select <- rep(0, ncol(adams_sc_dataset))
names(Scissor_select) <- colnames(adams_sc_dataset)
Scissor_select[infos3$Scissor_pos] <- 1
Scissor_select[infos3$Scissor_neg] <- 2
adams_sc_dataset <- AddMetaData(adams_sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(adams_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

#If you want to change the alpha (transparency of points)
DimPlot(adams_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = alpha(c('grey','indianred1','royalblue'), 0.3), pt.size = 1.2, order = c(2,1))

#For RAPTOR only because there are no Scissor+ cells
DimPlot(adams_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','royalblue','indianred1'), pt.size = 1.2, order = c(2,1))

#Plot gene of interest (the below are the same irrespective of which bulk dataset you give)
FeaturePlot(adams_sc_dataset, features = c("COL1A1"))
FeaturePlot(adams_sc_dataset, features = c("ACTA2"))
FeaturePlot(adams_sc_dataset, features = c("CTHRC1"))
FeaturePlot(adams_sc_dataset, features = c("TWIST1"))
FeaturePlot(adams_sc_dataset, features = c("POSTN"))

#Add disease identity to scissor object
adams_sc_dataset@meta.data$Disease_Identity <- adams@meta.data$Disease_Identity
DimPlot(adams_sc_dataset, group.by = "Disease_Identity", reduction = "umap")

#Add cell labels to scissor object
adams_sc_dataset@meta.data$c_cell_type <- adams@meta.data$c_cell_type
DimPlot(adams_sc_dataset, group.by = "c_cell_type", reduction = "umap")

#Use FindMarkers function in Seurat to look at the differentially expressed genes between Scissor - and scissor + cells. Metadata is called "scissor"
Idents(adams_sc_dataset) <- adams_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
tgf.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0.25)
ctrl.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(2), min.pct = 0.25, logfc.threshold = 0.25)

write.csv(tgf.markers, file = "crCTRL_TGF_scissorpos_markers.csv")
write.csv(ctrl.markers, file = "crCTRL_MC_scissorneg_markers.csv")

write.csv(tgf.markers, file = "crRAPTOR_TGF_scissorpos_markers.csv")
write.csv(ctrl.markers, file = "crRAPTOR_MC_scissorneg_markers.csv")

write.csv(tgf.markers, file = "crRICTOR_TGF_scissorpos_markers.csv")
write.csv(ctrl.markers, file = "crRICTOR_MC_scissorneg_markers.csv")


#Visualisation when using Habermann as scRNA-seq reference:
#Visualise on UMAP. Red are scissor + (TGF-B1 associated) and blue are scissor - (control associated)
Scissor_select <- rep(0, ncol(habermann_sc_dataset))
names(Scissor_select) <- colnames(habermann_sc_dataset)
Scissor_select[infos3$Scissor_pos] <- 1
Scissor_select[infos3$Scissor_neg] <- 2
habermann_sc_dataset <- AddMetaData(habermann_sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(habermann_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

#For RAPTOR only because there are no Scissor+ cells
DimPlot(habermann_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','royalblue','indianred1'), pt.size = 1.2, order = c(2,1))

#Plot gene of interest (the below are the same irrespective of which bulk dataset you give)
FeaturePlot(habermann_sc_dataset, features = c("COL1A1"))
FeaturePlot(habermann_sc_dataset, features = c("ACTA2"))
FeaturePlot(habermann_sc_dataset, features = c("CTHRC1"))
FeaturePlot(habermann_sc_dataset, features = c("TWIST1"))
FeaturePlot(habermann_sc_dataset, features = c("POSTN"))

#Add disease identity to scissor object
habermann_sc_dataset@meta.data$Disease_Identity <- habermann@meta.data$Disease_Identity
DimPlot(habermann_sc_dataset, group.by = "Disease_Identity", reduction = "umap")

#Add cell labels to scissor object
habermann_sc_dataset@meta.data$c_cell_type <- habermann@meta.data$c_cell_type
DimPlot(habermann_sc_dataset, group.by = "c_cell_type", reduction = "umap")

#Use FindMarkers function in Seurat to look at the differentially expressed genes between Scissor - and scissor + cells. Metadata is called "scissor"
Idents(habermann_sc_dataset) <- habermann_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
tgf.markers = FindMarkers(habermann_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0.25)
ctrl.markers = FindMarkers(habermann_sc_dataset, only.pos = FALSE, ident.1=c(2), min.pct = 0.25, logfc.threshold = 0.25)

write.csv(tgf.markers, file = "crCTRL_TGF_scissorpos_markers_habermann.csv")
write.csv(ctrl.markers, file = "crCTRL_MC_scissorneg_markers_habermann.csv")

write.csv(tgf.markers, file = "crRAPTOR_TGF_scissorpos_markers_habermann.csv")
write.csv(ctrl.markers, file = "crRAPTOR_MC_scissorneg_markers_habermann.csv")

write.csv(tgf.markers, file = "crRICTOR_TGF_scissorpos_markers_habermann.csv")
write.csv(ctrl.markers, file = "crRICTOR_MC_scissorneg_markers_habermann.csv")


##Tuning the parameters
#Try different alpha parameters
infos.alpha <- Scissor(TGF_bulk_dataset_crCTRL, adams_sc_dataset, phenotype_crCTRL, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03,
                       family = "binomial", Load_file = "Scissor_crCTRL_Adams.RData")
infos.alpha <- Scissor(TGF_bulk_dataset_crRAPTOR, adams_sc_dataset, phenotype_crRAPTOR, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03,
                       family = "binomial", Load_file = "Scissor_crRAPTOR_Adams.RData")
infos.alpha <- Scissor(TGF_bulk_dataset_crRICTOR, adams_sc_dataset, phenotype_crRICTOR, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03, #The alpha will stop searching when the total percentage of selected cells is less than 3%, usually set to 20% (0.2)
                       family = "binomial", Load_file = "Scissor_crRICTOR_Adams.RData")


#Test reliability of your scissor run. n = number of permutations. You have to first load the .RData object (in COSBI folder) to get the X and Y values and the network.
#Test for alpha 0.5
load("Scissor_crCTRL_Adams.RData")
load("Scissor_crRAPTOR_Adams.RData")
load("Scissor_crRICTOR_Adams.RData")
numbers <- length(infos3$Scissor_pos) + length(infos3$Scissor_neg)
result <- reliability.test(X, Y, network, alpha = 0.5, family = "binomial", cell_num = numbers, 
                           n = 2, nfold = 2#k-fold cross-validation parameters, n=number of permutations, nfold=number of splits in the data (2 means half data is compared to other half)
)

#Test for alpha 0.2
load("Scissor_crCTRL_Adams_alpha0.2.RData")
load("Scissor_crRAPTOR_Adams_alpha0.2.RData")
load("Scissor_crRICTOR_Adams_alpha0.2.RData")
numbers <- length(infos4$Scissor_pos) + length(infos4$Scissor_neg)
result <- reliability.test(X, Y, network, alpha = 0.2, family = "binomial", cell_num = numbers, 
                           n = 2, nfold = 2  #k-fold cross-validation parameters, n=number of permutations, nfold=number of splits in the data (2 means half data is compared to other half)
)



############### Bulk RNA-Seq from dual mTOR inhibited (AZD8055) or partial mTORC1 inhibited (rapamycin) pHLFs mapped to Habermann or Adams scRNA-seq #################
#Single cell matrix Adams
setwd("~/Documents/Chambers/COSBI")
adams <- LoadSeuratRds("GSE136831_fib.rds")
adams_sc_matrix <- as.matrix(adams@assays$RNA$data) #Take the normalised counts
adams_sc_dataset <- Seurat_preprocessing(adams_sc_matrix, verbose = F) #Pre-process in the way they want, add cell types later in metadata

#If you want to use the AZD/Rapa dataset
setwd("~/Documents/Chambers/AZD analysis YaleSummer2022")
TGF_bulk_dataset <- read.csv("diffexpr_results_CTRL_TGF_response.unfiltered.csv")
setwd("~/Documents/Chambers/COSBI")
TGF_bulk_dataset <- as.matrix(TGF_bulk_dataset)
rownames(TGF_bulk_dataset) <- TGF_bulk_dataset[,2]
TGF_bulk_dataset <- TGF_bulk_dataset[,10:33]
TGF_bulk_dataset.num <- matrix(as.numeric(TGF_bulk_dataset), ncol=ncol(TGF_bulk_dataset)) #Need to make it numeric
colnames(TGF_bulk_dataset.num) <- paste(1:ncol(TGF_bulk_dataset.num))
rownames(TGF_bulk_dataset.num) <- rownames(TGF_bulk_dataset)

#Prep phenotype
phenotype_DMSO <- matrix(
  c(rep(c(0, 1), each = 4, times = 1)),
  ncol = 1,
)
colnames(phenotype_DMSO) <- c("Phenotype")
tag <- c('Control', 'TGF-β1')

#DMSO Control only (first eight samples as we have n=4)
TGF_bulk_dataset_DMSO <- TGF_bulk_dataset.num[,1:8]
infos5 <- Scissor(TGF_bulk_dataset_DMSO, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                  alpha = 0.05, 
                  family = "binomial", Save_file = "Scissor_DMSO_Adams.RData")

#Rapamycin only (second eight samples)
TGF_bulk_dataset_Rapa <- TGF_bulk_dataset.num[,9:16]
infos5 <- Scissor(TGF_bulk_dataset_Rapa, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                  alpha = 0.05, 
                  family = "binomial", Save_file = "Scissor_Rapa_Adams.RData")

#AZD8055 only (final eight samples)
TGF_bulk_dataset_AZD <- TGF_bulk_dataset.num[,17:24]
infos5 <- Scissor(TGF_bulk_dataset_AZD, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                  alpha = 0.05, 
                  family = "binomial", Save_file = "Scissor_AZD_Adams.RData")


#Tuning the parameters
#Try different alpha parameters
infos.alpha <- Scissor(TGF_bulk_dataset_DMSO, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03,
                       family = "binomial", Load_file = "Scissor_DMSO_Adams.RData")
infos.alpha <- Scissor(TGF_bulk_dataset_Rapa, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03,
                       family = "binomial", Load_file = "Scissor_Rapa_Adams.RData")
infos.alpha <- Scissor(TGF_bulk_dataset_AZD, adams_sc_dataset, phenotype_DMSO, tag = tag, 
                       alpha = NULL, 
                       cutoff = 0.03, #The alpha will stop searching when the total percentage of selected cells is less than 3%, usually set to 20% (0.2)
                       family = "binomial", Load_file = "Scissor_AZD_Adams.RData")


#Test reliability of your scissor run. n = number of permutations. You have to first load the .RData object (in COSBI folder) to get the X and Y values and the network.
#Test for alpha 0.5
load("Scissor_crCTRL_Adams.RData")
load("Scissor_crRAPTOR_Adams.RData")
load("Scissor_crRICTOR_Adams.RData")
numbers <- length(infos3$Scissor_pos) + length(infos3$Scissor_neg)
result <- reliability.test(X, Y, network, alpha = 0.5, family = "binomial", cell_num = numbers, 
                           n = 2, nfold = 2#k-fold cross-validation parameters, n=number of permutations, nfold=number of splits in the data (2 means half data is compared to other half)
)


#Visualise on UMAP. Red are scissor + (TGF-B1 associated) and blue are scissor - (control associated)
Scissor_select <- rep(0, ncol(adams_sc_dataset))
names(Scissor_select) <- colnames(adams_sc_dataset)
Scissor_select[infos5$Scissor_pos] <- 1
Scissor_select[infos5$Scissor_neg] <- 2
adams_sc_dataset <- AddMetaData(adams_sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(adams_sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

#Use FindMarkers function in Seurat to look at the differentially expressed genes between Scissor - and scissor + cells. Metadata is called "scissor"
Idents(adams_sc_dataset) <- adams_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
tgf.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0.25)
ctrl.markers = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(2), min.pct = 0.25, logfc.threshold = 0.25)

write.csv(tgf.markers, file = "DMSO_TGF_scissorpos_markers_adams.csv")
write.csv(ctrl.markers, file = "DMSO_MC_scissorneg_markers_adams.csv")

write.csv(tgf.markers, file = "Rapa_TGF_scissorpos_markers_adams.csv")
write.csv(ctrl.markers, file = "Rapa_MC_scissorneg_markers_adams.csv")

write.csv(tgf.markers, file = "AZD_TGF_scissorpos_markers_adams.csv")
write.csv(ctrl.markers, file = "AZD_MC_scissorneg_markers_adams.csv")



##### Visualisation of CRISPR/Cas9 RNA-seq Scissor analysis for CTHRC1 pathogenic fibroblasts paper:
#Note the fold change cut off was changed from 0.25 above to 0.58 to be consistent with rest of paper

## Compare Scissor plus cells with all other cells
#Use FindMarkers function in Seurat to look at the differentially expressed genes between Scissor - and scissor + cells. Metadata is called "scissor"
Idents(adams_sc_dataset) <- adams_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
scissor.pos.markers.all = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0) #Scissor positive, get all genes for the volcano plot
write.csv(scissor.pos.markers.all, file = "crCTRL_TGF_scissorpos_markers_forvolcano.csv")

#Use FindMarkers function in Seurat to look at the differentially expressed genes of scissor + cells. Metadata is called "scissor"
Idents(adams_sc_dataset) <- adams_sc_dataset@meta.data$'scissor' #Change identity so it now uses scissor as the identity for each cell
scissor.pos = FindMarkers(adams_sc_dataset, only.pos = FALSE, ident.1=c(1), min.pct = 0.25, logfc.threshold = 0.58)
write.csv(scissor.pos, file = "crCTRL_TGF_scissorpos_markers_FC0.58.csv")


## Volcano plots of Scissor+ genes
#Load packages into R session
library(EnhancedVolcano)
pacman::p_load_gh("trinker/textshape")

#Load the Scissor pos genes (all including those with low or zero FC)
setwd("~/Documents/Chambers/COSBI/Volcanos-pathways-analysis")
scissor.pos.markers.all <- read.table("crCTRL_TGF_scissorpos_markers_forvolcano.txt", header=TRUE, comment.char="#")

#Make the gene symbols into row columns
scissor.pos.markers.all <- column_to_rownames(scissor.pos.markers.all, loc=1)


#Set publication theme
mytheme2 = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))


#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
  #Load data
  scissor.pos.markers.all,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  
  #Labelling
  lab = rownames(scissor.pos.markers.all),
  selectLab = c('ASPN','POSTN','CTHRC1','ADAM12','TGFBI','COL1A1'),
  cutoffLineType = 'dashed',
  legendPosition = 'right',
  drawConnectors = TRUE,
  arrowheads=F,
  labSize=3,
  
  title = '', 
  subtitle = '',
  #legendlabels = c('test', 'mTOR-dependent'),
  #legendPosition = 'right',
  caption = '',
  
  #Cutoffs
  pCutoff = 5*10e-03,
  FCcutoff = 0.58,
  
  #Size colour shape of points
  col = c('darkgrey','darkgrey','darkgrey','red'),
  colAlpha=1,
  #colCustom = keyvals,
  pointSize = 1.5
)



#Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-2.5, 2),ylim=c(0,20)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-2.5, 2, 0.5)) +
  mytheme2 +
  theme(legend.title = element_blank())




## Pathways analysis (dot plots)

#Load packages into R session
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#Load the Scissor pos genes (FC set above, should be 0.58)
setwd("~/Documents/Chambers/COSBI/Volcanos-pathways-analysis")
scissor.pos.fc0.58.q0.05 <- read.table("crCTRL_TGF_scissorpos_markers_FC0.58_padj0.05.txt", header=TRUE, comment.char="#")

de <- as.character(scissor.pos.fc0.58.q0.05[,1])
head(de)
entrez_ids <- mapIds(org.Hs.eg.db, keys = de, column = "ENTREZID", keytype = "SYMBOL")
head(entrez_ids)

#Pathway enrichment using ReactomePA
x <- enrichPathway(gene=entrez_ids, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE, minGSSize = 10, maxGSSize = 500)
head(x)

y <- enrichKEGG(gene=entrez_ids,
                organism = 'hsa',
                qvalueCutoff = 0.05,
                minGSSize = 10
)

#Save analysis file
write.csv(x, file="scissor-pos-reactome-pathways.csv")
write.csv(y, file="scissor-pos-kegg-pathways.csv")

#Visualise using bar plot
barplot(x, showCategory=20,
        label_format = function(x) stringr::str_wrap(x, width=60)) +
  scale_fill_gradient("Adjusted p-value", high = "#27214d", low = "red") +
  mytheme




## Violin plots of the key genes
library(ggplot2)
library(Seurat)

mytheme3 = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"), legend.position = 'none'
    ))

genes <- c("ASPN", "POSTN", "CTHRC1", "ADAM12", "TGFBI", "COL1A1")

VlnPlot(adams_sc_dataset, features = genes,
        cols=c("grey", "blue","red"),
        pt.size = 0) &
  scale_x_discrete(labels = c("0" = "Background",
                              "1" = "Scissor+",
                              "2" = "Scissor-")) &
  mytheme3





