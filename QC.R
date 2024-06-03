#Quality control of publicly available IPF scRNA-seq datasets (fibroblasts only) from Tsukui et al.
#In the paper they exclude cells which are five median absolute deviation distant
#from the median value of library size, number of detected genes or mitochondrial gene proportion.

#Define libraries to use
library(Seurat)
library(cowplot)

data <- readRDS("/home/jowilson/cosbi/tsukui/seurat_GSE132771.rds")

#Check the percentage of reads that map to mitochondrial genome, low-quality cells have mt contamination
percent.mt <- PercentageFeatureSet(data, pattern = "^MT-") 
data$percent.mt <- percent.mt #Add to the metadata of seurat object

#Violin plot raw data
vln_plot <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
title <- ggdraw() + draw_label("Quality Control GSE132771", fontface = 'bold')
vln_plot <- plot_grid(title, vln_plot, ncol = 1, rel_heights = c(0.1, 1))
vln_plot
pdf(file = 'QC-violin-plot.pdf', width = 12)
print(vln_plot)
dev.off()

#Scatter plot raw data
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter_plot <- CombinePlots(plots = list(plot1, plot2))
  title <- ggdraw() + draw_label("Quality Control GSE132771", fontface='bold')
  scatter_plot <- plot_grid(title, scatter_plot, ncol=1, rel_heights=c(0.1, 1))
  scatter_plot
pdf(file = 'QC-scatterplot.pdf', width = 12)
print(scatter_plot)
dev.off()




#Calculate the median absolute deviation (MAD) for library size, detected genes and mtDNA
mad_library_size <- mad(data$nCount_RNA)
mad_detected_genes <- mad(data$nFeature_RNA)
mad_mt_proportion <- mad(data$percent.mt)

#Calculate the median values
median_library_size <- median(data$nCount_RNA)
median_detected_genes <- median(data$nFeature_RNA)
median_percent_mt <- median(data$percent.mt)

# Apply filtering to exclude cells with extreme MAD values
filtered_cells <- which(
  abs(data$nCount_RNA - median_library_size) < 5 * mad_library_size &
  abs(data$nFeature_RNA - median_detected_genes) < 5 * mad_detected_genes &
  abs(data$percent.mt - median(data$percent.mt)) < 5 * mad(data$percent.mt)
)

#OR
filtered_cells <- which(
    data$nFeature_RNA > 200 &
    data$percent.mt < 20
)

# Subset the Seurat object to keep only the filtered cells
filtered_data <- data[filtered_cells, ]


# Print the number of cells before and after filtering
cat("Number of cells before filtering:", nrow(data), "\n")
cat("Number of cells after filtering:", nrow(filtered_data), "\n")

#Save the filtered object
saveRDS(filtered_data, file = "/home/jowilson/cosbi/tsukui/seurat_GSE132771_f.rds")



