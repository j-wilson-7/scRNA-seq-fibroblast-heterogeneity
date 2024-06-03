# scRNA-seq-fibroblast-heterogeneity
 Analysing publicly available scRNA-seq datasets and matching to bulk RNA-seq in pHLFs using Scissor tool to further understand fibroblast heterogeneity
1. Analysis of Tsukui et al. IPF scRNA-seq dataset (GSE132771) and extraction of fibroblasts and myofibroblasts to identify CTHRC1+ population
- Quality control and filtering
- Clustering, UMAP and DGE

2.  Scissor tool (https://github.com/sunduanchen/Scissor) was applied to map bulk RNA-seq datasets where mTORC1 and mTORC2 were depleted using CRISPR/Cas9 gene editing to pre-existing IPF scRNA-seq datasets. Fibroblasts from Habermann et al (GSE135893) and Adams et al. (GSE136831) were extracted and used as reference scRNA-seq datasets.
- Scissor analysis and visualisation of results
