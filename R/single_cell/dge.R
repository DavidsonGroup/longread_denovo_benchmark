library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(edgeR)

# setwd("~/lab_davidson/yan.a/pb/5m_stranded/bambu")

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
# dir <- 'isonform_corset_oarfish'

merged_gene <- readRDS(paste0(dir, '/merged_gene.rds'))

# clustering and integration
merged_gene <- AddMetaData(merged_gene,
                           metadata = data.frame(name = colnames(merged_gene)) %>%
                             separate(name, c('source','barcode'), sep = 'k_'))

## add original annotation using azimuth in full data
filtered_annotated <- readRDS("/vast/projects/lab_davidson/yan.a/pb/seurat/filtered_annotated.rds")

annotation <- filtered_annotated@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, c('rep','num','bc', 'suffix'), remove = F) %>%
  mutate(rowname2 = str_remove(rowname, '-1$'))
rownames(annotation) <- annotation$rowname2

annotation2 <- merged_gene@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, c('rep','num','bc'), remove = F)

newmeta <- annotation2 %>% 
  select(rowname) %>%
  left_join(annotation, by = c('rowname'='rowname2')) 

merged_gene <- AddMetaData(merged_gene,
                             metadata = newmeta)

# run standard anlaysis workflow
pdf(paste0(dir, '/gene_qc.pdf'), width = 8, height = 6)
FeatureScatter(merged_gene, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", log = T) +
  ggtitle('Scatter plot before filtering')
VlnPlot(merged_gene, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0, log = T) +
  ggtitle('Violin plot before filtering')

# filtered_gene <- subset(merged_gene, subset = nFeature_RNA > 500 & nCount_RNA < 25000)
# 1. Get the metadata from the failing object
meta_data <- merged_gene@meta.data

# 2. Find the cell barcodes that pass your filter
cells_to_keep <- rownames(meta_data)[meta_data$nFeature_RNA > 500 & meta_data$nCount_RNA < 25000]

# 3. Manually subset the Seurat object using the list of cell names
#    This uses a different (and more robust) subsetting method.
filtered_gene <- merged_gene[, cells_to_keep]

VlnPlot(filtered_gene, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = T) +
  ggtitle('Violin plot after filtering')
dev.off()

filtered_gene
# An object of class Seurat 
# 205239 features across 50793 samples within 1 assay 
# Active assay: RNA (22371 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

filtered_gene <- NormalizeData(filtered_gene)
filtered_gene <- FindVariableFeatures(filtered_gene)
filtered_gene <- ScaleData(filtered_gene)
filtered_gene <- RunPCA(filtered_gene, 
                        features = VariableFeatures(object = filtered_gene))

filtered_gene <- FindNeighbors(filtered_gene, dims = 1:10, reduction = "pca")
filtered_gene <- FindClusters(filtered_gene, resolution = 0.2, cluster.name = "unintegrated_clusters")

filtered_gene <- RunUMAP(filtered_gene, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

pdf(paste0(dir, '/gene_umaps.pdf'), width = 8, height = 6)
DimPlot(filtered_gene, reduction = "umap.unintegrated", group.by = 'source') +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
DimPlot(filtered_gene, reduction = "umap.unintegrated", group.by = 'unintegrated_clusters', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )

DimPlot(filtered_gene, reduction = "umap.unintegrated", group.by = 'SingleR.DICE.main', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
DimPlot(filtered_gene, reduction = "umap.unintegrated", group.by = 'predicted.celltype.l1', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
dev.off()

pdf(paste0(dir, '/gene_barplot.pdf'), width = 8, height = 5)
filtered_gene@meta.data %>%
  ggplot(aes(x = unintegrated_clusters, fill = SingleR.DICE.fine)) +
  geom_bar(position = 'fill')
filtered_gene@meta.data %>%
  ggplot(aes(x = unintegrated_clusters, fill = predicted.celltype.l1)) +
  geom_bar(position = 'fill')
dev.off()

# cluster 0 CD16Mo and cluster 2 CD14Mo
join <- filtered_gene
join[['RNA']] <- JoinLayers(join[['RNA']])

dge <- Seurat2PB(join, 
                 sample = 'source', 
                 cluster = "unintegrated_clusters")
rownames(dge) <- str_replace_all(rownames(dge), '-', '_')
dge$genes$gene <- str_replace_all(dge$genes$gene, '-', '_')

saveRDS(filtered_gene, paste0(dir, '/filtered_gene.rds'))
saveRDS(dge, paste0(dir, '/dgelist_gene.rds'))

cat('Finished generating gene DGE objects!\n')
