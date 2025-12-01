library(tidyverse)
library(Seurat)
library(Matrix)
library(patchwork)
library(edgeR)

# setwd("~/lab_davidson/yan.a/pb/5m_stranded/bambu")

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
# dir <- 'rnabloom2_corset_oarfish'

merged_tx <- readRDS(paste0(dir, '/merged_tx.rds'))
filtered_gene <- readRDS(paste0(dir, '/filtered_gene.rds'))

meta <- data.frame(name = colnames(merged_tx)) %>% 
  separate(name, c('source','barcode'), sep = 'k_', remove = F) %>%
  left_join(filtered_gene@meta.data %>% select(unintegrated_clusters, rowname),
            by = c('name'='rowname')) %>%
  rename(gene_cluster = unintegrated_clusters)

# clustering and integration
merged_tx <- AddMetaData(merged_tx,
                           metadata = meta)

## add original annotation using singleR in full data
filtered_annotated <- readRDS("/vast/projects/lab_davidson/yan.a/pb/seurat/filtered_annotated.rds")

annotation <- filtered_annotated@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, c('rep','num','bc', 'suffix'), remove = F) %>%
  mutate(rowname2 = str_remove(rowname, '-1$'))
rownames(annotation) <- annotation$rowname2

annotation2 <- merged_tx@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, c('rep','num','bc'), remove = F)

newmeta <- annotation2 %>% 
  select(rowname) %>%
  left_join(annotation, by = c('rowname'='rowname2')) 

merged_tx <- AddMetaData(merged_tx,
                           metadata = newmeta)

# run standard anlaysis workflow
pdf(paste0(dir, '/tx_qc.pdf'), width = 8, height = 6)

FeatureScatter(merged_tx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", log = T) +
  ggtitle('Scatter plot before filtering')
VlnPlot(merged_tx, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0, log = T) +
  ggtitle('Violin plot before filtering')

# filtered_tx <- subset(merged_tx, subset = nFeature_RNA > 500 & nCount_RNA < 25000)
# 1. Get the metadata from the failing object
meta_data <- merged_tx@meta.data

# 2. Find the cell barcodes that pass your filter
cells_to_keep <- rownames(meta_data)[meta_data$nFeature_RNA > 500 & meta_data$nCount_RNA < 25000]

# 3. Manually subset the Seurat object using the list of cell names
#    This uses a different (and more robust) subsetting method.
filtered_tx <- merged_tx[, cells_to_keep]

VlnPlot(filtered_tx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = T) +
  ggtitle('Violin plot before filtering')
dev.off()

filtered_tx
# An object of class Seurat 
# 205239 features across 50793 samples within 1 assay 
# Active assay: RNA (22371 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

filtered_tx <- NormalizeData(filtered_tx)
filtered_tx <- FindVariableFeatures(filtered_tx)
filtered_tx <- ScaleData(filtered_tx)
filtered_tx <- RunPCA(filtered_tx, 
                        features = VariableFeatures(object = filtered_tx))

filtered_tx <- FindNeighbors(filtered_tx, dims = 1:10, reduction = "pca")
filtered_tx <- FindClusters(filtered_tx, resolution = 0.2, cluster.name = "unintegrated_clusters")

filtered_tx <- RunUMAP(filtered_tx, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

pdf(paste0(dir, '/tx_umaps.pdf'), width = 8, height = 6)
DimPlot(filtered_tx, reduction = "umap.unintegrated", group.by = 'source') +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
DimPlot(filtered_tx, reduction = "umap.unintegrated", group.by = 'unintegrated_clusters', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )

DimPlot(filtered_tx, reduction = "umap.unintegrated", group.by = 'SingleR.DICE.main', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
DimPlot(filtered_tx, reduction = "umap.unintegrated", group.by = 'predicted.celltype.l1', label = T) +
  plot_layout(ncol = 1, nrow = 1,     
              widths = unit(3, "in"),   
              heights = unit(3, "in") )
dev.off()

pdf(paste0(dir, '/tx_barplot.pdf'), width = 8, height = 5)
filtered_tx@meta.data %>%
  ggplot(aes(x = unintegrated_clusters, fill = SingleR.DICE.fine)) +
  geom_bar(position = 'fill')
filtered_tx@meta.data %>%
  ggplot(aes(x = unintegrated_clusters, fill = predicted.celltype.l1)) +
  geom_bar(position = 'fill')
dev.off()

# add gene level clustering 


# cluster 0 CD16Mo and cluster 2 CD14Mo
join <- filtered_tx
join[['RNA']] <- JoinLayers(join[['RNA']])

dge <- Seurat2PB(join, 
                 sample = 'source', 
                 cluster = "gene_cluster")
rownames(dge) <- str_replace_all(rownames(dge), '-', '_')

rt_gene <- readRDS(paste0(dir, '/rt_gene.rds'))

dge$genes <- dge$genes %>%
  rownames_to_column(var = 'isoform') %>%
  left_join(rt_gene, by = 'isoform') %>%
  select(-gene)

saveRDS(filtered_tx, paste0(dir, '/filtered_tx.rds'))
saveRDS(dge, paste0(dir, '/dgelist_tx.rds'))

cat('Finished generating tx DGE objects!\n')
