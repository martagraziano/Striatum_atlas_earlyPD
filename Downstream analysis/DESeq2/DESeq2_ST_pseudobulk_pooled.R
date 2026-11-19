# =============================================================================
# Pseudobulk DESeq2 - spatial transcriptomics, mild 6-OHDA vs control
# Graziano et al., "A molecular atlas of cell-type specific signatures in the
# parkinsonian striatum"
#
# Raw Spatial-assay UMI counts from the dorsal striatal ST clusters 4, 5 and 8
# are summed per mouse (one pseudobulk sample per animal; n = 3 per group).
# DESeq2, design ~ condition (Control = reference), Wald test, BH correction.
# Significant: padj < 0.05 and |log2FC| > 0.25.
#
# Input : integrated Visium Seurat object (Idents = ST clusters)
# Output: DESeq2_ST_Pooled_striatal_clusters_all_genes.csv
#         DESeq2_ST_Pooled_striatal_clusters.csv  (significant genes;
#         Supplementary Table 5)
# =============================================================================

library(Seurat)
library(DESeq2)
library(dplyr)

# ---- Paths and parameters ---------------------------------------------------
seurat_rds_path <- "data/data.integrated.clusters_2022-08-22.rds"
out_dir         <- "results/DESeq2_ST"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

condition_col  <- "sample_type"      # "Treated" / "Control"
mouse_id_col   <- "mouse_id"         # one pseudobulk sample per mouse
treated_label  <- "Treated"
control_label  <- "Control"
pooled_clusters <- c("4", "5", "8")  # DLS (4) and DMS (5, 8) striatal clusters

padj_cut <- 0.05
lfc_cut  <- 0.25

# ---- Load object ------------------------------------------------------------
st_object <- readRDS(seurat_rds_path)
DefaultAssay(st_object) <- "Spatial"

# ---- Pseudobulk: sum counts of all spots in clusters 4, 5, 8 per mouse -------
spots      <- colnames(st_object)[as.character(Idents(st_object)) %in% pooled_clusters]
meta       <- st_object@meta.data[spots, ]
raw_counts <- if (packageVersion("SeuratObject") >= "5.0.0") {
  LayerData(st_object, assay = "Spatial", layer = "counts")[, spots]      # Seurat v5
} else {
  GetAssayData(st_object, assay = "Spatial", slot = "counts")[, spots]    # Seurat v4
}

mice <- unique(meta[[mouse_id_col]])
pb_counts <- sapply(mice, function(m) {
  spots_m <- rownames(meta)[meta[[mouse_id_col]] == m]
  rowSums(raw_counts[, spots_m, drop = FALSE])
})
colnames(pb_counts) <- mice

sample_meta <- data.frame(
  condition = sapply(mice, function(m)
    unique(meta[[condition_col]][meta[[mouse_id_col]] == m])[1]),
  row.names = mice
)
sample_meta$condition <- factor(sample_meta$condition,
                                levels = c(control_label, treated_label))
print(sample_meta)

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = pb_counts[, rownames(sample_meta)],
                              colData   = sample_meta,
                              design    = ~ condition)
dds <- dds[rowSums(counts(dds)) > 0, ]     # remove genes with zero counts
dds <- DESeq(dds)

res <- results(dds,
               contrast      = c("condition", treated_label, control_label),
               pAdjustMethod = "BH")

res_df <- as.data.frame(res) %>%
  mutate(gene = rownames(.), cluster = "Striatum_pooled") %>%
  filter(!is.na(padj)) %>%
  arrange(padj, desc(abs(log2FoldChange)))

sig <- res_df %>% filter(padj < padj_cut, abs(log2FoldChange) > lfc_cut)

write.csv(res_df, file.path(out_dir, "DESeq2_ST_Pooled_striatal_clusters_all_genes.csv"),
          row.names = FALSE)
write.csv(sig,    file.path(out_dir, "DESeq2_ST_Pooled_striatal_clusters.csv"),
          row.names = FALSE)

cat("Significant genes (padj <", padj_cut, ", |log2FC| >", lfc_cut, "):", nrow(sig), "\n")
sessionInfo()
