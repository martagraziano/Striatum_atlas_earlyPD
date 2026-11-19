# =============================================================================
# GO Molecular Function enrichment of cell-type DEGs - mouse snRNA-seq
# Graziano et al., "A molecular atlas of cell-type specific signatures in the
# parkinsonian striatum"
#
# Input : per-cluster MAST DE tables (one sheet per cluster) produced by the
#         snRNA-seq clustering/DE workflow (Mouse full-dataset Rmd):
#           Combined-6OHDA_wt_DEGs-cluster_res05_nospl.xlsx     (mild 6-OHDA vs control)
#           Combined-mp15_18_wt_DEGs-cluster_res05_nospl.xlsx   (MitoPark 15-18 w vs control)
# Method: clusterProfiler::enrichGO (ont = "MF", org.Mm.eg.db, SYMBOL;
#         one-sided hypergeometric test, BH correction, default cut-offs
#         pvalueCutoff = 0.05), run separately on up- and downregulated DEGs
#         (Bonferroni-adjusted p < 0.0001; |log2FC| >= 0.2 already applied by
#         FindMarkers logfc.threshold = 0.2). The 40 terms with the highest gene
#         count are kept per cluster/direction.
# Output: GO_mostfreq_DEGs.xlsx (all models), GO_common_across_models_DEGs.xlsx,
#         GO_dotplot_<cluster>_<model>.pdf (Figure 2 / Figure S5 GO panels)
# =============================================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tibble)
library(openxlsx)
library(readxl)
library(ggplot2)

# ---- Paths ------------------------------------------------------------------
data_dir <- "data"          # folder with the MAST DE tables
out_dir  <- "results/GO"    # output folder
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

de_files <- c(
  "6OHDA"         = file.path(data_dir, "Combined-6OHDA_wt_DEGs-cluster_res05_nospl.xlsx"),
  "Mitopark15_18" = file.path(data_dir, "Combined-mp15_18_wt_DEGs-cluster_res05_nospl.xlsx")
)

# ---- Helpers ----------------------------------------------------------------
read_de_list <- function(file) {
  sheets <- excel_sheets(file)
  setNames(lapply(sheets, function(s) as.data.frame(read_xlsx(file, sheet = s))), sheets)
}

# enrichGO on one gene set; keep the 40 terms with the highest gene count
run_go <- function(genes, clust, direction) {
  GO_results_Df <- as.data.frame(
    enrichGO(gene = genes, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "MF"))
  if (nrow(GO_results_Df) == 0) return(NULL)
  GO_results_Df %>%
    slice_max(n = 40, order_by = as.numeric(Count)) %>%
    add_column(cluster = clust, reg = direction)
}

# ---- GO per model, cluster and direction ------------------------------------
Gene_ontology <- data.frame()

for (experiment in names(de_files)) {
  de_genes <- read_de_list(de_files[[experiment]])

  for (clust in names(de_genes)) {
    d <- de_genes[[clust]]
    genes_down <- d$Gene[d$p_val_adj < 0.0001 & d$avg_log2FC < 0]
    genes_up   <- d$Gene[d$p_val_adj < 0.0001 & d$avg_log2FC > 0]

    res <- bind_rows(run_go(genes_down, clust, "downreg"),
                     run_go(genes_up,   clust, "upreg"))
    if (nrow(res) > 0) {
      res$Experiment <- experiment
      Gene_ontology  <- rbind(Gene_ontology, res)
    }
  }
}
write.xlsx(Gene_ontology, file = file.path(out_dir, "GO_mostfreq_DEGs.xlsx"))

# GO terms found in both models
common_ids <- intersect(Gene_ontology$Description[Gene_ontology$Experiment == "6OHDA"],
                        Gene_ontology$Description[Gene_ontology$Experiment == "Mitopark15_18"])
Gene_ontology_common <- Gene_ontology[Gene_ontology$Description %in% common_ids, ]
write.xlsx(Gene_ontology_common, file = file.path(out_dir, "GO_common_across_models_DEGs.xlsx"))

# ---- Bubble plots (top 10 terms by gene count per cluster and model) --------
clusters    <- c("Cortex", "Striatum D1", "oligodendrocytes", "Striatum D2")
experiments <- unique(Gene_ontology$Experiment)

for (clust in clusters) {
  subset_cl <- subset(Gene_ontology, cluster == clust)
  if (nrow(subset_cl) == 0) next
  global_count_range <- range(subset_cl$Count)   # same bubble scale across models

  for (i in experiments) {
    subset_exp <- subset(subset_cl, Experiment == i) %>%
      slice_max(n = 10, order_by = as.numeric(Count))

    plot <- ggplot(subset_exp, aes(x = reorder(Description, -log10(p.adjust)),
                                   y = -log10(p.adjust), size = Count,
                                   color = reg == "downreg")) +
      geom_point() +
      scale_size_continuous(limits = global_count_range, range = c(3, 6)) +
      scale_y_continuous(limits = c(0, 7.5)) +
      scale_colour_manual(name = "Reg",
                          values = setNames(c("#31869bff", "#963634ff"), c(TRUE, FALSE))) +
      labs(title = paste("Gene ontology on cluster:", clust, i), size = "Number of Genes") +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_blank(),
            legend.position = "right")

    ggsave(file.path(out_dir, paste0("GO_dotplot_", gsub("/", "_", clust), "_", i, ".pdf")),
           plot = plot, width = 10, height = 6)
  }
}

sessionInfo()
