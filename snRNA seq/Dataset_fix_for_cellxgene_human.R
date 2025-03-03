
rna.combined<-readRDS("Z:/dmclab/Marta/PD/snRNA human putamen/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_final.rds")
rna.combined<-readRDS("Z:/dmclab/Marta/PD/snRNA human putamen/data/rna.combined_res1_sub_annotated_allen_noduplicates.rds")

metadata<- rna.combined@meta.data

# Rename columns
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "SampleID"] <- "donor_id"
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "Sample_group"] <- "batch_condition" 


#add clinical data
library(readxl)
clinical_data<- read_xlsx("Z:/dmclab/Marta/PD/snRNA human putamen/clinical_data.xlsx")



# Extract the metadata from the Seurat object
rna_metadata <- rna.combined@meta.data

# Create a subset of clinical_data that includes only relevant columns
# Ensure clinical_data and rna_metadata have a common sample identifier
clinical_subset <- clinical_data[, c("SampleID", "Sex","Age", "Braak", "PMI")]



# Merge the clinical data with Seurat metadata based on 'SampleID'
# This will match rows where SampleID matches between the two data frames
merged_metadata <- merge(rna_metadata, clinical_subset, by = "SampleID", all.x = TRUE)

# Add the updated metadata back to the Seurat object
rna.combined <- AddMetaData(rna.combined, metadata = merged_metadata)

#save with new names
saveRDS(rna.combined, "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_allen_newcols.rds")

rna.combined<- readRDS("Z:/dmclab/Marta/PD/Combined_snRNAseq/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_allen_newcols.rds")






rna.combined$development_stage_ontology_term_id <- colnames(rna.combined)
rna.combined$development_stage_ontology_term_id<- "HsapDv_0000227"


rna.combined$sex_ontology_term_id <- colnames(rna.combined)

rna.combined$sex_ontology_term_id[rna.combined$Sex %in% c("F")] <- "PATO:0000383"
rna.combined$sex_ontology_term_id[rna.combined$Sex %in% c("M")] <- "PATO:0000384"


rna.combined$tissue_type <- colnames(rna.combined)
rna.combined$tissue_type <- "tissue"

rna.combined$tissue_ontology_term_id <- colnames(rna.combined)
rna.combined$tissue_ontology_term_id <-"UBERON_0001874"


rna.combined$suspension_type <- colnames(rna.combined)
rna.combined$suspension_type <-"nucleus"


rna.combined$disease_ontology_term_id <- colnames(rna.combined)
rna.combined$disease_ontology_term_id[rna.combined$batch_condition %in% c( "PD")] <- "MONDO_0008199"
rna.combined$disease_ontology_term_id[rna.combined$batch_condition %in% c("Control")] <- "PATO_0000461"

rna.combined$ organism_ontology_term_id  <- colnames(rna.combined)
rna.combined$ organism_ontology_term_id  <- "NCBITaxon:9606"

#reduce the dimension of object deleting the reductions and normalizations that are not necessary
rna.combined[["integrated"]] <- NULL
rna.combined[["harmony_sid_1"]] <- NULL
rna.combined[["harmony_sid_2"]] <- NULL

rna.combined[["umap_after_harmony"]] <- NULL
rna.combined[["umap_after_harmony_2"]] <- NULL

rna.combined[["harmony_sid_3"]] <- NULL
rna.combined[["umap_after_harmony_3"]] <- NULL

rna.combined[["harmony_sid"]] <- rna.combined[["harmony_sid_4"]]
rna.combined[["umap_after_harmony"]] <- rna.combined[["umap_after_harmony_4"]]

rna.combined[["harmony_sid_4"]] <- NULL
rna.combined[["umap_after_harmony_4"]] <- NULL

rna.combined$clusters <- Idents(rna.combined)
rna.combined$clusters <- as.character(rna.combined$clusters)

# add cell onology terms for the clusters main object 
rna.combined$cell_type_ontology_term_id <- colnames(rna.combined)
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Oligodendrocytes")] <- "CL_0000128"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("COPs-Oligo")] <- "CL_4023059"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("SPNs")] <- "CL_1001474"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Splatter")] <- "unknown"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Astrocytes")] <- "CL_0000127"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("eSPNs")] <- "CL_4030057"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("INT Ache-Npy")] <- "CL_4042013"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("INT Ache-Npy-Pvalb-SSt")] <- "CL_4042013"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Upper Layer IT")] <- "CL_0011005"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Microglia")] <- "CL_0000129"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Midbrain INH")] <- "CL_4023079"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("OPCs")] <- "CL_0002453"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Vascular")] <- "CL_0002139"


# add cell onology terms for the clusters subclustered object 
rna.combined$cell_type_ontology_term_id <- colnames(rna.combined)
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("dSPNs Mtx")] <- "CL_4030043"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("dSPNs CHST9+")] <- "CL_4030054"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("dSPNs Patch")] <- "CL_4030048"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("iSPNs Patch")] <- "CL_4030049"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("eSPNs")] <- "CL_4030057"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("SPNs Mixed")] <- "CL_1001474"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("iSPNs Mtx")] <- "CL_4030047"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("ICj")] <- "CL_4030053"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Splatter")] <- "unknown"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Ctx IT")] <- "CL_4023008"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("INs")] <- "CL_0000099"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("iSPNs Patch HTR7+")] <- "CL_4030049"


#save rna assay as .h5ad
# Load the raw Seurat object
raw_data <- readRDS("Z:/dmclab/Marta/PD/snRNA human putamen/data/raw_combined_snRNAseq_2.rds")
duplicates<- c("Seq176_1", "Seq176_10"   ,"Seq176_11" ,"Seq176_12","Seq176_13" , "Seq176_14","Seq176_2","Seq176_3","Seq176_4","Seq176_5","Seq176_6" , "Seq176_7","Seq176_8", "Seq176_9","Seq176_15", "Seq176_16","Seq176_17","AA_ASAP143_ctrl_NP22-75_PUT")                                               
raw_data[["RNA3"]] <- as(object = raw_data[["RNA"]], Class = "Assay")
DefaultAssay(raw_data)<- "RNA3"
raw_data[['RNA']] = NULL


raw_data<-subset(raw_data, SampleID %in% duplicates,invert = TRUE)

raw_data$batch_condition <- colnames(raw_data)


# Dynamically assign groups based on substrings in SampleID
raw_data$batch_condition[grepl("ctrl", raw_data$SampleID, ignore.case = TRUE)] <- "Control"
raw_data$batch_condition[grepl("PD", raw_data$SampleID, ignore.case = TRUE)] <- "PD"

length(unique(raw_data$SampleID))
metadata<- raw_data@meta.data

# Filter the raw object to keep only cells in rna.combined
cells_to_keep <- colnames(rna.combined)
raw_data <- subset(raw_data, cells = cells_to_keep)

# Check number of cells in both objects
cat("Number of cells in rna.combined (SCT-normalized):", ncol(rna.combined), "\n")
cat("Number of cells in raw_data (filtered raw):", ncol(raw_data), "\n")

# Check if the cell barcodes are identical and in the same order
same_cells <- identical(colnames(rna.combined), colnames(raw_data))

if (same_cells) {
  cat("Cell barcodes are identical and in the same order.\n")
} else {
  cat("WARNING: Cell barcodes do NOT match or are in a different order.\n")
}


# Install biomaRt if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Connect to Ensembl BioMart for Mouse
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get mapping of Gene Symbols to Ensembl IDs
gene_symbols <- unique(c(rownames(raw_data), rownames(rna.combined)))
mapping <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                 filters = "hgnc_symbol",
                 values = gene_symbols,
                 mart = ensembl)

# Check the mapping
head(mapping)

# Create a named vector for mapping
gene_map <- setNames(mapping$ensembl_gene_id, mapping$hgnc_symbol)

# Check if all symbols have an Ensembl ID
no_match_raw <- setdiff(rownames(raw_data), names(gene_map))
cat("Number of gene symbols in raw_data without Ensembl ID:", length(no_match_raw), "\n")

# Keep the original symbols as a new column
raw_data[["RNA3"]]@meta.features$original_gene_symbol <- rownames(raw_data[["RNA3"]])

# Add Ensembl IDs as a new column called "EnsemblID"
raw_data[["RNA3"]]@meta.features$EnsemblID <- gene_map[rownames(raw_data[["RNA3"]])]



cat("Mapped gene symbols to Ensembl IDs in raw_data, keeping original names.\n")

# Check if all symbols have an Ensembl ID
no_match_sct <- setdiff(rownames(rna.combined), names(gene_map))
cat("Number of gene symbols in rna.combined without Ensembl ID:", length(no_match_sct), "\n")

# Keep the original symbols as a new column
rna.combined[["RNA"]]@meta.features$original_gene_symbol <- rownames(rna.combined[["RNA"]])
rna.combined[["SCT"]]@meta.features$original_gene_symbol <- rownames(rna.combined[["SCT"]])

# Add Ensembl IDs as a new column called "EnsemblID"
rna.combined[["RNA"]]@meta.features$EnsemblID <- gene_map[rownames(rna.combined[["RNA"]])]
rna.combined[["SCT"]]@meta.features$EnsemblID <- gene_map[rownames(rna.combined[["SCT"]])]

cat("Mapped gene symbols to Ensembl IDs in rna.combined, keeping original names.\n")
# Check for the columns in raw_data
cat("Meta features in raw_data:\n")
print(head(raw_data[["RNA3"]]@meta.features))

# Check specifically for original_gene_symbol and EnsemblID
if (all(c("original_gene_symbol", "EnsemblID") %in% colnames(raw_data[["RNA3"]]@meta.features))) {
  cat("Both 'original_gene_symbol' and 'EnsemblID' exist in raw_data meta.features.\n")
} else {
  cat("WARNING: One or both columns are missing in raw_data!\n")
}

# Repeat for rna.combined (SCT object)
cat("\nMeta features in rna.combined (RNA):\n")
print(head(rna.combined[["RNA"]]@meta.features))

cat("\nMeta features in rna.combined (SCT):\n")
print(head(rna.combined[["SCT"]]@meta.features))

# Check specifically for original_gene_symbol and EnsemblID
if (all(c("original_gene_symbol", "EnsemblID") %in% colnames(rna.combined[["RNA"]]@meta.features))) {
  cat("Both 'original_gene_symbol' and 'EnsemblID' exist in rna.combined (RNA) meta.features.\n")
} else {
  cat("WARNING: One or both columns are missing in rna.combined (RNA)!\n")
}

if (all(c("original_gene_symbol", "EnsemblID") %in% colnames(rna.combined[["SCT"]]@meta.features))) {
  cat("Both 'original_gene_symbol' and 'EnsemblID' exist in rna.combined (SCT) meta.features.\n")
} else {
  cat("WARNING: One or both columns are missing in rna.combined (SCT)!\n")
}



# Export the filtered raw object as h5Seurat
DefaultAssay(raw_data)<- "RNA3"
SeuratDisk::SaveH5Seurat(object = raw_data, 
                         filename = "Z:/dmclab/Marta/PD/snRNA human putamen/data/seurat_obj_sub_human_forcellxegene_raw.h5Seurat")
SeuratDisk::Convert(source = "Z:/dmclab/Marta/PD/snRNA human putamen/data/seurat_obj_sub_human_forcellxegene_raw.h5Seurat", 
                    dest = "h5ad")


#save sct assay as .h5ad

DefaultAssay(rna.combined) <- "SCT"
SeuratDisk::SaveH5Seurat(object = rna.combined, filename = "Z:/dmclab/Marta/PD/snRNA human putamen/data/seurat_obj_sub_human_forcellxegene_sct_.h5Seurat")
SeuratDisk::Convert(source = "Z:/dmclab/Marta/PD/snRNA human putamen/data/seurat_obj_sub_human_forcellxegene_sct_.h5Seurat", dest = "h5ad")

#merge the two objects in python, from the file called: Conversion from anndata for cellxgene.ipynb

library(patchwork)

# Create the plot
myplot <- FeaturePlot(
  rna.combined, 
  reduction = "umap_after_harmony", feature="SLC32A1",
  pt.size = 0.6, slot="data",
  label = T, 
  raster = FALSE, 
)+ scale_color_gradientn(colors = c("gray","#B34B50", "#5F070C","#240705"))


# Arrange the plots in 2 rows and 3 columns
myplot + plot_layout(ncol = 4, nrow = 5)


# Run t-SNE using the PCA reduction
rna.combined <- RunTSNE(rna.combined, reduction = "harmony_sid", dims = 1:42)

# Confirm that t-SNE has been added
DimPlot(rna.combined, reduction = "tsne", cols= cl.colors)


