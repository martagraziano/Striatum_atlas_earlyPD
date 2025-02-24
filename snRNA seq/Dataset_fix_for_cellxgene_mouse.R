metadata<- rna.combined@meta.data

# Step 2: Update the metadata in the Seurat object
rna.combined@meta.data <- metadata

# View the updated Seurat object metadata to confirm
meta_head<- head(rna.combined@meta.data)

samples<- as.data.frame(cbind(metadata$Sample, metadata$SampleID, metadata$Sample_group, metadata$Renamed_Sample))
colnames(samples)<- c("Sample", "SampleID","Sample_group")
unique_samples <- samples %>%
    distinct()

clusters_allen<- as.data.frame(cbind(as.character(metadata$clusters), metadata$cluster_label, metatdata$class))
unique_samples <- clusters_allen %>%
  distinct()

colnames(unique_samples)<- c("Sample", "SampleID","Sample_group")
write.csv(meta_head,"Z:/dmclab/Marta/PD/Combined_snRNAseq/data/Samples_metadata_head.csv")

# Create a mapping of the new names
rename_mapping <- c(
  "15a" = "MP_15A",
  "15b" = "MP_15B",
  "15c" = "MP_15C",
  "15d" = "MP_15D",
  "18a" = "MP_18A",
  "18d" = "MP_18B",
  "C4" = "Ctr_4",
  "C5" = "Ctr_5",
  "C1" = "Ctr_1",
  "C2" = "Ctr_2",
  "C3" = "Ctr_3",
  "MP1" = "MP_11A",
  "MP2" = "MP_11B",
  "MP3" = "MP_11C",
  "6OHDA_1" = "6OHDA_1",
  "6OHDA_2" = "6OHDA_2",
  "6OHDA_3" = "6OHDA_3",
  "Control_1" = "Ctr_6",
  "Control_2" = "Ctr_7",
  "Control_3" = "Ctr_8"
)

# Add a new column with renamed samples
samples <- samples %>%
  mutate(Renamed_Sample = rename_mapping[Sample])
# View the updated dataframe
View(samples)


# Step 3: Merge MapMyCells annotations into Seurat metadata
rna.combined <- AddMetaData(
  object = rna.combined,
  metadata = samples
)




# Rename columns
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "SampleID"] <- "Experiment_ID"
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "Sample_name"] <- "Sample"
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "Sample"] <- "Sample_name_old"


# Specify the columns to remove
columns_to_remove <- c("Sample_name_old") # Replace with actual column names

# Remove the specified columns
rna.combined@meta.data <- rna.combined@meta.data %>%
  select(-all_of(columns_to_remove))


#save with new names
saveRDS(rna.combined, "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_allen_newcols.rds")

rna.combined<- readRDS("Z:/dmclab/Marta/PD/Combined_snRNAseq/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_allen_newcols.rds")




colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "Sample"] <- "donor_id"
colnames(rna.combined@meta.data)[colnames(rna.combined@meta.data) == "Sample_group"] <- "batch_condition" 

rna.combined$development_stage_ontology_term_id <- colnames(rna.combined)

rna.combined$development_stage_ontology_term_id[rna.combined$donor_id %in% c("MP_15A", "MP_15B", "MP_15C","MP_15D")] <- "MmusDv_0000161"
rna.combined$development_stage_ontology_term_id[rna.combined$donor_id %in% c("MP_18A", "MP_18B")] <- "MmusDv_0000164"
rna.combined$development_stage_ontology_term_id[rna.combined$donor_id %in% c("MP_11A", "MP_11B", "MP_11C")] <- "MmusDv_0000156"
rna.combined$development_stage_ontology_term_id[rna.combined$donor_id %in% c("Ctr_1", "Ctr_2", "Ctr_3", "Ctr_4", "Ctr_5", "Ctr_6", "Ctr_7" ,"Ctr_8")] <- "MmusDv_0000137"
rna.combined$development_stage_ontology_term_id[rna.combined$donor_id %in% c("6OHDA_1", "6OHDA_2", "6OHDA_3")] <- "MmusDv_0000137"


rna.combined$sex_ontology_term_id <- colnames(rna.combined)

rna.combined$sex_ontology_term_id[rna.combined$donor_id %in% c("MP_15A", "MP_15B", "MP_15C","MP_18A", "MP_18B","MP_11A","Ctr_2","Ctr_5")] <- "PATO:0000383"
rna.combined$sex_ontology_term_id[rna.combined$donor_id %in% c("6OHDA_1", "6OHDA_2", "6OHDA_3","Ctr_1","Ctr_3","Ctr_6", "Ctr_4", "Ctr_7" ,"Ctr_8","MP_15D","MP_11B", "MP_11C")] <- "PATO:0000384"

rna.combined$tissue_type <- colnames(rna.combined)
rna.combined$tissue_type <- "tissue"

rna.combined$tissue_ontology_term_id <- colnames(rna.combined)
rna.combined$tissue_ontology_term_id <-"UBERON_0002435"


rna.combined$suspension_type <- colnames(rna.combined)
rna.combined$suspension_type <-"nucleus"

#reduce the dimension of object deleting the reductions and normalizations that are not necessary
rna.combined[["integrated"]] <- NULL
rna.combined[["harmony_sid_1"]] <- NULL
rna.combined[["harmony_sid_2"]] <- NULL

rna.combined[["umap_after_harmony"]] <- NULL
rna.combined[["umap_after_harmony_2"]] <- NULL

rna.combined[["harmony_sid"]] <- rna.combined[["harmony_sid_3"]]
rna.combined[["umap_after_harmony"]] <- rna.combined[["umap_after_harmony_3"]]

rna.combined[["harmony_sid_3"]] <- NULL
rna.combined[["umap_after_harmony_3"]] <- NULL

rna.combined$clusters <- Idents(rna.combined)

# add cell onology terms for the clusters


rna.combined$cell_type_ontology_term_id <- colnames(rna.combined)
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("oligodendrocytes")] <- "CL_0000128"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Striatum D1")] <- "PCL_0110121"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Striatum D2")] <- "PCL_0113356"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Splatter")] <- "unknown"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Astrocytes/Epen")] <- "PCL_0110030"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("eSPN")] <- "PCL_0051482"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Cortex CT−L6b")] <- "PCL_0051474"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Cortex IT−ET")] <- "PCL_0110001"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Interneurons")] <- "PCL_0113215"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Microglia")] <- "PCL_0117692"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Chat/Npy/Pval+ Interneurons")] <- "PCL_0112998"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("OPCs")] <- "PCL_0110031"
rna.combined$cell_type_ontology_term_id[rna.combined$clusters %in% c("Vascular")] <- "PCL_0110033"

rna.combined$clusters <- as.character(rna.combined$clusters)
#save rna assay as .h5ad
DefaultAssay(rna.combined) <- "RNA"
SeuratDisk::SaveH5Seurat(object = rna.combined, filename = "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/seurat_obj_full_mouse_forcellxegene_rna.h5Seurat")
SeuratDisk::Convert(source = "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/seurat_obj_full_mouse_forcellxegene_rna.h5Seurat", dest = "h5ad")

#save sct assay as .h5ad

DefaultAssay(rna.combined) <- "SCT"
SeuratDisk::SaveH5Seurat(object = rna.combined, filename = "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/seurat_obj_full_mouse_forcellxegene_sct.h5Seurat")
SeuratDisk::Convert(source = "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/seurat_obj_full_mouse_forcellxegene_sct.h5Seurat", dest = "h5ad")

#merge the two objects in python, from the file called: Conversion from anndata for cellxgene.ipynb


#save with new names
saveRDS(rna.combined, "Z:/dmclab/Marta/PD/Combined_snRNAseq/data/rna.combined.snRNA_seq_merged_42_res1.5_annotated_allen_cellxgene.rds")


library(patchwork)

# Create the plot
myplot <- FeaturePlot(
  rna.combined, 
  features = "Lhx6", 
  reduction = "umap_after_harmony", 
  slot = "data", 
  pt.size = 0.6, 
  label = TRUE, 
  raster = TRUE
)

# Arrange the plots in 2 rows and 3 columns
myplot + plot_layout(ncol = 4, nrow = 5)