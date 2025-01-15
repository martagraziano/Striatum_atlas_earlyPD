#Trajectory analysis
suppressPackageStartupMessages({
  library(Seurat)
  library(plotly)
  options(rgl.printRglwidget = TRUE)
  library(Matrix)
  install.packages("sparseMatrixStats")
  library(sparseMatrixStats)
  BiocManager::install("slingshot")
  library(slingshot)
  BiocManager::install("tradeSeq")
  library(tradeSeq)
  library(patchwork)
  BiocManager::install("DelayedMatrixStats")
  library(DelayedMatrixStats)
  BiocManager::install("scran")
  library(scran)
  library(readxl)
})

plot_directory<- "Z:/dmclab/Marta/PD/snRNA human putamen/new_allen_clustering/Pseudotime"
dir.create(plot_directory)
# Calculate cluster centroids (for plotting the labels later)
rna.combined_sub<- readRDS( "Z:/dmclab/Marta/PD/snRNA human putamen/data/rna.combined_res1_annotated_2_clinicaldata_ann_allen.rds")
rna.combined_sub <- subset(x = rna.combined_sub, idents = setdiff(unique(Idents(rna.combined_sub)), c("MGE Int","Deep layer CT","LAMP5-LHX6 Int","CGE Int","Splatter")))
rna.combined_sub$clusters <- Idents(rna.combined_sub)
rna.combined_sub <- subset(rna.combined_sub, subset = Braak != 3)


clinical_data<- read_xlsx("Z:/dmclab/Marta/PD/snRNA human putamen/clinical_data.xlsx")

# Extract the metadata from the Seurat object
rna_metadata <- rna.combined_sub@meta.data

# Create a subset of clinical_data that includes only relevant columns
# Ensure clinical_data and rna_metadata have a common sample identifier
clinical_subset <- clinical_data[, c("SampleID", "Sex","Age", "Braak", "PMI")]



# Merge the clinical data with Seurat metadata based on 'SampleID'
# This will match rows where SampleID matches between the two data frames
merged_metadata <- merge(rna_metadata, clinical_subset, by = "SampleID", all.x = TRUE)

# Add the updated metadata back to the Seurat object
rna.combined_sub <- AddMetaData(rna.combined_sub, metadata = merged_metadata)

mm <- sparse.model.matrix(~ 0 + factor(rna.combined_sub$Braak))
colnames(mm) <- levels(factor(rna.combined_sub$Braak))
centroids2d <- as.matrix(t(t(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings) %*% mm) / Matrix::colSums(mm))

cl.colors<- c("#0a4347ff","#e080b9","#bc7ed4","#99B9B5","#f3c9d9","#b6dddb","#a2cae9","#f69e92" )


cl.colors <-  c("#F5D3C8","#c29e93ff","#8d7d93ff","#558877ff","#5566aaff","#fccde5ff","#8d7d93ff","#484444ff","#856fa7ff","#b89a01ff","#5f070cff","#102c52ff","#b6deeaff","#481518ff")
braak_colors <- c(
  "0" = "gray", 
  
  "3" = "#D9BFBF", 
  "4" = "#a34246ff", 
  "5" = "#79090fff", 
  "6" = "#330000"
)

vars <- c("Braak", "Sample_group", "clusters")
pl <- list()

for (i in vars) {
  if (i == "Braak"){
  pl[[i]] <- DimPlot(rna.combined_sub, group.by = i, reduction = "umap_after_harmony_4",   label = T) + NoLegend()
}else{
  if (i=="clusters"){
    pl[[i]] <- DimPlot(rna.combined_sub, group.by = i, reduction = "umap_after_harmony_4",  label = T) + NoLegend()
  }else{
    pl[[i]] <- DimPlot(rna.combined_sub, group.by = i, reduction = "umap_after_harmony_4", label = T) + NoLegend()
  }
 } 
}

wrap_plots(pl)


ENDS<- c("6")



set.seed(1)
lineages <- as.SlingshotDataSet(getLineages(
  data           =rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings,
  clusterLabels  = rna.combined_sub$Braak,
  dist.method    = "slingshot", # It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
  end.clus       = ENDS, # You can also define the ENDS!
  start.clus     = "0"
)) # define where to START the trajectories


# IF NEEDED, ONE CAN ALSO MANULALLY EDIT THE LINEAGES, FOR EXAMPLE:
# sel <- sapply( lineages@lineages, function(x){rev(x)[1]} ) %in% ENDS
# lineages@lineages <- lineages@lineages[ sel ]
# names(lineages@lineages) <- paste0("Lineage",1:length(lineages@lineages))
# lineages


# Change the reduction to our "fixed" UMAP2d (FOR VISUALISATION ONLY)
# Specify the PDF file name and path
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired

# Generate the plot
plot(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings, 
     col=adjustcolor(braak_colors, alpha.f = 0.5), 
     cex = .5, pch = 16)

# Add lines
lines(lineages, lwd = 2, col = "black")

# Add text labels
text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 5, col = "white")

# Close the PDF device to finalize the file
dev.off()


# Define curves
curves <- as.SlingshotDataSet(getCurves(
  data          = lineages,
  thresh        = 1e-1,
  stretch       = 1e-1,
  allow.breaks  = F,
  approx_points = 100
))


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_curves.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired

  plot(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings,col=cl.colors[rna.combined_sub$clusters], cex = .5, pch = 16)
  lines(curves, lwd = 2, col = "black")
  text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)
dev.off()

pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

# Convert pseudotime to a dataframe for easier manipulation
pseudotime_df <- as.data.frame(pseudotime)

x <- rowMeans(pseudotime)
x <- x / max(x)
o <- order(x)

# Extract curves from the Slingshot object
curves1 <- slingshot::slingCurves(curves)  # Ensure this is your Slingshot object with curves
library(ggplot2)
library(dplyr)

# Prepare data for ggplot
curves_data <- do.call(rbind, lapply(seq_along(curves1), function(i) {
  data.frame(
    x = curves1[[i]]$s[, 1],
    y = curves1[[i]]$s[, 2],
    Lineage = as.factor(i)
  )
}))

# Define colors for the curves
curve_colors <- c("red", "purple4")  # Add more colors if needed


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_gray.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6) 

# Define stages for each lineage
stages <- c("Lineage1")

# Plot UMAP with points colored by pseudotime
plot(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings,
     main = paste0("Pseudotime"), pch = 16, cex = 0.4, axes = TRUE, 
     xlab = "UMAP 1", ylab = "UMAP 2",
     col = "gray" , alpha=0.5 # Keep the points gray
)

# Overlay Slingshot curves with specified colors
for (i in seq_along(curves1)) {
  lines(curves1[[i]]$s, col = curve_colors[i %% length(curve_colors) + 1], lwd = 2)  # Use colors cyclically
}

# Create labels for the legend combining lineage and stages
legend_labels <- paste("Lineage", 1:length(curves1), "(", stages, ")", sep = " ")

# Add legend for curves
legend("bottomright", legend = legend_labels, 
       col = curve_colors[1:length(curves1)], lwd = 2, cex = 0.8, box.col = "transparent")
dev.off()




# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_pseudotime_.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6) 



#plotGeneCount(curves, clusters = rna.combined_sub$clusters, models = sceGAM)
  plot(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings[o, ],
       main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = T, xlab = "umap_after_harmony_4 1", ylab = "umap_after_harmony_4 2",
       col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]   
  )
  #points(centroids2d, cex = 2.5, pch = 16, col = "#FFFFFF99")
  #text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)
dev.off()






# Convert pseudotime to a dataframe for easier manipulation
pseudotime_df <- as.data.frame(pseudotime)
pseudotime_df$mean<- rowMeans(pseudotime)
pseudotime_df$cells<- rownames(pseudotime_df)



# Extract UMAP coordinates
umap_coordinates <- rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings
umap_df <- as.data.frame(umap_coordinates)
umap_df$Cluster <- as.factor(rna.combined_sub$clusters)  # Add cluster information if needed
umap_df$cells<- rownames(umap_df)

merged_df <- merge(umap_df, pseudotime_df, by = "cells", all = TRUE)
#merged_df<- merge(merged_df, pseudo_df, by = "cells", all = TRUE)


merged_df_long<-merged_df  %>% pivot_longer(cols=c('Lineage1',"Lineage2"),
                                            names_to='Lineage')



# Plots
plot_filename <- paste0(plot_directory, "/pseudotime_percluster.pdf")


# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6) 

# Create the jittered point plot with specified color gradient
ggplot(merged_df_long,aes(x =value , y = Cluster, color = value)) +
  geom_jitter(size = 2, width = 0.2, height = 0.2, alpha = 0.6) +  # Jittered points
  scale_color_gradientn(colors = c("grey70", "orange3", "firebrick", "purple4")) +  # Custom gradient
  labs(title = "Cell Distribution Across Lineages and Clusters",
       x = "Lineage",
       y = "Cluster",
       color = "Pseudotime") +
  facet_grid(.~Lineage)+
  theme_classic() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Angled x-axis text for
dev.off()


# Step 1: Aggregate the data
# This code aggregates the merged dataset by grouping it based on Lineage and Cluster.
# For each unique combination of Lineage and Cluster, it computes:
# - CellCount: The total number of cells that belong to that combination, calculated using n().
# - MeanPseudotime: The average pseudotime value of the cells in that combination, ignoring any missing values (NA).
# The result is a new data frame (aggregated_df) that contains the number of cells and the average pseudotime 
# for each lineage-cluster combination, facilitating visualization and further analysis.


aggregated_df <- merged_df_long %>%
  group_by(Lineage, Cluster) %>%
  summarize(
    CellCount = n(),  # Count the number of cells
    MeanPseudotime = mean(value, na.rm = TRUE)  # Calculate the mean pseudotime
  )


# Step 2: Create the plot
# Plots
plot_filename <- paste0(plot_directory, "/pseudotime_percluster_norm_nolineage_bar.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=7, height=4) 
ggplot(aggregated_df, aes(x = MeanPseudotime, y = Cluster, fill=MeanPseudotime)) +
  geom_bar(stat="identity") + 
  
  #Use points to represent each cluster-lineage combination
  scale_fill_gradientn(colors =c("grey70", "orange3", "firebrick", "purple4")) +  # Custom color gradient
  labs(title = "Cell Distribution Across Lineages and Clusters",
       x = "Lineage",
       y = "Cluster",
       size = "Number of Cells",
       color = "Mean Pseudotime") +
  #xlim(0,30)+
  #scale_size(range = c(1,5))+
  facet_grid(.~Lineage)+
  theme_classic() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()
#tradeSeq is a recently proposed algorithm to find trajectory differentially expressed genes. 
#It works by smoothing the gene expression along the trajectory by fitting a smoother using generalized 
#additive models (GAMs), and testing whether certain coefficients are statistically different between 
#points in the trajectory.
BiocParallel::register(BiocParallel::SnowParam())



sel_cells <- split(colnames(rna.combined_sub@assays$RNA@data), rna.combined_sub$Braak)

sel_cells <- sel_cells[names(sel_cells) != "3"]
#randomly selecting the nmber of cells that is in the smallest cluster
sel_cells <- unlist(lapply(sel_cells, function(x) {
  set.seed(1)
  return(sample(x, 1113))
}))


gv <- as.data.frame(na.omit(scran::modelGeneVar(rna.combined_sub@assays$RNA@data[, sel_cells])))
gv <- gv[order(gv$bio, decreasing = T), ]
sel_genes <- sort(rownames(gv)[1:3000])

path_file <- paste0(plot_directory,"/rna.combined_subclustering_trajectory_scegam.rds")
sceGAM<- readRDS(paste0(plot_directory,"/rna.combined_subclustering_trajectory_scegam.rds"))
# fetch_data is defined at the top of this document

  sceGAM <- fitGAM(
    counts = drop0(rna.combined_sub@assays$RNA@data[sel_genes, sel_cells]),
    pseudotime = pseudotime[sel_cells, ],
    cellWeights = cellWeights[sel_cells, ],
    nknots = 10, verbose = T, parallel = T, sce = TRUE,
    BPPARAM = BiocParallel::SnowParam()
  )
  saveRDS(sceGAM, file="Z:/dmclab/Marta/PD/snRNA human putamen/new_allen_clustering/subclustering/Pseudotime/rna.combined_subclustering_trajectory_scegam.rds")

  
  # Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_curves_2.pdf")
  
  # Open the PDF device to save the plot
pdf(plot_filename, width=7, height=6)  # Adjust width and height as desired

plotGeneCount(curves, clusters = rna.combined_sub$clusters,  models = sceGAM)+
  scale_color_manual(values=cl.colors) 
dev.off()

lc <- sapply(lineages@lineages, function(x) {
  rev(x)[1]
})
names(lc) <- gsub("Lineage", "L", names(lc))

# Define some color palette



# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_lineage.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired
  plot(rna.combined_sub@reductions$umap_after_harmony_3@cell.embeddings, col = cl.colors[rna.combined_sub$clusters], pch = 16)
  lines(curves, lwd = 2, col =c("red","blue"))
  points(centroids2d[lc, ], col =c("red","blue"), pch = 16, cex = 4)
  text(centroids2d[lc, ], labels = names(lc), cex = 1, font = 2, col = "white")
dev.off()



#Genes that change with pseudotime
#We can first look at general trends of gene expression across pseudotime.

set.seed(8)
res <- na.omit(associationTest(sceGAM, contrastType = "consecutive"))
res <- res[res$pvalue < 1e-3, ]
res <- res[res$waldStat > mean(res$waldStat), ]
res <- res[order(res$waldStat, decreasing = T), ]
write.csv(res, file=paste0(plot_directory, "/Genes_pseudotime.csv"))
res$Gene<- rownames(res)
res1<- res[!grepl("^GM", res$Gene) & !grepl("RIK$", res$Gene)  & !grepl("PTGDS", res$Gene) & !grepl("^AC", res$Gene)& !grepl("^AL", res$Gene) & !grepl("^AP", res$Gene)& !grepl("^LINC", res$Gene),]
res[1:10, ]
# Set up the plotting window

# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_top15_trajectory_de.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=16, height=12) 

par(mfrow = c(4, 4), mar = c(.1, .1, 2, 1)) 


  plot(rna.combined_sub@reductions$umap_after_harmony_4@cell.embeddings, 
       col = cl.colors[rna.combined_sub$clusters], 
       cex = .5, pch = 16, axes = F, xlab = "", ylab = "")
  
  lines(curves, lwd = 2, col = "black")
  
  # Plot centroids as black squares
  points(centroids2d[lc, ], col = "black", pch = 15, cex = 3, xpd = T)
  
  # Add labels to the centroids
  text(centroids2d[lc, ], labels = names(lc), cex = 1, font = 2, col = "white", xpd = F)


# Extract and filter the variable names
vars <- rownames(res1[1:30, ])
vars <- na.omit(vars[vars != "NA"])





# Loop through each variable (gene) to plot expression levels
for (i in vars) {
  # Extract the gene expression vector and normalize it
  x <- drop0(rna.combined_sub@assays$RNA3@data)[i, ]
  x <- (x - min(x)) / (max(x) - min(x))  # Normalize expression to range [0, 1]
  
  # Order cells based on normalized expression
  o <- order(x)
  
  # Plot UMAP embedding with the color-coded gene expression
  plot(rna.combined_sub@reductions$umap_after_harmony_3@cell.embeddings[o, ],
       main = paste0(i), pch = 16, cex = 0.5, axes = F, xlab = "", ylab = "",
       col = colorRampPalette(c("lightgray", "#B34B50", "#5F070C","#240705"))(99)[x[o] * 98 + 1])
}
dev.off()


res_plot<- subset(res1, rownames(res1) %in% vars)
umap_df$Cluster <- as.factor(rna.combined_sub$clusters) 
res_plot$genes<- rownames(res_plot)

lower_quantile <- quantile(res_plot$waldStat, 0.10, na.rm = TRUE)
upper_quantile <- quantile(res_plot$waldStat, 0.5, na.rm = TRUE)
plot_filename <- paste0(plot_directory, "/top_30_pseudotimeDE_WALD2.pdf")
# Open the PDF device to save the plot
pdf(plot_filename, width=4, height=6) 
ggplot(res_plot, aes(meanLogFC, reorder(Gene,meanLogFC), fill=waldStat))+
geom_bar(stat="identity")+
  scale_fill_gradientn(colors = c("orange3", "firebrick", "purple4"),limits = c(min(res_plot$waldStat), max(res_plot$waldStat))) +
xlim(0,1.5)+# Using RdPu palette
theme_classic()
dev.off()

#genes that change between 2 pseudotime

set.seed(8)
res <- startVsEndTest(sceGAM)
#Genes that change with pseudotime
#We can first look at general trends of gene expression across pseudotime.
res <- res[res$pvalue < 1e-3, ]
res <- res[res$waldStat > mean(res$waldStat), ]
res <- res[order(res$waldStat , decreasing = T), ]
write.csv(res, file="Z:/dmclab/Marta/PD/snRNA human putamen/new_allen_clustering/subclustering/Pseudotime/Genes_pseudotime_STARTVSEND.csv")
res<- read.csv("Z:/dmclab/Marta/PD/snRNA human putamen/new_allen_clustering/subclustering/Pseudotime/Genes_pseudotime_STARTVSEND.csv")
res$Gene<- rownames(res)
res1<- res[!grepl("^GM", res$Gene) & !grepl("RIK$", res$Gene)  & !grepl("PTGDS", res$Gene) & !grepl("^AC", res$Gene)& !grepl("^AL", res$Gene) & !grepl("^AP", res$Gene)& !grepl("^LINC", res$Gene),]

res1[1:10, ]
# Set up the plotting window
install.packages("fields")  # Install the package if you haven't already
library(fields) 

# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_top15_trajectory_de_startvsend.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=15, height=12) 

par(mfrow = c(6, 6), mar = c(.1, .1, 2, 1)) 


plot(rna.combined_sub@reductions$umap_after_harmony_3@cell.embeddings, 
     cex = .5, pch = 16, axes = F, xlab = "", ylab = "")

lines(curves, lwd = 2, col = "black")

# Plot centroids as black squares
points(centroids2d[lc, ], col = "black", pch = 15, cex = 3, xpd = T)

# Add labels to the centroids
text(centroids2d[lc, ], labels = names(lc), cex = 1, font = 2, col = "white", xpd = T)



# Extract and filter the variable names
vars <- rownames(res1[1:40, ])
vars <- na.omit(vars[vars != "NA"])





# Loop through each variable (gene) to plot expression levels
for (i in vars) {
  # Extract the gene expression vector and normalize it
  x <- drop0(rna.combined_sub@assays$RNA3@data)[i, ]
  x <- (x - min(x)) / (max(x) - min(x))  # Normalize expression to range [0, 1]
  
  # Order cells based on normalized expression
  o <- order(x)
  
  color_palette <- colorRampPalette(c("lightgray", "#B34B50", "#5F070C", "#240705"))(200)
  
  # Plot UMAP embedding with the color-coded gene expression
  plot(rna.combined_sub@reductions$umap_after_harmony_3@cell.embeddings[o, ],
       main = paste0(i), pch = 16, cex = 0.5, axes = F, xlab = "", ylab = "",
       col = color_palette[x[o] * 149 + 1])
  # Scale the gradient by specifying `legend.only = TRUE`
  
}
# image.plot(
#   legend.only = TRUE,
#   zlim = range(x),  # Adjust based on the range of values in `x`
#   col = color_palette,
#   legend.width = 1,  # Adjust width as needed
#   legend.shrink = 0.6,  # Adjust legend size as needed
#   horizontal = FALSE  # Horizontal or vertical gradient legend
# )
dev.off()




res_plot<- subset(res1, rownames(res1) %in% vars)
res_plot$genes<- rownames(res_plot)
res_plot <- merge(umap_df, res_plot, by = "cells", all = TRUE)

lower_quantile <- quantile(res_plot$waldStat, 0.25, na.rm = TRUE)
upper_quantile <- quantile(res_plot$waldStat, 0.75, na.rm = TRUE)
plot_filename <- paste0(plot_directory, "/top_30_pseudotimeDE_2_lineage1.pdf")
# Open the PDF device to save the plot
pdf(plot_filename, width=4, height=5) 

ggplot(res_plot, aes(logFClineage1, reorder(Gene,-logFClineage1 ), fill=waldStat))+
  geom_bar(stat="identity")+
  scale_fill_gradientn(
    colors = c( "orange3", "firebrick", "purple4"),  # More colors for more contrast
    values = scales::rescale(c(min(res_plot$waldStat), lower_quantile, median(res_plot$waldStat), upper_quantile, max(res_plot$waldStat)))  # Rescale values for color distribution
  ) +
  xlim(-10,8)+# Using RdPu palette
  theme_classic()
dev.off()
plot_filename <- paste0(plot_directory, "/top_30_pseudotimeDE_2_lineage2.pdf")
# Open the PDF device to save the plot
pdf(plot_filename, width=4, height=5)
ggplot(res_plot, aes(logFClineage2, reorder(Gene,-logFClineage2 ), fill=waldStat))+
  geom_bar(stat="identity")+
  scale_fill_gradientn(
    colors = c( "orange3", "firebrick", "purple4"),  # More colors for more contrast
    values = scales::rescale(c(min(res_plot$waldStat), lower_quantile, median(res_plot$waldStat), upper_quantile, max(res_plot$waldStat)))  # Rescale values for color distribution
  ) +
  xlim(-10,8)+# Using RdPu palette
  theme_classic()
dev.off()
genes<- c("DRD2","SLC24A4","ERBB4")

DefaultAssay(rna.combined_sub)<- "SCT"
lapply(unique(genes), function(gene){
  tryCatch({
    myplot <- FeaturePlot(rna.combined_sub, features = gene, reduction = "umap_after_harmony_4", slot = "data", pt.size=0.6,
                          label=FALSE, raster=FALSE) & 
      theme(legend.position = "right") &
      #scale_color_gradientn(colors = c("gray","#60795bff", "#354332ff"), limits=c(0,3))
      scale_color_gradientn(colors = c("gray","#B34B50", "#5F070C","#240705"),, limits=c(0,6.5))
    
    
    plot_filename <- paste0(plot_directory, "/marker_umap_HD_norm_fig_sized_", gene, ".pdf")
    pdf(plot_filename, width = 10, height = 9)  # Adjust width and height as desired
    print(myplot)
    dev.off()

  }, error = function(e) {
    message(paste("Error with gene:", gene, " - Skipping to next..."))
    # Continue with the next gene if there's an error
  })
  

})

