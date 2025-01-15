
quartz<-function(width,height){windows(width, height)}

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
})

#open clusters mouse
rna.combined_sub<- readRDS("Z:/dmclab/Marta/Combined_snRNAseq/data/rna.combined.subclusters_res1.5_annotated_allen_1122.rds")

rna.combined_sub <- subset(x = rna.combined, idents = setdiff(unique(Idents(rna.combined)), "Splatter"))
rna.combined_sub$clusters <- Idents(rna.combined_sub)


#only mitopark 
rna.filtered <- subset(rna.combined_sub, subset = Sample_group != "6OHDA")

# rna.combined$sample_group_2 <- colnames(rna.combined)
# rna.combined$sample_group_2[rna.combined$Sample %in% c("6OHDA_1", "6OHDA_2", "6OHDA_3")] <- "6OHDA"
# rna.combined$sample_group_2[rna.combined$Sample %in% c( "C4"      ,  "C5"   ,     "C1"   ,     "C2"     ,   "C3" ,"Control_1", "Control_2", "Control_3")] <- "Control"
# rna.combined$sample_group_2[rna.combined$Sample %in% c( "MP1"   ,    "MP2"     ,  "MP3"  )] <- "Mitopark_11"
# rna.combined$sample_group_2[rna.combined$Sample %in% c("15a"    ,   "15b"     ,  "15c"    ,   "15d" )] <- "Mitopark_15"
# rna.combined$sample_group_2[rna.combined$Sample %in% c("18a"  ,     "18d")] <- "Mitopark_18"

plot_directory<- "Z:/dmclab/Marta/Combined_snRNAseq/new_allen_based_clustering/Trajectory_analysis_mp"

dir.create(plot_directory)
cl.colors <-   c("#F5D3C8", "#558877ff","#0a4347ff","#b89a01ff","#B95F89","#484444ff","#dace7eff","#977482ff","#061c8dff","#FCCDE5","#856fa7ff","#102c52","#5F070C")

cl.colors <-  c("#0a4347ff","#B95F89","#fab9d7","#7fddc4ff","#856fa7ff")



library(RColorBrewer)

# Determine the number of clusters
n_clusters <- length(unique(rna.combined_sub$clusters))

# Get the "Paired" palette from RColorBrewer
palette <- brewer.pal(n = min(n_clusters, 12), name = "Paired")  # Adjust for max palette size

# Extend the palette if more than 12 clusters (optional)
if (n_clusters > 12) {
  palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_clusters)
}

vars <- c("Sample", "Sample_group", "clusters")
pl <- list()

for (i in vars) {
  if (i=="clusters"){
    pl[[i]] <- DimPlot(rna.filtered, group.by = i, reduction = "umap_after_harmony_3",  label = T) + NoLegend()
  }else{
    pl[[i]] <- DimPlot(rna.filtered, group.by = i, reduction = "umap_after_harmony_3", label = T) + NoLegend()
  }
 } 


wrap_plots(pl)



mm <- sparse.model.matrix(~ 0 + factor(rna.filtered$Sample_group))
colnames(mm) <- levels(factor(rna.filtered$Sample_group))
centroids2d <- as.matrix(t(t(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings) %*% mm) / Matrix::colSums(mm))

ENDS<- c("Mitopark15_18")



set.seed(1)
# lineages_mnn_mp <- as.SlingshotDataSet(getLineages(
#   data           =rna.filtered@reductions$umap_after_harmony_3@cell.embeddings,
#   clusterLabels  = rna.filtered$Sample_group,
#   dist.method    = "mnn", # It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
#   end.clus       = ENDS, # You can also define the ENDS!
#   start.clus     = "Control"
# )) # define where to START the trajectories


# Slingshot for Slingshot distance
lineages_slingshot <- as.SlingshotDataSet(getLineages(
  data = rna.filtered@reductions$umap_after_harmony_3@cell.embeddings,
  clusterLabels = rna.filtered$Sample_group,
  dist.method = "slingshot",
  start.clus = "Control",
  end.clus = "Mitopark15_18"
))


# If needed, manually assign clusters along the expected path
selected_clusters <- c("Control", "Mitopark11", "6OHDA", "Mitopark15_18")


# Assuming `lineages_mnn` and `lineages_slingshot` are your Slingshot datasets
# For MNN
plot(reducedDims(lineages_mnn_mp),
     pch = 16, asp = 1, main = "MNN Distance Method")
lines(SlingshotDataSet(lineages_mnn_mp), lwd = 2, col = "blue")

# For Slingshot distance
plot(reducedDims(lineages_slingshot), col = rna.combined_sub$Sample_group, 
     pch = 16, asp = 1, main = "Slingshot Distance Method")
lines(SlingshotDataSet(lineages_slingshot), lwd = 2, col = "red")


# Change the reduction to our "fixed" UMAP2d (FOR VISUALISATION ONLY)
# Specify the PDF file name and path
plot_filename <- paste0(plot_directory, "/UMAP_40PC_RES0.4_trajectory_1.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired

# Generate the plot
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings, 
     col=palette[rna.combined_sub$clusters], 
     cex = .5, pch = 16)

# Add lines
lines(lineages, lwd = 2, col = "black")

# Add text labels
#text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 4, col = "black")

# Close the PDF device to finalize the file
dev.off()


# Define curves
curves <- as.SlingshotDataSet(getCurves(
  data          = lineages_slingshot,
  thresh        = 1e-1,
  stretch       = 1e-1,
  allow.breaks  = F,
  approx_points = 100
))


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_40PC_RES0.4_trajectory_curves_mp.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired
# Plot the UMAP embedding
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings, 
     col = cl.colors[rna.filtered$clusters], cex = 0.5, pch = 16)

# Extract curve coordinates from Slingshot
curve_coords <- curves@curves[[1]]$s

# Draw the principal curve as a line
lines(curve_coords, col = "black", lwd = 2)

# Add a single arrow at the end of the curve
end_idx <- nrow(curve_coords)
arrows(curve_coords[end_idx - 1, 1], curve_coords[end_idx - 1, 2], 
       curve_coords[end_idx, 1], curve_coords[end_idx, 2], 
       col = "black", length = 0.2, lwd = 2)

# Overlay centroids
text(centroids2d, labels = rownames(centroids2d), cex = 0.5, font = 2)

dev.off()





pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)



x <- rowMeans(pseudotime)
x <- x / max(x)
o <- order(x)



# Extract curves from the Slingshot object
curves1 <- slingshot::slingCurves(curves)  # Ensure this is your Slingshot object with curves
library(ggplot2)
library(dplyr)

# Prepare data for ggplot (if needed for visualization)
curves_data <- do.call(rbind, lapply(seq_along(curves1), function(i) {
  data.frame(
    x = curves1[[i]]$s[, 1],
    y = curves1[[i]]$s[, 2],
    Lineage = as.factor(i)
  )
}))

# Define the gradient color palette
curve_gradient <- colorRampPalette(c("gray", "orange3", "firebrick", "purple4"))

# Number of gradient steps for the line segments
n_colors <- 100
gradient_colors <- curve_gradient(n_colors)

# Open the PDF device to save the plot
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_gradient_line.pdf")
pdf(plot_filename, width = 6, height = 6)

# Plot UMAP with points colored in gray
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, ],
     main = "Pseudotime", pch = 16, cex = 0.4, axes = TRUE,
     xlab = "UMAP 1", ylab = "UMAP 2",
     col = "gray"  # Keep the points gray
)

# Overlay Slingshot curves with gradient colors
for (i in seq_along(curves1)) {
  curve_coords <- curves1[[i]]$s
  n_points <- nrow(curve_coords)
  
  # Loop through each segment of the curve
  for (j in 1:(n_points - 1)) {
    # Determine the color for the segment
    color_index <- ceiling(j / n_points * n_colors)
    lines(curve_coords[j:(j + 1), ], col = gradient_colors[color_index], lwd = 8)
  }
}

# Close the PDF device
dev.off()




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
curve_colors <- c( "firebrick", "purple4")  # Add more colors if needed


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_gray_2.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6) 

# Define stages for each lineage
stages <- c( "Ctrl, MP15 ","Control and MP11")

# Plot UMAP with points colored by pseudotime
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, ],
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
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_pseudotime_1_b.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6) 



#plotGeneCount(curves, clusters = rna.combined_sub$clusters, models = sceGAM)
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, ],
     main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = T, xlab = "umap_after_harmony_4 1", ylab = "umap_after_harmony_4 2",
     col = colorRampPalette(c( "gray", "orange3", "firebrick", "purple4"))(110)[x[o] * 109 + 1]   
)
# points(centroids2d, cex = 1, pch = 16, col = "black")
# lines(curves, lwd = 2, col = "black")
# text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)
dev.off()




# Convert pseudotime to a dataframe for easier manipulation
pseudotime_df <- as.data.frame(pseudotime)
pseudotime_df$mean<- rowMeans(pseudotime)
pseudotime_df$cells<- rownames(pseudotime_df)



# Extract UMAP coordinates
umap_coordinates <- rna.filtered@reductions$umap_after_harmony_3@cell.embeddings
umap_df <- as.data.frame(umap_coordinates)
umap_df$Cluster <- as.factor(rna.filtered$clusters)  # Add cluster information if needed
umap_df$cells<- rownames(umap_df)

merged_df <- merge(umap_df, pseudotime_df, by = "cells", all = TRUE)
#merged_df<- merge(merged_df, pseudo_df, by = "cells", all = TRUE)


merged_df_long<-merged_df  %>% pivot_longer(cols=c('Lineage1',"Lineage2"),
                      names_to='Lineage')

# Plots
plot_filename <- paste0(plot_directory, "/pseudotime_percluster_lineage_2.pdf")

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
plot_filename <- paste0(plot_directory, "/pseudotime_percluster_norm_lineage_bar_2.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=4) 
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


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_pseudotime_.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired


#plotGeneCount(curves, clusters = rna.combined_sub$clusters, models = sceGAM)
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, ],
     main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = T, xlab = "umap_after_harmony_3 1", ylab = "umap_after_harmony_3 2",
     col = colorRampPalette(c("gray","orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]   
)
#points(centroids2d, cex = 2.5, pch = 16, col = "#FFFFFF99")
#text(centroids2d, labels = rownames(centroids2d), cex = 0.5, font = 2)
dev.off()




lc <- sapply(lineages_slingshot@lineages, function(x) {
  rev(x)[1]
})
names(lc) <- gsub("Lineage", "L", names(lc))


# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_lineage.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=6, height=6)  # Adjust width and height as desired
plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings, col = cl.colors[rna.filtered$clusters], pch = 16)
lines(curves, lwd = 2, col = "black")
points(centroids2d[lc, ], col = "black", pch = 16, cex = 4)
text(centroids2d[lc, ], labels = names(lc), cex = 1, font = 2, col = "white")
dev.off()

#tradeSeq is a recently proposed algorithm to find trajectory differentially expressed genes. 
#It works by smoothing the gene expression along the trajectory by fitting a smoother using generalized 
#additive models (GAMs), and testing whether certain coefficients are statistically different between 
#points in the trajectory.


BiocParallel::register(BiocParallel::SnowParam())

sel_cells <- split(colnames(rna.filtered@assays$RNA@data), rna.filtered$Sample_group)
#randomly selecting the nmber of cells that is in the smallest cluster
sel_cells <- unlist(lapply(sel_cells, function(x) {
  set.seed(1)
  return(sample(x, 10639))
}))




gv <- as.data.frame(na.omit(scran::modelGeneVar(rna.filtered@assays$RNA@data[, sel_cells])))
gv <- gv[order(gv$bio, decreasing = T), ]
sel_genes <- sort(rownames(gv)[1:6000])

path_file <- paste0(plot_directory, "/rna.combined_trajecotry_curves.rds")

# fetch_data is defined at the top of this document



sceGAM <- fitGAM(
  counts = drop0(rna.filtered@assays$RNA@data[sel_genes, sel_cells]),
  pseudotime = pseudotime[sel_cells, ],
  cellWeights = cellWeights[sel_cells, ],
  nknots = 10, verbose = T, parallel = T, sce = TRUE,
  BPPARAM = BiocParallel::SnowParam()
)
saveRDS(curves, path_file)

table(rowData(sceGAM)$tradeSeq$converged)



#A first exploration of the data analysis may consist of checking whether gene expression is associated with a particular lineage. 
#The statistical test performed here, implemented in the associationTest function, is testing the null hypothesis that all smoother 
#coefficients are equal to each other. This can be interpreted as testing whether the average gene expression is significantly changing along pseudotime.

assoRes <- associationTest(sceGAM)
head(assoRes)

startRes <- startVsEndTest(sceGAM)
#We can visualize the estimated smoothers for the third most significant gene.

oStart <- order(startRes$waldStat, decreasing = TRUE)
for (i in 1:15){
  sigGeneStart <- names(sceGAM)[oStart[i]]
  
  p<-plotSmoothers(sceGAM, drop0(rna.filtered@assays$RNA@data[sel_genes, sel_cells]), gene = sigGeneStart) +scale_color_manual(values= curve_colors)+ggtitle(paste0("Gene: ",sigGeneStart ))
  plot_filename <- paste0(plot_directory, "/estimated smoothers for_",sigGeneStart,".pdf")
  # Open the PDF device to save the plot
  pdf(plot_filename, width=8, height=6) 
  print(p)
  dev.off()
}

patternRes <- patternTest(sceGAM)
oPat <- order(patternRes$waldStat, decreasing = TRUE)

# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_top15_trajectory_patterndiff.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=15, height=12) 

par(mfrow = c(6, 6), mar = c(.1, .1, 2, 1)) 



# Loop through each variable (gene) to plot expression levels
for (i in 16:30) {
  sigGeneStart <- names(sceGAM)[oPat[i]]
  # Extract the gene expression vector and normalize it
  x <- drop0(rna.filtered@assays$RNA@data)[sigGeneStart, ]
  x <- (x - min(x)) / (max(x) - min(x))  # Normalize expression to range [0, 1]
  
  # Order cells based on normalized expression
  o <- order(x)
  
  color_palette <- colorRampPalette(c("lightgray", "#B34B50", "#5F070C", "#240705"))(200)
  
  umap_data <- data.frame(
    x = rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, 1],
    y = rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, 2],
    color = color_palette[x[o] * 149 + 1]
  )
  
  p1 <- ggplot(umap_data, aes(x = x, y = y, color = color)) +
    geom_point(size = 0.5) +
    scale_color_identity() +
    ggtitle(paste0(sigGeneStart)) +
    theme_void()
  # Scale the gradient by specifying `legend.only = TRUE`
  p<-plotSmoothers(sceGAM, drop0(rna.filtered@assays$RNA@data[sel_genes, sel_cells]), gene = sigGeneStart, alpha=0) +ggtitle(paste0("Gene: ",sigGeneStart ))
  p2<- p1+p
  
  plot_filename <- paste0(plot_directory, "/pattern_test_for_",sigGeneStart,".pdf")
  pdf(plot_filename, width=22, height=8) 
  print(p2)
  dev.off()
}



res <- startVsEndTest(sceGAM)
#Genes that change with pseudotime
#We can first look at general trends of gene expression across pseudotime.
res <- na.omit(associationTest(sceGAM, contrastType = "consecutive"))
res <- res[res$pvalue < 1e-3, ]
res <- res[res$waldStat > mean(res$waldStat), ]
res <- res[order(res$waldStat , decreasing = T), ]
write.csv(res, file=paste0(plot_directory, "/Genes_pseudotime.csv"))
res$Gene<- rownames(res)
res1<- res[!grepl("^Gm", res$Gene) & !grepl("Rik$", res$Gene)  & !grepl("Ptgds", res$Gene) & !grepl("^Ac", res$Gene)& !grepl("^Al", res$Gene) & !grepl("^Ap", res$Gene),]
res1[1:10, ]
# Set up the plotting window
install.packages("fields")  # Install the package if you haven't already
library(fields) 

# Plots
plot_filename <- paste0(plot_directory, "/UMAP_sub_40PC_RES0.4_trajectory_top15_trajectory_de_startvsend.pdf")

# Open the PDF device to save the plot
pdf(plot_filename, width=15, height=12) 

par(mfrow = c(6, 6), mar = c(.1, .1, 2, 1)) 


  plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings, 
       cex = .5, pch = 16, axes = F, xlab = "", ylab = "")
  
  lines(curves, lwd = 2, col = "black")
  
  # Plot centroids as black squares
  points(centroids2d[lc, ], col = "black", pch = 15, cex = 3, xpd = T)
  
  # Add labels to the centroids
  text(centroids2d[lc, ], labels = names(lc), cex = 1, font = 2, col = "white", xpd = T)


  
  # Extract and filter the variable names
  vars <- rownames(res1[1:20, ])
  vars <- na.omit(vars[vars != "NA"])
  




# Loop through each variable (gene) to plot expression levels
for (i in vars) {
  # Extract the gene expression vector and normalize it
  x <- drop0(rna.filtered@assays$RNA@data)[i, ]
  x <- (x - min(x)) / (max(x) - min(x))  # Normalize expression to range [0, 1]
  
  # Order cells based on normalized expression
  o <- order(x)
  
  color_palette <- colorRampPalette(c("lightgray", "#B34B50", "#5F070C", "#240705"))(200)
  
  # Plot UMAP embedding with the color-coded gene expression
  plot(rna.filtered@reductions$umap_after_harmony_3@cell.embeddings[o, ],
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
  xlim(0,5)+# Using RdPu palette
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
  xlim(0,5)+
  theme_classic()
dev.off()






