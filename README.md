<div align="center">
    <h2><strong>A molecular atlas of cell types in the mouse and human parkinsonian striatum</strong></h2>
    <img src="snRNA seq/Githhub_image.png" alt="Description of the image" width="900" />
</div>

Parkinson’s disease (PD) is defined by progressive dopaminergic neuron loss, disrupting the function of dopamine receptor-expressing neurons in the basal ganglia. Beyond this, the transcriptional and spatial changes in diverse striatal cell types during PD progression remain unclear. Our study combines spatial transcriptomics and single-nucleus RNA sequencing in mild and progressive PD mouse models and human postmortem tissue to map the cell-type-specific effects of dopamine degeneration. We found that both neurons and glial cells exhibit significant transcriptomic changes, with glial cells showing early and prominent effects. In mice, oligodendrocytes are the most affected glial population, whereas astrocytes display the greatest alterations in human tissue. A distinct subtype of striatal projection neurons, eccentric SPNs (eSPNs), along with Chst9+ direct-pathway spiny projection neurons (dSPNs), exhibit notable resilience to dopamine loss in mouse models and during early Braak pathology in human patients. Interestingly, while Chst9+ dSPNs maintain resilience across all Braak stages and mouse models, eSPNs show increased vulnerability as Braak stages progress and as dopamine depletion advances in the dopamine-specific mitochondrial deficiency model (MitoPark mice). This study highlights the dynamic shifts in cell states within the parkinsonian striatum, offering valuable insights into disease progression and neuronal resilience.This resource establishes the dynamic landscape of cell states in the parkinsonian striatum.


[Project Overview](#Project_Overview)

[System Requirements](#System_Requirements)

[Hardware Requirements](#Hardware_Requirements)

[Software Requirements](#Software_Requirements)

[Folder Structure](#Folder_Structure)

[License and Citation](#License_and_Citation)

[Acknowledgements](#Acknowledgements)

## Project Overview

Analyze and plot single-nucleus RNA sequencing (snRNA-seq) data for cell-specific transcriptomic changes in the striatum from Parkinson's Disease (PD) mouse models and human tissue.

Process and visualize spatial transcriptomics (ST) data to map regional vulnerabilities within the striatum.

Perform differential expression and pseudotime trajectory analyses for specific cell types, particularly focusing on astrocytes, oligodendrocytes, and spiny projection neurons (SPNs).

Apply computational models for trajectory inference and pseudotime analyses.

This repository was used in the study titled "A molecular atlas of cell types in the mouse and human parkinsonian striatum."
## System Requirements

The code in this repository is designed to run on systems capable of supporting:

- R (version 4.1.3 or later)
- Python (version 3.8 or later)
- Seurat (version 4.0.1 or later) for downstream data analysis.
- Required R Packages: libraries in the file library.R in the folder **snRNA seq/:**
## Hardware Requirements

The computational requirements to run the analyses are well within the capacity of modern mid-range laptops or desktops, with the following minimum specifications:

- 32 GB of RAM
- 4-core processor, 2.5+ GHz
- Storage: At least 50 GB free disk space for temporary files and intermediate results.
## Folder Structure
The repository is organized as follows:

- snRNA_seq/: Scripts for analyzing single-nucleus RNA sequencing data.
    - Conversion from anndata for cellxgene.ipynb and Dataset_fix_for_cellxgene_mouse.R: Python and R script used to format and convert the data for CellXGene resource.
    - gene_matrix_creation_for_mapmycell.ipynb: Python script used to convert the expression matrix extracted from the seurat object to map the cells to Allen tool MapMyCells (https://portal.brain-map.org/atlases-and-data/bkp/mapmycells)
    - snRNA_seq/Mouse datasets: Code to analize and to generate figures for mouse datasets. Divided in two subfolders, for analysis of the full dataset and subclustering.

       - snRNA_seq/Human datasets/Analysis_full_dataset
       - snRNA_seq/Human datasets/Sub_clustering_SPNs_mouse
          
    - snRNA_seq/Human datasets: Code to analize and to generate figures for human datasets. Divided in two subfolders, for analysis of the full dataset and subclustering.

       - snRNA_seq/Human datasets/Analysis_full_dataset
       - snRNA_seq/Human datasets/Sub_clustering_SPNs_human

- Spatial_Transcriptomics/: Code to analize and to generate figures of spatial transcriptomics datasets.

- Trajectory_Inference/: R scripts for pseudotime and trajectory analyses.


- Datasets/: Instructions for accessing raw and processed data.


##  License and Citation

If you use this code in your research, please cite the following:

A molecular atlas of cell types in the mouse and human parkinsonian striatum. (2025). Authors: Marta Graziano, Ioannis Mantas, Yuvarani Masarapu, Solène Frapard, Raquel
Garza, (..), Johan Jakobsson, Stefania Giacomello and Konstantinos Meletis

## Acknowledgements

This research was funded in whole by Aligning Science Across Parkinson’s **[ASAP-020370]** through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).
