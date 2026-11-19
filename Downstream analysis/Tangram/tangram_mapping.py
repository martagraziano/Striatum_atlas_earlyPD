"""
Tangram deconvolution of mouse Visium spatial transcriptomics (control and
mild 6-OHDA sections) with the mouse snRNA-seq SPN subclusters as reference.

Used for: Figure 4B (Tangram-based annotation of spatial clusters),
Figure S8A (cell-type enrichment z-scores across striatal spatial clusters)
and Figure S8B (spot-level Tangram probabilities in clusters 4, 5, 8, 27).

Inputs (exported from the Seurat objects):
  - st_object_complete.h5ad   Visium object, raw UMI counts in X and raw.X,
                              obs: sample_id, sample_type, seurat_clusters,
                              obsm: X_X_spatial
  - snRNA_reference.h5ad      snRNA-seq reference, raw counts,
                              obs: cluster_name
Outputs are written to OUT_DIR.

Run as a script or cell by cell (# %% cells).
"""

# %% ============================================================
# Configuration
# ================================================================
import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import tangram as tg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ST_PATH = Path("data/st_object_complete.h5ad")
SC_PATH = Path("data/snRNA_reference.h5ad")
OUT_DIR = Path("results/tangram")
OUT_H5AD = OUT_DIR / "st_object_with_tangram_results_subsn.h5ad"
OUT_DIR.mkdir(parents=True, exist_ok=True)

section_col = "sample_id"
condition_col = "sample_type"

# Spatial clusters shown in Figure S8B and used as the striatal baseline
target_clusters = ["4", "5", "8", "27"]
striatal_clusters = ["4", "5", "8", "27", "24", "31", "6", "25"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# %% ============================================================
# Load data
# ================================================================
adata_st = ad.read_h5ad(ST_PATH)
adata_sc = ad.read_h5ad(SC_PATH)

# Both objects must contain raw UMI counts (needed for density_prior="rna_count_based")
print(adata_st.X.max(), adata_st.raw.X.max())

adata_st.obsm["spatial"] = adata_st.obsm["X_X_spatial"]

# %% ============================================================
# Harmonise gene names (symbols, capitalised) and check overlap
# ================================================================
if "_index" in adata_st.var.columns:
    adata_st.var_names = adata_st.var["_index"].astype(str)
else:
    logger.info("'_index' column not found in adata_st.var - using existing var_names.")

adata_st.var_names = adata_st.var_names.str.strip().str.capitalize()
adata_st.var_names_make_unique()

adata_sc.var_names = adata_sc.var_names.str.strip().str.capitalize()
adata_sc.var_names_make_unique()

overlap = set(adata_sc.var_names).intersection(set(adata_st.var_names))
print(f"Gene overlap snRNA-seq / ST: {len(overlap)}")
if len(overlap) < 500:
    logger.warning(f"Only {len(overlap)} overlapping genes found - check gene symbol conventions.")

# %% ============================================================
# Reference preparation and marker selection
# ================================================================
min_cells = 3
counts = adata_sc.obs["cluster_name"].value_counts()
valid_types = counts[counts >= min_cells].index
adata_sc_filtered = adata_sc[adata_sc.obs["cluster_name"].isin(valid_types)].copy()
if hasattr(adata_sc_filtered.obs["cluster_name"], "cat"):
    adata_sc_filtered.obs["cluster_name"] = adata_sc_filtered.obs["cluster_name"].cat.remove_unused_categories()

# Keep raw counts for mapping; normalise a copy for marker ranking only
adata_sc_filtered.raw = adata_sc_filtered
sc.pp.normalize_total(adata_sc_filtered, target_sum=1e4)
sc.pp.log1p(adata_sc_filtered)

# Top 100 marker genes per snRNA-seq cluster (training genes)
sc.tl.rank_genes_groups(adata_sc_filtered, groupby="cluster_name", use_raw=False)
markers_df = pd.DataFrame(adata_sc_filtered.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
markers = [str(g).strip().capitalize() for g in markers]

adata_sc_map = adata_sc_filtered.raw.to_adata()

# %% ============================================================
# Tangram mapping (cluster mode) and projection of cell-type probabilities
# ================================================================
tg.pp_adatas(adata_sc_map, adata_st, genes=markers)

ad_map = tg.map_cells_to_space(
    adata_sc_map,
    adata_st,
    mode="clusters",
    cluster_label="cluster_name",
    density_prior="rna_count_based",
    num_epochs=500,
    device="cpu",
)

tg.project_cell_annotations(ad_map, adata_st, annotation="cluster_name")

# Store per-spot probabilities also in obs (saved with the h5ad)
df_pred = adata_st.obsm["tangram_ct_pred"]
for col in df_pred.columns:
    adata_st.obs[col] = df_pred[col]

# %% ============================================================
# Validation: gene projection, training vs test gene scores
# ================================================================
adata_sc_avg = tg.mapping_utils.adata_to_cluster_expression(adata_sc_map, "cluster_name")
common_clusters = ad_map.obs.index.intersection(adata_sc_avg.obs.index)
ad_map_aligned = ad_map[common_clusters].copy()
adata_sc_aligned = adata_sc_avg[common_clusters].copy()

ad_ge = tg.project_genes(adata_map=ad_map_aligned, adata_sc=adata_sc_aligned)
df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc_aligned)

train_scores = df_all_genes.loc[df_all_genes["is_training"], "score"]
test_scores = df_all_genes.loc[~df_all_genes["is_training"], "score"]

summary_rows = []
for label, scores in [("training", train_scores), ("test", test_scores)]:
    summary_rows.append({
        "gene_set": label,
        "n_genes": len(scores),
        "mean_score": scores.mean(),
        "median_score": scores.median(),
        "pct_above_0.3": 100 * (scores > 0.3).mean(),
        "pct_above_0.5": 100 * (scores > 0.5).mean(),
        "pct_above_0.7": 100 * (scores > 0.7).mean(),
    })
pd.DataFrame(summary_rows).to_csv(OUT_DIR / "train_vs_test_score_summary.csv", index=False)

# Seaborn compatibility shim: tangram.plot_auc passes x, y positionally,
# which newer seaborn versions no longer accept.
_orig_scatterplot = sns.scatterplot


def _scatterplot_compat(*args, **kwargs):
    if len(args) >= 2:
        kwargs["x"], kwargs["y"] = args[0], args[1]
        args = args[2:]
    return _orig_scatterplot(*args, **kwargs)


sns.scatterplot = _scatterplot_compat
try:
    tg.plot_training_scores(ad_map, bins=20, alpha=.5)
    plt.savefig(OUT_DIR / "tangram_training_scores.pdf", bbox_inches="tight")
    plt.close("all")

    tg.plot_auc(df_all_genes)
    plt.savefig(OUT_DIR / "tangram_auc.pdf", bbox_inches="tight")
    plt.close("all")

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.histplot(train_scores, color="tab:blue", label=f"Training (n={len(train_scores)})",
                 stat="density", kde=True, alpha=0.5, ax=ax)
    sns.histplot(test_scores, color="tab:orange", label=f"Test (n={len(test_scores)})",
                 stat="density", kde=True, alpha=0.5, ax=ax)
    ax.set_xlabel("Score (cosine similarity)")
    ax.set_title("Training vs. test gene scores")
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "tangram_train_vs_test_scores.pdf", bbox_inches="tight")
    plt.close(fig)
finally:
    sns.scatterplot = _orig_scatterplot

# %% ============================================================
# Per-section probability maps, all cell types, shared colour scale
# ================================================================
sample_ids = sorted(adata_st.obs[section_col].dropna().unique())
annotation_list = list(adata_sc_map.obs["cluster_name"].unique())
n_annot = len(annotation_list)
coord_key = "X_X_spatial" if "X_X_spatial" in adata_st.obsm else "spatial"

n_cols = 4
n_rows = int(np.ceil(n_annot / n_cols))
SHARED_VMIN = 0.0
SHARED_VMAX = 0.001
map_cmap = sns.blend_palette(["#e0e0e0", "#B34B50", "#5F070C", "#240705"], as_cmap=True)

with PdfPages(OUT_DIR / "all_sections_tangram_all_14_single_page_shared_scale_custompalette.pdf") as pdf:
    for s_id in sample_ids:
        adata_sub = adata_st[adata_st.obs[section_col] == s_id]
        cond_val = str(adata_sub.obs[condition_col].iloc[0])
        coords = adata_sub.obsm[coord_key]
        x, y = coords[:, 0], coords[:, 1]

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 12))
        axes = axes.flatten()
        scatter_obj = None
        for i, ct in enumerate(annotation_list):
            ax = axes[i]
            val = df_pred.loc[adata_sub.obs_names, ct].values
            scatter_obj = ax.scatter(x, y, c=val, cmap=map_cmap, s=3.5,
                                     vmin=SHARED_VMIN, vmax=SHARED_VMAX, edgecolor="none")
            ax.set_title(ct, fontsize=9, fontweight="bold", pad=3)
            ax.axis("off")
            ax.set_aspect("equal")
            ax.invert_yaxis()
        for j in range(n_annot, len(axes)):
            axes[j].axis("off")

        fig.suptitle(f"Sample ID: {s_id}   |   Sample Type: {cond_val}",
                     fontsize=14, fontweight="bold", y=0.97)
        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.04, right=0.90, hspace=0.18, wspace=0.18)
        cbar_ax = fig.add_axes([0.92, 0.20, 0.015, 0.60])
        cbar = fig.colorbar(scatter_obj, cax=cbar_ax)
        cbar.set_label("Tangram Prediction Probability", fontsize=10, labelpad=8)
        cbar.ax.tick_params(labelsize=8)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

# %% ============================================================
# Cluster-level enrichment: mean probability per spatial cluster,
# z-scored across striatal clusters ("Splatter" excluded)
# ================================================================
df_pred = adata_st.obsm["tangram_ct_pred"].copy()
if not isinstance(df_pred, pd.DataFrame):
    df_pred = pd.DataFrame(df_pred, index=adata_st.obs_names,
                           columns=list(adata_sc_map.obs["cluster_name"].unique()))

splatter_cols = [c for c in df_pred.columns if "splatter" in str(c).lower()]
if splatter_cols:
    print(f"Excluding cell type(s): {splatter_cols}")
    df_pred.drop(columns=splatter_cols, inplace=True)

df_pred_cluster = df_pred.copy()
df_pred_cluster["seurat_clusters"] = adata_st.obs["seurat_clusters"].astype(str).values


def zscore_cols(df):
    return df.apply(lambda col: (col - col.mean()) / col.std() if col.std() > 0 else col * 0, axis=0)


df_striatal = df_pred_cluster[df_pred_cluster["seurat_clusters"].isin(striatal_clusters)]
striatal_cluster_means = df_striatal.groupby("seurat_clusters").mean().loc[striatal_clusters]
striatal_cluster_zscores = zscore_cols(striatal_cluster_means)

# --- Tables for the four target clusters ---
target_means = striatal_cluster_means.loc[target_clusters]
target_zscores = striatal_cluster_zscores.loc[target_clusters]
target_means.to_csv(OUT_DIR / "target_clusters_mean_tangram_prob.csv")
target_zscores.to_csv(OUT_DIR / "target_clusters_enrichment_zscore.csv")

long_means = target_means.reset_index().melt(
    id_vars="seurat_clusters", var_name="cell_type", value_name="mean_probability")
long_zscores = target_zscores.reset_index().melt(
    id_vars="seurat_clusters", var_name="cell_type", value_name="enrichment_zscore")
summary_table = pd.merge(long_means, long_zscores, on=["seurat_clusters", "cell_type"])
summary_table.sort_values(by=["seurat_clusters", "enrichment_zscore"], ascending=[True, False], inplace=True)
summary_table.to_csv(OUT_DIR / "target_clusters_scores_summary_table.csv", index=False)

# --- Heatmap of z-scores across striatal clusters (Figure S8A) ---
heat_cmap = sns.blend_palette(["#2c7bb6", "#ffffff", "#5F070C"], as_cmap=True)
norm = TwoSlopeNorm(vcenter=0.0, vmin=-2.0, vmax=2.0)

plt.figure(figsize=(10, 5))
sns.heatmap(striatal_cluster_zscores, cmap=heat_cmap, norm=norm, annot=False,
            linewidths=0.5, linecolor="#e0e0e0",
            cbar_kws={"label": "Enrichment Z-Score (Striatal Baseline)"})
plt.title("Deconvolved Cell-Type Specificity Across Striatal Spatial Clusters",
          fontsize=12, fontweight="bold", pad=12)
plt.ylabel("Seurat Spatial Cluster", fontsize=11, fontweight="bold")
plt.xlabel("Mapped snRNA-seq Cell Type", fontsize=11, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()
plt.savefig(OUT_DIR / "spatial_clusters_zscore_heatmap_allstriatum.pdf", format="pdf", bbox_inches="tight")
plt.close("all")

# %% ============================================================
# Top-5 enriched cell types per target cluster: bar plots and
# per-section spot probability maps (Figure S8B)
# ================================================================
Z_AXIS_LIMIT = 3.0

top_z_dict = {}
selected_cell_types = set()
for cl in target_clusters:
    top_z = target_zscores.loc[cl].sort_values(ascending=False).head(5)
    top_z_dict[cl] = top_z
    selected_cell_types.update(top_z.index.tolist())

# Shared spot-probability scale: 98th percentile of non-zero probabilities
all_target_spots = adata_st.obs["seurat_clusters"].astype(str).isin(target_clusters)
probs_target = df_pred.loc[adata_st.obs_names[all_target_spots], list(selected_cell_types)].values
nonzero_probs = probs_target[probs_target > 0]
spot_vmax_global = float(np.percentile(nonzero_probs, 98)) if len(nonzero_probs) > 0 else 1e-4

spot_cmap = LinearSegmentedColormap.from_list(
    "custom_tangram", ["#f0f0f0", "#ff6b6b", "#c91818", "#700505", "#240000"])

# --- Bar plots ---
fig, axs = plt.subplots(len(target_clusters), 1, figsize=(9, 2.8 * len(target_clusters)), sharex=True)
for i, cl in enumerate(target_clusters):
    top5_z = top_z_dict[cl]
    vals = top5_z.values
    bar_colors = ["#5F070C" if v >= 0 else "#2c7bb6" for v in vals]
    y_pos = np.arange(len(top5_z))
    bars = axs[i].barh(y_pos, vals, color=bar_colors, edgecolor="black", linewidth=0.5, align="center")
    axs[i].set_yticks(y_pos)
    axs[i].set_yticklabels(top5_z.index, fontsize=9.5, fontweight="bold")
    axs[i].invert_yaxis()
    axs[i].set_xlim(-Z_AXIS_LIMIT, Z_AXIS_LIMIT)
    axs[i].axvline(0, color="black", linestyle="-", linewidth=0.8, alpha=0.9)
    axs[i].grid(axis="x", linestyle="--", alpha=0.35)
    axs[i].set_title(f"Spatial Cluster {cl} — Top Cell Types (Striatal Baseline)",
                     fontsize=10.5, fontweight="bold", pad=6)
    axs[i].tick_params(labelsize=9)
    offset = Z_AXIS_LIMIT * 0.02
    for bar, val in zip(bars, vals):
        axs[i].text(val + offset if val >= 0 else val - offset,
                    bar.get_y() + bar.get_height() / 2, f"{val:.2f}",
                    va="center", ha="left" if val >= 0 else "right",
                    fontsize=8, fontweight="bold", color="#333333")
axs[-1].set_xlabel("Enrichment Z-score (vs. Striatal Baseline)", fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT_DIR / "Target_Clusters_Top5_Enriched_Barplots.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- One page per section: z-score bars, cluster domain, top-1/top-2 spot probabilities ---
with PdfPages(OUT_DIR / "Zscore_Summary_All_Sections_Striatal_Baseline_SharedScale.pdf") as pdf:
    for target_section in sample_ids:
        adata_sample = adata_st[adata_st.obs[section_col] == target_section].copy()
        stype = str(adata_sample.obs[condition_col].iloc[0])
        coords = adata_sample.obsm["spatial"] if "spatial" in adata_sample.obsm else adata_sample.obsm["X_X_spatial"]
        x, y = coords[:, 0], coords[:, 1]

        fig, axs = plt.subplots(len(target_clusters), 4, figsize=(22, 18),
                                gridspec_kw={"width_ratios": [1.4, 1.0, 1.0, 1.0]})
        for row_idx, cl in enumerate(target_clusters):
            top5_z = top_z_dict[cl]
            st_names = top5_z.index.tolist()
            st_zscores = top5_z.values
            cluster_mask = (adata_sample.obs["seurat_clusters"].astype(str) == cl).values

            # Col 1: z-score bars
            y_pos = np.arange(len(st_names))
            bar_colors = ["#5F070C" if v >= 0 else "#2c7bb6" for v in st_zscores]
            axs[row_idx, 0].barh(y_pos, st_zscores, color=bar_colors, edgecolor="black",
                                 linewidth=0.4, align="center")
            axs[row_idx, 0].set_yticks(y_pos)
            axs[row_idx, 0].set_yticklabels([s[:26] + "..." if len(s) > 26 else s for s in st_names],
                                            fontsize=8.5, fontweight="bold")
            axs[row_idx, 0].invert_yaxis()
            axs[row_idx, 0].set_xlim(-Z_AXIS_LIMIT, Z_AXIS_LIMIT)
            axs[row_idx, 0].axvline(0, color="black", linestyle="-", linewidth=0.8, alpha=0.8)
            axs[row_idx, 0].set_xlabel("Enrichment Z-score (vs. striatal clusters)", fontsize=8.5)
            axs[row_idx, 0].set_title(f"Cluster {cl}: Top Enriched Cell Types\n[Striatal Cohort Z-Score]",
                                      fontsize=9.5, fontweight="bold")
            axs[row_idx, 0].grid(axis="x", linestyle="--", alpha=0.4)

            # Col 2: spots of the spatial cluster
            axs[row_idx, 1].scatter(x[~cluster_mask], y[~cluster_mask], c="#e8e8e8", s=6, edgecolors="none")
            axs[row_idx, 1].scatter(x[cluster_mask], y[cluster_mask], c="#5F070C", s=12,
                                    edgecolors="black", linewidths=0.15)
            axs[row_idx, 1].set_title(f"Cluster {cl} Domain\n[n={cluster_mask.sum()} spots on slice]",
                                      fontsize=9.5, fontweight="bold")
            axs[row_idx, 1].axis("off")
            axs[row_idx, 1].set_aspect("equal")
            axs[row_idx, 1].invert_yaxis()

            # Cols 3-4: spot-level probabilities of the top-1 and top-2 cell types
            for idx, col_idx in enumerate([2, 3]):
                st = st_names[idx]
                probs = df_pred.loc[adata_sample.obs_names, st].values
                axs[row_idx, col_idx].scatter(x[~cluster_mask], y[~cluster_mask], c="#e8e8e8", s=6,
                                              edgecolors="none")
                sc_p = axs[row_idx, col_idx].scatter(x[cluster_mask], y[cluster_mask], c=probs[cluster_mask],
                                                     cmap=spot_cmap, vmin=0.0, vmax=spot_vmax_global, s=12,
                                                     edgecolors="black", linewidths=0.15)
                short_title = st[:20] + "..." if len(st) > 20 else st
                axs[row_idx, col_idx].set_title(f"#{idx+1} Enriched: {short_title}\n[Spot-Level Tangram Probability]",
                                                fontsize=8.5, fontweight="bold")
                axs[row_idx, col_idx].axis("off")
                axs[row_idx, col_idx].set_aspect("equal")
                axs[row_idx, col_idx].invert_yaxis()
                cbar = fig.colorbar(sc_p, ax=axs[row_idx, col_idx], fraction=0.035, pad=0.02)
                cbar.set_label(f"Spot Prob (0 to {spot_vmax_global:.1e})", fontsize=6.5)
                cbar.ax.tick_params(labelsize=6.5)
                cbar.formatter.set_powerlimits((0, 0))

        fig.suptitle(
            f"Spatial Cluster Deconvolution Summary — Section: {target_section} | Condition: {stype}\n"
            f"[Col 1: Striatal Z-Score Enrichment (±{Z_AXIS_LIMIT:.1f}) | Col 2: Cluster Domain | "
            f"Cols 3 & 4: Spot Probabilities (0–{spot_vmax_global:.1e})]",
            fontsize=13, fontweight="bold", y=0.99)
        plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.97])
        pdf.savefig(fig, bbox_inches="tight", dpi=300)
        plt.close(fig)

# %% ============================================================
# Full striatal table (mean probability and z-score for every cell type)
# and per-cluster bar plots. Baseline here includes cluster 26 (as in the
# original analysis), so z-scores differ slightly from the heatmap above.
# ================================================================
striatal_clusters_full = ["4", "5", "8", "27", "24", "31", "6", "25", "26"]

df_striatal_full = df_pred_cluster[df_pred_cluster["seurat_clusters"].isin(striatal_clusters_full)]
means_full = df_striatal_full.groupby("seurat_clusters").mean().loc[striatal_clusters_full]
zscore_full = zscore_cols(means_full)

long_means = means_full.reset_index().melt(
    id_vars="seurat_clusters", var_name="cell_type", value_name="mean_probability")
long_zscores = zscore_full.reset_index().melt(
    id_vars="seurat_clusters", var_name="cell_type", value_name="enrichment_zscore")
merged_summary = pd.merge(long_means, long_zscores, on=["seurat_clusters", "cell_type"])
merged_summary.sort_values(by=["seurat_clusters", "enrichment_zscore"], ascending=[True, False], inplace=True)
merged_summary.to_csv(OUT_DIR / "striatal_clusters_mean_prob_and_zscore.csv", index=False)

global_prob_max = means_full.values.max() * 1.10
with PdfPages(OUT_DIR / "striatal_clusters_vertical_barplots_prob_and_zscore.pdf") as pdf:
    for cl in striatal_clusters_full:
        sorted_cell_types = zscore_full.loc[cl].sort_values(ascending=False).index.tolist()
        z_vals = zscore_full.loc[cl, sorted_cell_types].values
        prob_vals = means_full.loc[cl, sorted_cell_types].values
        z_colors = ["#5F070C" if z >= 0 else "#2c7bb6" for z in z_vals]
        xpos = np.arange(len(sorted_cell_types))

        fig, axs = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
        bars_z = axs[0].bar(xpos, z_vals, color=z_colors, edgecolor="black", linewidth=0.5, width=0.7)
        axs[0].set_ylim(-3.0, 3.0)
        axs[0].axhline(0, color="black", linestyle="-", linewidth=0.8, alpha=0.9)
        axs[0].grid(axis="y", linestyle="--", alpha=0.35)
        axs[0].set_ylabel("Enrichment Z-Score", fontsize=10, fontweight="bold")
        axs[0].set_title(f"Spatial Cluster {cl} — Enrichment Z-Score (Heatmap Metric)",
                         fontsize=12, fontweight="bold", pad=8)
        axs[0].tick_params(labelsize=9)
        for bar, val in zip(bars_z, z_vals):
            axs[0].text(bar.get_x() + bar.get_width() / 2, val + (0.12 if val >= 0 else -0.22), f"{val:.2f}",
                        ha="center", va="bottom" if val >= 0 else "top",
                        fontsize=7.5, fontweight="bold", color="#222222")

        axs[1].bar(xpos, prob_vals, color="#5F070C", edgecolor="black", linewidth=0.5, width=0.7)
        axs[1].set_ylim(0, global_prob_max)
        axs[1].grid(axis="y", linestyle="--", alpha=0.35)
        axs[1].set_ylabel("Mean Tangram Probability", fontsize=10, fontweight="bold")
        axs[1].set_title(f"Spatial Cluster {cl} — Mean Prediction Probability",
                         fontsize=12, fontweight="bold", pad=8)
        axs[1].yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter("%.1e"))
        axs[1].tick_params(labelsize=9)
        axs[1].set_xticks(xpos)
        axs[1].set_xticklabels(sorted_cell_types, rotation=45, ha="right", fontsize=9.5, fontweight="bold")
        axs[1].set_xlabel("Mapped snRNA-seq Cluster", fontsize=10.5, fontweight="bold")

        plt.suptitle(f"Striatal Cluster {cl} Deconvolution Metrics", fontsize=14, fontweight="bold", y=0.99)
        plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.97])
        pdf.savefig(fig, bbox_inches="tight", dpi=300)
        plt.close(fig)

# %% ============================================================
# Save the ST object with Tangram results
# ================================================================
# Drop the '_index' column that blocks h5ad writing
for df_attr in ["obs", "var"]:
    df = getattr(adata_st, df_attr)
    if "_index" in df.columns:
        df.drop(columns=["_index"], inplace=True)

if adata_st.raw is not None:
    raw_var = adata_st.raw.var.copy()
    if "_index" in raw_var.columns:
        raw_var.drop(columns=["_index"], inplace=True)
    adata_st.raw = ad.AnnData(X=adata_st.raw.X, var=raw_var, obs=adata_st.obs)

adata_st.write_h5ad(OUT_H5AD)
