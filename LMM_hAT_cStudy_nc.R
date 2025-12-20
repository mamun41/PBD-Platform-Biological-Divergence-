# -------- 1. Load data & pre-processing ------------------------------------
# Bulk vs Pseudo-bulk (PB) RNA-seq cross-platform divergence identification  
# Author: Md Mamunur Rashid
# Date: 2025-12-20

# Inputs :
#   - Bulk RNA-seq raw counts
#   - snRNA-seq Seurat objects


# Outputs:
#   - Normalized bulk & PB matrices
#   - Matched metadata
#   - QC-filtered scRNA objects
#   - PBD (Platform Biological Divergence) 
###########################################################################

rm(list = ls())


# --------------------------------------------------------------------------
# 1.0 Environment & paths
# --------------------------------------------------------------------------

## Set working directory to script location (RStudio only)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  dp <- dirname(rstudioapi::getSourceEditorContext()$path)
  if (dir.exists(dp)) setwd(dp)
}

pdir <- file.path(dp, "hAT_fig")
dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

## Helper functions
source("pbulk_nc_utils.R")

# --------------------------------------------------------------------------
# 1.1 Load BULK/PB RNA-seq data
# --------------------------------------------------------------------------

bulk_data <- readRDS("hATdata/bulk_count_data.rds")
pb_data <- readRDS("hATdata/pb_count_data.rds")


# --------------------------------------------------------------------------
# 1.2 Metadata
# --------------------------------------------------------------------------

metaInfo <- readRDS("hATdata/metaInfo.rds")


# --------------------------------------------------------------------------
# 1.3 Load snRNA-seq datasets
# --------------------------------------------------------------------------

single.cell.hSAT <- readRDS("hATdata/sc_hSAT_seurat.rds")
single.cell.hVAT <- readRDS("hATdata/sc_hVAT_seurat.rds")

cat("Cells | hSAT:", ncol(single.cell.hSAT), "| hVAT:", ncol(single.cell.hVAT), "\n")


# --------------------------------------------------------------------------
# 1.4 Normalization & gene filtering
# --------------------------------------------------------------------------

common_genes <- intersect(rownames(bulk_data), rownames(pb_data)); cat("nCommon genes: ", length(common_genes))
norm_res <- filter_and_normalize(
  count = cbind(bulk_data[common_genes, ], pb_data[common_genes, ]),
  meta  = metaInfo)

# log2CPM data 
logCPM_all <- norm_res$logcpm
dim(logCPM_all)



# ============================================================================
# 2. Single-cell QC & visualization
# ============================================================================

fig_1a <- DimPlot(single.cell.hSAT, group.by = "Cell_type") + 
              DimPlot(single.cell.hVAT, group.by = "Cell_type")

ggsave("fig_1a_sc_celltype_vis.png", fig_1a, width = 14, height = 6, dpi = 300, path = pdir)




# ============================================================================
# 3. Bulk vs Pseudobulk global evaluation (Figure 1)
# ============================================================================

# ---------------------------------------------------------------------------
# 3.1 Global PCA structure
# ---------------------------------------------------------------------------

pca_all <- run_pca(logCPM_all, do.center = TRUE, do.scale = FALSE)

fig_1b <- plot_pca(pca_all, metaInfo, title = " ")
ggsave("fig_1b_PCA_orig.png", fig_1b, width = 5, height = 4, path = pdir)


# ---------------------------------------------------------------------------
# 3.2 Mean expression agreement (scatter)
# ---------------------------------------------------------------------------

res_scatter <- plot_bulk_pb_scatter(
  matX = logCPM_all, meta = metaInfo, method_cor = "pearson", title = " ")

ggsave("fig_1c_Bulk_PB_scatter.png", res_scatter$plot, width = 6, height = 4, path = pdir)


# ---------------------------------------------------------------------------
# 3.3 Marker gene concordance
# ---------------------------------------------------------------------------

marker_genes <- c("WDPCP", "NEGR1", "CAV2", "ALCAM", "FABP4")

fig_1c <- plot_marker_slopes(
  logCPM_all, metaInfo, marker_genes,
  save_path = file.path(pdir, "fig_1c_marker_genes_scatter.png"),
  width = 6, height = 4
  )


# ---------------------------------------------------------------------------
# 3.4 Correlation-based similarity
# ---------------------------------------------------------------------------

cor_mat <- corr_matrix(logCPM_all, method = "pearson")

cor_bar <- plot_corr_summary(cor_mat, metaInfo); cor_bar
ggsave("fig_1d_cor_barplot.png", cor_bar$plot, width = 5, height = 5, path = pdir)

fig_1d_violin <- plot_corr_distribution(cor_mat, metaInfo, "pearson", title = " ")
ggsave("fig_1d_cor_violin.png",fig_1d_violin, width = 5, height = 5, path = pdir)


# ============================================================================
# 4. Feature divergence & enrichment (Figure 1e, S1)
# ============================================================================

# ---------------------------------------------------------------------------
# 4.1 Sparse PCA for divergent gene discovery
# ---------------------------------------------------------------------------

spca_all <- run_spca(logCPM_all, k = 5, alpha.val = 1e-4, beta.val = 1e-3)

g_pc1 <- rownames(logCPM_all)[spca_all$loadings[,1] != 0]


# ---------------------------------------------------------------------------
# 4.2 GSVA enrichment (PC1 genes)
# ---------------------------------------------------------------------------

res_gsva_PC1 <- run_gsva_enrichment(
  expr_mat = logCPM_all,
  genes_of_interest = g_pc1,
  metaInfo = metaInfo,
  title_prefix = " "
)

ggsave("fig_1h_GO_pc1.png", res_gsva_PC1$barplot + NoGrid(), width = 9, height = 7, path = pdir)


# ---------------------------------------------------------------------------
# 4.3 Global distribution diagnostics 
# ---------------------------------------------------------------------------

fig_1e <- plot_expression_density(logCPM_all, metaInfo)
ggsave("fig_1e_KDE_sim.png", fig_1e, width = 6, height = 4, path = pdir)



# ============================================================================
# 5. Platform Bias Detection (PBD) & Validation
# ============================================================================
vp <- readRDS("manuscript_hAT_ncomm_fig/varPar_res.rds")

# do not run : can delay execution
# param <- SnowParam(3, "SOCK", progressbar = TRUE)
# 
# metadata <- metaInfo
# metadata$Platform <- ifelse(metadata$platform == "bulk", 1, 0)
# 
# form <- ~ Platform + (1 | Sample) + (1 | Tissue)
# vp <- fitExtractVarPartModel(logCPM_all, form, metadata, BPPARAM = param)



# ---------------------------------------------------------------------------
# 5.1 PBD Score Definition & Thresholding
# ---------------------------------------------------------------------------

# PBD = Platform component - mean(random effect components)
bio_factors <- c("Sample", "Tissue")
vp$PBD <- vp$Platform - rowMeans(vp[, bio_factors, drop = FALSE])

# cutoff selection from bimodal curve
pbd_result <- find_peak_threshold(vp$PBD, bandwidth = 0.1, select = "last")
fig_pbd <- visualize_peak_detection(vp$PBD, pbd_result, xname = "PBD score"); fig_pbd

# ggsave("fig_2a_PBD_cutoff.png", fig_pbd + NoGrid(), width = 6, height = 4, path = pdir)


# ---------------------------------------------------------------------------
# 5.2 Platform-sensitive Gene Sets
# ---------------------------------------------------------------------------

pbd_gene <- rownames(vp)[vp$PBD > pbd_result$threshold_high]
pbd_gene2 <- rownames(vp)[vp$PBD > pbd_result$threshold_low]


# ---------------------------------------------------------------------------
# 5.3 Global Structure Validation (PCA, Correlation, Density) for PBD genes
# ---------------------------------------------------------------------------

# PCA (PBD highly batch sensitive)
plot_pca(run_pca(logCPM_all[pbd_gene, ]), metaInfo) 

# Correlation 
cor_low  <- plot_corr_summary(
  corr_matrix(logCPM_all[!rownames(logCPM_all) %in% pbd_gene2, ]),
  metaInfo
  ); cor_low$results

ggsave("fig_1g_cor_lowPBD.png", cor_low$plot, width = 5, height = 5, path = pdir)




# ---------------------------------------------------------------------------
# 5.4 Integrated DE + sPCA + PBD Volcano
# ---------------------------------------------------------------------------

fig_2a <- plot_volcano_multi_method(
  exMAT = logCPM_all,
  meta = metaInfo,
  spca_all = spca_all,
  vp_res = vp,
  vp_col = "PBD",
  pbd_cutoff = "threshold",
  top_n = 20,
  bandwidth = .1
)

ggsave("fig_1f_volcano.png", fig_2a$volcano_plot, width = 8, height =6, path = pdir)
res_df <- fig_2a$summary_df


genes_limma <- res_df$Gene[res_df$limma_flag == TRUE]
genes_sPCA  <- res_df$Gene[res_df$sPCA_flag == TRUE]
genes_PBD  <- res_df$Gene[res_df$PBD_flag == "High-PBD"]

all_intersect_gene <- Reduce(intersect, list(genes_limma, genes_sPCA, genes_PBD))


# ---------------------------------------------------------------------------
# 5.5 Prepare VP table and annotations
# ---------------------------------------------------------------------------

vc_df <- as.data.frame(vp)
vc_df$Gene <- rownames(vp)

# Join PBD / PBD annotations
df_vp <- left_join(vc_df, res_df, by = "Gene")
rownames(df_vp) <- df_vp$Gene
df_vp <- df_vp[, c(bio_factors, "Platform", "PBD_flag")]
colnames(df_vp)[colnames(df_vp) == "PBD_flag"] <- "Condition"


# ---------------------------------------------------------------------------
# 5.2 VP: High-PBD vs Low-PBD & Com-DE vs Non-Com-DE
# ---------------------------------------------------------------------------

df_vp_com <- df_vp
df_vp_com$Condition <- "Non.com-DE"
df_vp_com$Condition[rownames(df_vp_com) %in% all_intersect_gene] <- "Com-DE"

vp_long <- reshape2::melt(df_vp, id.vars = "Condition")
vp_long_com <- reshape2::melt(df_vp_com, id.vars = "Condition")
vp_long2 <- rbind(vp_long, vp_long_com)

vp_long2$Condition <- factor(
  vp_long2$Condition,
  levels = c("Low-PBD", "Non.com-DE", "High-PBD", "Com-DE")
)

plot_vp_box <- function(df_long) {
  ggplot(df_long, aes(x = variable, y = value, fill = Condition)) +
    geom_boxplot(width = 0.5, outlier.size = 0.5) +
    theme_classic() +
    labs(x = " ", y = "Proportion of Variance") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

fig_2f_vp.pbd <- plot_vp_box(vp_long)
ggsave("fig_2f_VariancePartition_PBD.png", fig_2f_vp.pbd, width = 6, height = 4, path = pdir)

fig_2f_vp.pbd_comDE <- plot_vp_box(vp_long2)
ggsave("fig_2f_VariancePartition_PBD_ComDE.png", fig_2f_vp.pbd_comDE, width = 6, height = 4, path = pdir)



# ============================================================================
# 6. PCA on gene subsets
# ============================================================================

# High / Low PBD
pca_h <- run_pca(logCPM_all[res_df$Gene %in% pbd_gene, ])
fig_2f_pbd_h <- plot_pca(pca_h, metaInfo, title = "High-PBD genes")

pca_l <- run_pca(logCPM_all[!(res_df$Gene %in% pbd_gene), ], do.scale = TRUE)
fig_2f_pbd_l <- plot_pca(pca_l, metaInfo, title = "Low-PBD genes", leg.pos = "none")

ggsave("fig_2f_PBD_stableProxy_vs_divergent_PCA_comparison.png",
       (fig_2f_pbd_h + fig_2f_pbd_l), width = 8, height = 6, dpi = 300, path = pdir)



# ============================================================================
# 7. Single-cell marker validation 
# ============================================================================

# ---------------------------
# hSAT
# ---------------------------
top5_genes <- single.cell.hSAT@misc$DE_sig$hSAT %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(avg_log2FC, n = 5) %>%
  dplyr::ungroup()


DotPlot(
  single.cell.hSAT,
  features = unique(top5_genes$gene),
  group.by = "Cell_type",
  cols = c("lightblue", "darkred"),
  dot.scale = 8
) + ggtitle("hSAT") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 18)) +
  labs(x = " ", y = " ")



# ---------------------------
# hVAT
# ---------------------------
top5_genes <- single.cell.hVAT@misc$DE_sig$hVAT %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(avg_log2FC, n = 5) %>%
  dplyr::ungroup()


DotPlot(
  single.cell.hSAT,
  features = unique(top5_genes$gene),
  group.by = "Cell_type",
  cols = c("lightblue", "darkred"),
  dot.scale = 8
) + ggtitle("hSAT") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 18)) +
  labs(x = " ", y = " ")



# ============================================================================
# 8. Enrichment analysis for high PBD gene sets
# ============================================================================

res_gsva_pbd <- run_gsva_enrichment2(
  expr_mat = logCPM_all,
  genes_of_interest = pbd_gene2,
  metaInfo = metaInfo,
  heatmap_option = "ComplexHeatmap"
)

res_gsva_pbd$barplot
save_plot(res_gsva_pbd$heatmap, "fig_1h_GO-Enrichment_PBD_BarPlot.png", path = pdir, width = 12, height = 7)


#---------------------- End ---------------------------------




