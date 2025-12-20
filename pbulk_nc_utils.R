########################################################################


# -------- 0. Setup & packages ----------------------------------------------
# # Auto-install missing packages (CRAN + Bioconductor) and load quietly
# 
# ## CRAN packages
# cran_pkgs <- c(
#   "ggplot2","patchwork","scales","RColorBrewer","reshape2","philentropy",
#   "cluster","mclust","Matrix","tidyr","purrr","tidyverse","pheatmap",
#   "sparsepca","entropy","e1071","matrixStats","gridExtra","ggrepel","grid",
#   "randomcoloR","stringr","cowplot","eulerr","UpSetR",
#   "parallel","doParallel","data.table"
# )
# 
# ## Bioconductor packages
# bioc_pkgs <- c(
#   "edgeR","limma","Seurat","clusterProfiler","org.Hs.eg.db","scater",
#   "sva","harmony","ComplexHeatmap","circlize",
#   "GSVA","GSEABase","msigdbr","variancePartition","Biobase"
# )
# 
# ## Install missing CRAN packages
# cran_missing <- setdiff(cran_pkgs, rownames(installed.packages()))
# if (length(cran_missing) > 0) {
#   install.packages(cran_missing, repos = "https://cloud.r-project.org")
# }
# 
# ## Install missing Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager", repos = "https://cloud.r-project.org")
# }
# bioc_missing <- setdiff(bioc_pkgs, rownames(installed.packages()))
# if (length(bioc_missing) > 0) {
#   BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
# }
# 
# ## Load all packages quietly
# suppressPackageStartupMessages({
#   invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))
# })



# -------- 0. Setup & packages ------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(reshape2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(scater)
  library(philentropy) # JSD
  library(cluster)     # silhouette
  library(mclust)      # ARI
  library(Matrix)
  library(tidyr)
  library(purrr)
  library(tidyverse)
  library(pheatmap)
  library(sva)
  library(sparsepca)
  library(entropy)
  library(e1071)
  library(harmony)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(gridExtra)
  library(GSVA);
  library(GSEABase)
  library(msigdbr); 
  library(dplyr)
  library(ggrepel)
  library(grid)
  library(randomcoloR)
  library(stringr)
  library(cowplot)
  library(eulerr)
  library(UpSetR)
  library(variancePartition)
  library(parallel)
  library(doParallel)
  library(data.table)
  library(Biobase)
  
})


# -------- 1. MAIN ANALYSIS FUNCTIONS ----------------------------------------

# ======== Data Processing Functions ========

# Normalization: TMM -> CPM -> log2 
filter_and_normalize <- function(counts, meta) {
  dge <- edgeR::DGEList(counts = counts, group = meta$platform)
  keep <- filterByExpr(dge, group = meta$platform)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dgeCPM <- edgeR::cpm(dge, log = FALSE)
  logCPM <- log2(dgeCPM + 1)
  
  return(list(logcpm = logCPM, cpm = dgeCPM))
}



# ======== Dimensionality Reduction Functions ========

# PCA wrapper that removes zero-variance genes
run_pca <- function(mat, do.center = TRUE, do.scale = FALSE) {
  # mat: genes x samples
  nz <- apply(mat, 1, var, na.rm = TRUE) > 0
  if (sum(!nz) > 0) message(sprintf("Removing %d zero-variance genes", sum(!nz)))
  
  pr <- prcomp(t(mat[nz, , drop = FALSE]), center = do.center, scale. = do.scale)
  # add "importance" slot like in summary(prcomp)
  sdev <- pr$sdev
  var  <- sdev^2
  prop <- var / sum(var)
  cum  <- cumsum(prop)
  pr$importance <- rbind(
    "Std. Dev."              = sdev,
    "Proportion of Variance" = prop,
    "Cumulative Proportion"  = cum
  )
  
  return(pr)
}


# ========
run_spca <- function(mat, k = 10, do.center = TRUE, do.scale = FALSE, alpha.val = 1e-4, beta.val = 1e-4) {
  # mat: genes x samples
  nz <- apply(mat, 1, var, na.rm = TRUE) > 0
  if (sum(!nz) > 0) message(sprintf("Removing %d zero-variance genes", sum(!nz)))
  
  X <- t(mat[nz, , drop = FALSE])  # samples x genes
  
  if (ncol(mat) < 100) {
    pr <- sparsepca::spca(X, k = k, center = do.center, scale = do.scale, alpha = alpha.val, beta = beta.val)
  } else {
    pr <- sparsepca::rspca(X, k = k, center = do.center, scale = do.scale, alpha = alpha.val, beta = beta.val)
  }
  
  # Align with prcomp-like output
  pr$x <- pr$scores
  pr$rotation <- pr$loadings
  rownames(pr$x) <- rownames(X)
  rownames(pr$rotation) <- colnames(X)
  
  # Compute variance explained
  vars <- pr$sdev^2
  pr$importance <- rbind(
    "Standard deviation" = pr$sdev,
    "Proportion of Variance" = vars / sum(vars),
    "Cumulative Proportion" = cumsum(vars) / sum(vars)
  )
  
  class(pr) <- c("spca", class(pr))
  return(pr)
}



# ======== Correlation Analysis Functions ========

# correlation matrix
corr_matrix <- function(expr_matrix, method = "pearson") {
  cor(expr_matrix, method = method)
}

# Analyze correlation patterns
corr_patterns <- function(corr_matrix, sample_info) {
  results <- list()
  
  bulk_mask <- sample_info$platform == "bulk"
  pb_mask   <- sample_info$platform == "pb"
  
  # Within-bulk
  bulk_corr <- corr_matrix[bulk_mask, bulk_mask]
  results$within_bulk <- mean(bulk_corr[upper.tri(bulk_corr)])
  
  # Within-PB
  pb_corr <- corr_matrix[pb_mask, pb_mask]
  results$within_pb <- mean(pb_corr[upper.tri(pb_corr)])
  
  # Between-groups
  between_corr <- corr_matrix[bulk_mask, pb_mask]
  results$between_groups <- mean(between_corr)
  
  return(results)
}



# ======== JSD (Distribution) Analysis Functions ========

calculate_distribution_jsd <- function(mat, meta, nbins = 15) {
  genes <- rownames(mat)
  
  temp <- sapply(1:length(genes), function(g) {
    x <- mat[g, meta$platform == "bulk", drop = FALSE]
    y <- mat[g, meta$platform == "pb", drop = FALSE]
    # Equal-quantile binning for robust distribution comparison
    all_vals <- c(x, y)
    probs <- seq(0, 1, length.out = nbins + 1)
    breaks <- quantile(all_vals, probs = probs, na.rm = TRUE)
    
    # Handle identical distributions
    if (diff(range(breaks)) < 1e-8) return(0)
    
    # Bin values and calculate probabilities
    counts_x <- hist(x, breaks = breaks, plot = FALSE)$counts
    counts_y <- hist(y, breaks = breaks, plot = FALSE)$counts
    
    P <- (counts_x + 1e-8) / sum(counts_x + 1e-8)
    Q <- (counts_y + 1e-8) / sum(counts_y + 1e-8)
    
    # JSD calculation
    M <- 0.5 * (P + Q)
    jsd <- 0.5 * sum(P * log2(P / M)) + 0.5 * sum(Q * log2(Q / M))
    
  })
  
  names(temp) <- genes
  return(temp)
}



# ======== Differential Expression Analysis Functions ========

run_DE <- function(mat, meta) {
  # Create design matrix
  design <- model.matrix(~ 0 + platform , data = meta)
  colnames(design) <- make.names(colnames(design))
  
  # Combine data for limma
  combined_data <- mat
  
  # limma- pipeline
  fit <- limma::lmFit(combined_data, design)
  contrast <- makeContrasts(PBvsBulk = platformpb - platformbulk, levels=design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  # Extract results 
  tt <- limma::topTable(fit2, number=nrow(logCPM_all), adjust="BH", sort.by = "none")
  tt$gene <- rownames(tt)
  tt$DE_flag <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1
  
  return(tt)
}


# 
identify_divergent_genes <- function(mat, meta, fdr_threshold = 0.05, nbins = 10) {
  
  # Calculate divergence measures
  jsd_results <- calculate_divergence_measures(mat = mat,meta = meta, nbins = nbins)
  jsd_results <- calculate_distribution_jsd(mat, meta, nbins)
  # Perform DE analysis
  de_results <- run_DE(mat = mat, meta = meta)
  de_results$pSig <- de_results$adj.P.Val <= fdr_threshold
  de_results$fcSig <- abs(de_results$logFC) > 1
  
  
  genes <- names(jsd_results)
  # Create results data frame
  results_df <- data.frame(gene = genes, dist_jsd = jsd_results, stringsAsFactors = FALSE)
  results_df <- left_join(results_df, de_results, by = "gene")
  
  # Integrate results
  return(list(dist_scores = jsd_results, de_results = de_results, jsd_res = results_df))
}






# -------- 2. VISUALIZATION FUNCTIONS ----------------------------------------

# ======== PCA Visualization ========

plot_pca <- function(pr, meta_df, title = "PCA", color = "group", shape = "platform", leg.pos = "top") {
  
  df <- data.frame(PC1 = pr$x[, 1], PC2 = pr$x[, 2], sample_id = rownames(pr$x))
  df <- dplyr::left_join(df, meta_df, by = "sample_id")
  
  pct1 <- 100 * pr$importance["Proportion of Variance", 1]
  pct2 <- 100 * pr$importance["Proportion of Variance", 2]
  
  p <- ggplot(df, aes_string("PC1", "PC2", color = color, shape = shape)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_minimal(base_size = 12) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", pct1),
      y = sprintf("PC2 (%.1f%%)", pct2)) +
    theme(
      axis.text.x   = element_text(size = 14),   # larger X-axis text
      axis.text.y   = element_text(size = 14),   # larger Y-axis text
      axis.title.x  = element_text(size = 16),   # larger X-axis title
      axis.title.y  = element_text(size = 16),   # larger Y-axis title
      plot.title    = element_text(size = 16), # larger centered title
      panel.grid = element_blank(),
      axis.line     = element_line(color = "black", linewidth = 0.6), # add axis lines
      panel.border  = element_blank()                   # no border box
    )
  
  # Handle legend position
  if (tolower(leg.pos) == "none") {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = leg.pos)
  }
  
  return(p)
}



# ======== Correlation Visualization ========

plot_corr <- function(cor_matrix, palette = c("firebrick3", "white", "royalblue3"),
                      show_numbers = FALSE, title = NULL, v_methods = "shade", 
                      return = c("both", "matrix", "plot", "none")) {
  # Load required package
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("Package 'corrplot' is required. Please install it with install.packages('corrplot').")
  }
  
  return  <- match.arg(return)
  # Correlation range (exclude diagonal)
  cor_range <- range(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)   
  
  # Define color palette
  col_fun <- colorRampPalette(palette)
  # Plot correlation matrix
  par(mfrow = c(1,1), mar = c(0, 0, 2, 0))
  plot_result <- corrplot::corrplot(cor_matrix, method = v_methods, type = "upper",
                                    col = col_fun(30), is.corr = FALSE,       # ✅ prevents auto-scaling to [-1,1]
                                    cl.lim = cor_range,    # ✅ adaptive scale
                                    number.cex = ifelse(show_numbers, 0.6, 0),
                                    tl.cex = 0.6, mar = c(0, 0, 0, 0),
                                    title = title
                                    )
  
  # Return logic
  if (return == "matrix") {
    return(cor_matrix)
  } else if (return == "plot") {
    return(invisible(plot_result))
  } else if (return == "both") {
    return(list(cor_matrix = cor_matrix, plot = plot_result))
  } else {
    invisible(NULL)
  }
}




# ======== Correlation Summary Visualization ========

plot_corr_summary <- function(corr_matrix, sample_info) {
  
  bulk_mask <- sample_info$platform == "bulk"
  pb_mask   <- sample_info$platform == "pb"
  
  # Within-bulk
  bulk_corr <- corr_matrix[bulk_mask, bulk_mask]
  within_bulk <- mean(bulk_corr[upper.tri(bulk_corr)], na.rm = TRUE)
  
  # Within-PB
  pb_corr <- corr_matrix[pb_mask, pb_mask]
  within_pb <- mean(pb_corr[upper.tri(pb_corr)], na.rm = TRUE)
  
  # Between-groups
  between_corr <- corr_matrix[bulk_mask, pb_mask]
  between_groups <- mean(between_corr, na.rm = TRUE)
  
  # Data frame
  corr_df <- data.frame(
    category = c("Between", "Within-bulk", "Within-pb"),
    value = c(between_groups, within_bulk, within_pb)
  )
  
  # Barplot
  p <- ggplot(corr_df, aes(x = category, y = value, fill = category)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_minimal(base_size = 14) +
    labs(title = " ", x = "",y = "Correlation") +
    scale_fill_brewer(palette = "Set2") +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),  # larger X-axis text
      axis.text.y   = element_text(size = 14),                         # larger Y-axis text
      axis.title.x  = element_text(size = 16),    # face = "bold"      # larger X-axis title
      axis.title.y  = element_text(size = 16),          # larger Y-axis title
      plot.title    = element_text(size = 16, hjust = 0.5), # larger centered title
      panel.grid    = element_blank(),
      legend.position = "none"                                      
    )
  # Return both results and plot
  return(list(
    results = list(
      within_bulk = within_bulk,
      within_pb = within_pb,
      between_groups = between_groups
    ),
    plot = p
  ))
}




# ======== Correlation Distribution Visualization ========
plot_corr_distribution <- function(corr_matrix, sample_info, method_name, title="Correlation Distribution") {
  corr_melt <- reshape2::melt(corr_matrix)
  colnames(corr_melt) <- c("Sample1", "Sample2", "Correlation")
  
  corr_melt$Type1 <- sample_info$platform[match(corr_melt$Sample1, sample_info$sample_id)]
  corr_melt$Type2 <- sample_info$platform[match(corr_melt$Sample2, sample_info$sample_id)]
  
  corr_melt <- corr_melt[corr_melt$Sample1 != corr_melt$Sample2, ]
  corr_melt$Comparison <- ifelse(corr_melt$Type1 == corr_melt$Type2, 
                                 paste0("Within-", corr_melt$Type1),
                                 "Between")
  
  ggplot(corr_melt, aes(x = Comparison, y = Correlation, fill = Comparison)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_minimal() +
    labs(title =  title,
         y = "Correlation",
         x = " ") +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),  # larger X-axis text
      axis.text.y   = element_text(size = 14),                         # larger Y-axis text
      axis.title.x  = element_text(size = 16),    # face = "bold"      # larger X-axis title
      axis.title.y  = element_text(size = 16),          # larger Y-axis title
      plot.title    = element_text(size = 16, hjust = 0.5), # larger centered title
      panel.grid    = element_blank(),
      legend.position = "none"                                         # remove legend
    )
}






# ======== Correlation Distribution Comparison Functions ========

plot_corr_distribution_match <- function(expMat, meta, jsd_idx, title = "Correlation Before and After PBD Filter") {
  # --- Load required libraries ---
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Please install 'ggplot2' and 'ggpubr' packages before running this function.")
  }
  library(ggplot2)
  library(ggpubr)
  
  # --- Helper function to get upper triangle (no duplicates) ---
  get_upper_tri <- function(mat) mat[upper.tri(mat)]
  
  # --- 1. Between correlations: matched samples only (diagonal) ---
  cor.b <- corr_matrix(expMat)  # all genes
  cor.a <- corr_matrix(expMat[!(rownames(expMat) %in% jsd_idx), ])  # Low JSD genes 
  
  cors_between_before <- diag(cor.b[meta$platform=="bulk", meta$platform=="pb"])
  cors_between_after  <- diag(cor.a[meta$platform=="bulk", meta$platform=="pb"])
  
  # --- 2. Within correlations (sample-sample) ---
  cors_within_bulk_before <- get_upper_tri(cor.b[meta$platform=="bulk", meta$platform=="bulk"])
  cors_within_bulk_after  <- get_upper_tri(cor.a[meta$platform=="bulk", meta$platform=="bulk"])
  cors_within_pb_before   <- get_upper_tri(cor.b[meta$platform=="pb", meta$platform=="pb"])
  cors_within_pb_after    <- get_upper_tri(cor.a[meta$platform=="pb", meta$platform=="pb"])
  
  # --- 3. Combine all into a single data frame ---
  df_cors <- data.frame(
    Correlation = c(
      cors_between_before, cors_within_bulk_before, cors_within_pb_before,
      cors_between_after,  cors_within_bulk_after,  cors_within_pb_after
    ),
    Type = factor(
      rep(c("Between", "Within Bulk", "Within PB",
            "Between", "Within Bulk", "Within PB"),
          times = c(
            length(cors_between_before), length(cors_within_bulk_before), length(cors_within_pb_before),
            length(cors_between_after),  length(cors_within_bulk_after),  length(cors_within_pb_after)
          )),
      levels = c("Between", "Within Bulk", "Within PB")
    ),
    Condition = factor(
      rep(c(rep("All genes", 
                length(cors_between_before) + length(cors_within_bulk_before) + length(cors_within_pb_before)),
            rep("PBD filter",  
                length(cors_between_after) + length(cors_within_bulk_after) + length(cors_within_pb_after))),
          levels = c("All genes", "PBD filter"))
    )
  )
  
  # --- 4. Plot correlation comparison ---
  p <- ggplot(df_cors, aes(x = Type, y = Correlation, fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.5) +
    scale_fill_manual(values = c("#969696", "#238443")) +
    theme_classic(base_size = 12) + 
    labs(
      # title = "Bulk–PB correlation improves after JSD", 
      x = NULL, y = "Correlation"
    ) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),  # larger X-axis text
      axis.text.y   = element_text(size = 14),                         # larger Y-axis text
      axis.title.x  = element_text(size = 16),    # face = "bold"      # larger X-axis title
      axis.title.y  = element_text(size = 16),          # larger Y-axis title
      plot.title    = element_text(size = 16, hjust = 0.5), # larger centered title
      panel.grid    = element_blank(),
      legend.position = "right"                                         # remove legend
    ) +
    stat_compare_means(
      aes(group = Condition),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    )
  
  return(p)
}




# --- Compare global correlation improvement --- 
plot_corr_distribution_cross <- function(expMat, meta, jsd_idx, common_de_idx,
                                         title = "Bulk–PB Correlation Comparison") {
  library(ggplot2)
  library(ggpubr)
  
  # --- Helper ---
  get_upper_tri <- function(mat) mat[upper.tri(mat)]
  
  # --- Correlations ---
  cor.b    <- corr_matrix(expMat)
  cor.a    <- corr_matrix(expMat[!(rownames(expMat) %in% jsd_idx), ])
  cor.comm <- corr_matrix(expMat[!(rownames(expMat) %in% common_de_idx), ])
  
  bulk_idx <- meta$platform == "bulk"
  pb_idx   <- meta$platform == "pb"
  
  cors_between_before <- diag(cor.b[bulk_idx, pb_idx])
  cors_between_after  <- diag(cor.a[bulk_idx, pb_idx])
  cors_between_comm   <- diag(cor.comm[bulk_idx, pb_idx])
  
  # --- Long data frame ---
  df_cors <- data.frame(
    Correlation = c(cors_between_before, cors_between_after, cors_between_comm),
    Condition = factor(
      c(
        rep("All genes", length(cors_between_before)),
        rep("High-PBD filter", length(cors_between_after)),
        rep("Common DE filter", length(cors_between_comm))
      ),
      levels = c("All genes", "High-PBD filter", "Common DE filter")
    )
  )
  
  # --- Colors (publication-friendly) ---
  cols <- c("All genes" = "#969696", "High-PBD filter" = "#1B9E77","Common DE filter" = "#66A61E")
  
  # --- Plot ---
  p <- ggplot(df_cors, aes(x = Condition, y = Correlation, fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.65, linewidth = 0.6) +
    geom_boxplot(width = 0.14, alpha = 0.85, linewidth = 0.6) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = "Correlation") +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 17),
      legend.position = "none",
      plot.title = element_text(size = 18, hjust = 0.5)
    ) +
    stat_compare_means(
      aes(group = Condition),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 6
    )
  
  return(p)
}





# ======== Average Cross Correlation Visualization ========

plot_corr_distribution_ave <- function(p_before, p_after){
  # Correlation comparison 
  cors_before <- p_before$data 
  cors_after <- p_after$data
  cors_before$Condition <- "All genes"
  cors_after$Condition <- "PBD filter"
  
  df_cors <- rbind(cors_before, cors_after)
  p <- ggplot(df_cors, aes(x = Comparison, y = Correlation, fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.6, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), alpha = 0.5) +
    scale_fill_manual(values = c("#969696", "#238443")) +
    theme_classic(base_size = 12) + 
    labs(
      # title = "Bulk–PB correlation improves after JSD", 
      x = NULL, y = "Correlation"
    ) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),  # larger X-axis text
      axis.text.y   = element_text(size = 14),                         # larger Y-axis text
      axis.title.x  = element_text(size = 16),    # face = "bold"      # larger X-axis title
      axis.title.y  = element_text(size = 16),          # larger Y-axis title
      plot.title    = element_text(size = 16, hjust = 0.5), # larger centered title
      panel.grid    = element_blank(),
      legend.position = "right"                                         # remove legend
    ) +
    stat_compare_means(
      aes(group = Condition),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    )
  return(p)
}



# ======== Correlation Barplot Visualization ========
plot_corr_barplot <- function(corr_matrix, sample_info, method_name, title = "Correlation Distribution") {
  
  corr_melt <- melt(corr_matrix)
  colnames(corr_melt) <- c("Sample1", "Sample2", "Correlation")
  
  corr_melt$Type1 <- sample_info$platform[match(corr_melt$Sample1, sample_info$sample_id)]
  corr_melt$Type2 <- sample_info$platform[match(corr_melt$Sample2, sample_info$sample_id)]
  
  corr_melt <- corr_melt[corr_melt$Sample1 != corr_melt$Sample2, ]
  
  corr_melt$Comparison <- ifelse(
    corr_melt$Type1 == corr_melt$Type2,
    paste0("Within-", corr_melt$Type1),
    "Between"
  )
  
  # --- Summary statistics for barplot ---
  summary_df <- corr_melt %>%
    group_by(Comparison) %>%
    summarise(
      mean_corr = mean(Correlation, na.rm = TRUE),
      se_corr   = sd(Correlation, na.rm = TRUE) / sqrt(n())
    )
  
  # --- Barplot with error bars ---
  ggplot(summary_df, aes(x = Comparison, y = mean_corr, fill = Comparison)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean_corr - se_corr,
                      ymax = mean_corr + se_corr),
                  width = 0.2) +
    theme_minimal() +
    labs(
      title = title,
      y = "Mean Correlation",
      x = ""
    ) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 16),
      axis.title.y  = element_text(size = 16),
      plot.title    = element_text(size = 16, hjust = 0.5),
      panel.grid    = element_blank(),
      legend.position = "none"
    )
}




# ======== Expression Visualization ========

# per-gene/per-sample density
plot_expression_density <- function(mat, meta_df, title = "KDE of mean expression", 
                                    mode = c("gene", "sample"), legend_pos = "top") {   # "top", "bottom", "right", "left", "none"
  
  mode <- match.arg(mode)
  if (mode == "gene") {
    bulk_rm <- as.numeric(mat[, meta_df$platform == "bulk", drop = FALSE])
    pb_rm   <- as.numeric(mat[, meta_df$platform == "pb", drop = FALSE])
    edata <- data.frame(Bulk = bulk_rm, PB = pb_rm) %>% reshape2::melt()
    
    p <- ggplot(edata, aes(x = value, fill = variable)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = c("Bulk" = "#ef3b2c", "PB" = "#20A39E")) +
      labs(x = expression(Log[2]~CPM), y = "Density", title = title, fill = NULL) + 
      theme_minimal(base_size = 12)
    
  } else {
    edata <- data.frame(sample_id = colnames(mat), Expression = colMeans(mat), platform = meta_df$platform)
    p <- ggplot(edata, aes(x = Expression, fill = platform)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = c("bulk" = "#ef3b2c", "pb" = "#20A39E")) +
      labs(x = expression(Log[2]~CPM), y = "Density", title = title, fill = NULL) + 
      theme_minimal(base_size = 12)
  }
  
  # Legend handling
  if (tolower(legend_pos) == "none") {
    p <- p + theme(
      axis.text.x   = element_text(size = 14),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 16),
      axis.title.y  = element_text(size = 16),
      plot.title    = element_text(size = 16, hjust = 0.5),
      panel.grid = element_blank(),
      axis.line     = element_line(color = "black", linewidth = 0.6), # add axis lines
      panel.border  = element_blank(),                   # no border box
      legend.position = "none"
    )
  } else if (legend_pos %in% c("top", "bottom")) {
    p <- p + theme(
      axis.text.x   = element_text(size = 14),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 16),
      axis.title.y  = element_text(size = 16),
      plot.title    = element_text(size = 16, hjust = 0.5),
      panel.grid = element_blank(),
      axis.line     = element_line(color = "black", linewidth = 0.6), # add axis lines
      panel.border  = element_blank(),                   # no border box
      legend.position = legend_pos,
      legend.direction = "horizontal"
    )
  } else {
    p <- p + theme(
      axis.text.x   = element_text(size = 14),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 16),
      axis.title.y  = element_text(size = 16),
      plot.title    = element_text(size = 16, hjust = 0.5),
      panel.grid = element_blank(),
      axis.line     = element_line(color = "black", linewidth = 0.6), # add axis lines
      panel.border  = element_blank(),                   # no border box
      legend.position = legend_pos
    )
  }
  return(p)
}



# ======== Scatter Plot Functions ========

plot_bulk_pb_scatter <- function(matX, meta, name1 = "Bulk", name2 = "PB", method_cor = c("pearson", "spearman"), alpha = 0.3, 
                                 point_size = 0.6, line_color = "red", abline_color = "skyblue", abline_size = 1.1, show_lm = TRUE, 
                                 title = paste0(name1,"vs",name2, "mean expression")
) {
  
  mat1 <- matX[, meta$platform=="bulk"]; mat2 <- matX[, meta$platform=="pb"]; 
  # Aggregate means per gene
  df_means <- data.frame(mean_x = as.numeric(mat1), mean_y = as.numeric(mat2))
  # df_means <- data.frame(mean_x = rowMeans(mat1), mean_y = rowMeans(mat2))
  # Compute correlations
  method_cor <- match.arg(method_cor)
  cor_val <- mean(cor(mat1, mat2, method = method_cor))
  
  # Plot
  p <- ggplot(df_means, aes(x = mean_x, y = mean_y)) +
    geom_point(alpha = alpha, size = point_size) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = abline_color,linewidth = abline_size) + 
    theme_minimal() + labs(
      title =  title,
      subtitle = paste(toupper(method_cor), "correlation =", round(cor_val, 3)),
      x = paste("Expression (", name1, ")", sep = ""),
      y = paste("Expression (", name2, ")", sep = "")) +
    theme(
      axis.text.x   = element_text(size = 14),              # larger X-axis text
      axis.text.y   = element_text(size = 14),              # larger Y-axis text
      axis.title.x  = element_text(size = 16),              # larger X-axis title
      axis.title.y  = element_text(size = 16),              # larger Y-axis title
      plot.title    = element_text(size = 16, hjust = 0.5), # larger centered title
      panel.grid    = element_blank(),
      legend.position = "none"                             
    )
  
  if (show_lm) {
    p <- p + geom_smooth(method = "lm", color = line_color, se = FALSE)
  }
  p + coord_equal()
  
  return(list(plot = p, correlation = cor_val, df = df_means))
}



# ======== Heatmap Functions (Bulk vs PB) ========

plot_heatmap <- function(expr_mat, meta_df, 
                         anno_cols = c("platform", "group"),
                         scale_data = TRUE, 
                         main_title = "Heatmap (scaled)",
                         top_n_genes = NULL,
                         backend = c("pheatmap", "ComplexHeatmap"),
                         color_palette = c("blue", "white", "red"),
                         color_range = c(-4, 0, 4)) {
  
  backend <- match.arg(backend)
  stopifnot(all(anno_cols %in% colnames(meta_df)))
  
  # Optional: top variable genes
  if (!is.null(top_n_genes)) {
    gene_vars <- matrixStats::rowVars(as.matrix(expr_mat))
    top_genes <- head(order(gene_vars, decreasing = TRUE), top_n_genes)
    expr_mat <- expr_mat[top_genes, , drop = FALSE]
  }
  
  # Annotation dataframe
  anno_df <- meta_df[, anno_cols, drop = FALSE]
  anno_df[] <- lapply(anno_df, factor)
  rownames(anno_df) <- colnames(expr_mat)
  
  # Annotation colors
  all_levels <- unique(unlist(lapply(anno_df, levels)))
  palette <- setNames(rainbow(length(all_levels)), all_levels)
  anno_colors <- lapply(anno_df, function(x) palette[levels(x)])
  
  # Scale data if requested
  mat_scaled <- if (scale_data) scale(t(expr_mat), center = TRUE, scale = TRUE) else t(expr_mat)
  mat_scaled <- t(mat_scaled)  # back to genes x samples
  
  # Clip to color range
  mat_scaled <- pmax(pmin(mat_scaled, color_range[3]), color_range[1])
  
  if (backend == "ComplexHeatmap") {
    ha <- HeatmapAnnotation(df = anno_df, col = anno_colors, annotation_name_side = "left")
    ht <- Heatmap(
      mat_scaled,
      name = "Expr.",
      col = colorRamp2(color_range, color_palette),
      top_annotation = ha,
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_gp = grid::gpar(fontsize = 8),
      column_title = main_title,
      column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      use_raster = TRUE
    )
  } else {
    heat_colors <- colorRampPalette(color_palette)(100)
    ht <- pheatmap::pheatmap(
      mat_scaled,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      show_colnames = TRUE,
      annotation_col = anno_df,
      annotation_colors = anno_colors,
      main = main_title,
      color = heat_colors,
      silent = TRUE
    )
  }
  return(ht)
}




# ======== Marker Gene Visualization Functions ========

plot_marker_slopes <- function(dgeCPM, metaInfo, marker_genes, width = 8, height = 4, save_path = NULL) 
{
  plot_df <- lapply(marker_genes, function(g) {
    data.frame(
      Gene = g,
      Bulk = as.numeric(dgeCPM[g, ][metaInfo$platform == "bulk"]),
      PB   = as.numeric(dgeCPM[g, ][metaInfo$platform == "pb"]),
      Sample = colnames(dgeCPM)
    )
  }) %>% bind_rows()
  
  slope_data <- plot_df %>%
    group_by(Gene) %>%
    do({
      model <- lm(PB ~ Bulk, data = .)
      data.frame(Slope = coef(model)[2])
    }) %>%
    mutate(Slope = round(Slope, 3))
  # color plate
  gg_colors <- scales::hue_pal()(length(marker_genes))
  p_main <- ggplot(plot_df, aes(x = Bulk, y = PB, color = Gene)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_color_manual(values = gg_colors) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none") +
    labs(x = "Bulk expression (CPM)", y = "Pseudo-bulk expression (CPM)")
  
  table_grob <- tableGrob(slope_data, rows = NULL, theme = ttheme_minimal(
    base_size = 11, core = list(
      fg_params = list(
        col = c(gg_colors, rep("black", nrow(slope_data))),
        fontface = c(rep("plain", nrow(slope_data)), rep("plain", nrow(slope_data)))
      )
    )
  ))
  
  combined <- arrangeGrob(p_main, table_grob, nrow = 1, widths = c(3, 1))
  
  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = width, height = height)
  }
  
  return(combined)
}




# ======== Mean Expression Shift Visualization ========
 
plot_mean_shift <- function(logCPM_all, metaInfo, bins = 50) {
  # Check input
  if (!"platform" %in% colnames(metaInfo))
    stop("metaInfo must contain a column named 'platform' with 'bulk' and 'pb' values.")
  
  # Compute mean shift per gene (bulk - pb)
  mean_shift <- rowMeans(logCPM_all[, metaInfo$platform == "bulk", drop = FALSE]) -
    rowMeans(logCPM_all[, metaInfo$platform == "pb", drop = FALSE])
  
  # Prepare dataframe
  df_shift <- data.frame(mean_shift = mean_shift)
  
  # Plot histogram
  p <- ggplot(df_shift, aes(x = mean_shift)) +
    geom_histogram(
      bins = bins,
      color = "white",
      fill = "#1f77b4",
      alpha = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    labs(
      x = "Mean expression shift (log2 CPM)",
      y = "Number of genes"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid = element_blank()
    )
  
  return(p)
}





# ======== Volcano Plot Functions ========

plot_volcano_multi_method <- function(exMAT, meta, spca_all, vp_res, vp_col = "PBD", 
                                      top_n = 10, highlight_genes = TRUE, pbd_cutoff = "threshold", bandwidth = .05
) {

  # --- Differential expression ---
  res_de <- run_DE(mat = exMAT, meta = meta)
  
  # --- Variance partition ---
  res_vp <- vp_res
  res_vp$Gene <- rownames(vp_res)
  
  # Peak-based threshold
  vp_cut <- find_peak_threshold(res_vp[, vp_col], bandwidth = bandwidth, select = "last")
  cut <- vp_cut[[pbd_cutoff]]
  res_vp$group <- factor(ifelse(res_vp[, vp_col] >= cut, "High-PBD", "Low-PBD"),
                         levels = c("Low-PBD", "High-PBD"))
  
  # --- sPCA ---
  sPCA_df <- data.frame(Gene = rownames(exMAT), PC1_loading = spca_all[["loadings"]][,1])
  
  # --- Merge results ---
  res_df <- data.frame(
    Gene = rownames(res_de),
    logFC = res_de$logFC,
    pval = res_de$P.Value,
    limma_flag = res_de$DE_flag,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      spca_loading = sPCA_df$PC1_loading[match(Gene, sPCA_df$Gene)],
      sPCA_flag = abs(spca_loading) != 0,
      PBD = res_vp[, vp_col][match(Gene, res_vp$Gene)],
      PBD_flag = res_vp$group[match(Gene, res_vp$Gene)],
      neglogP = -log10(pval),
      method_label = "None"
    )
  
  # Assign method labels
  res_df$method_label[res_df$limma_flag == TRUE] <- "limma-DE"
  res_df$method_label[res_df$sPCA_flag == TRUE] <- "sPCA"
  res_df$method_label[res_df$PBD_flag == "High-PBD"] <- "High-PBD"
  res_df$method_label <- factor(res_df$method_label, levels = c("High-PBD", "sPCA", "limma-DE", "None"))
  
  # Colors
  method_colors <- c("High-PBD" = "red", "sPCA" = "orange", "limma-DE" = "steelblue", "None" = "grey80")
  
  
  # Summary of method percentages
  pct_df <- data.frame(
    method = c("limma-DE", "sPCA", "High-PBD"),
    pct = c(
      mean(res_df$limma_flag) * 100,
      mean(res_df$sPCA_flag) * 100,
      mean(res_df$PBD_flag == "High-PBD") * 100
    )
  )
  
  pct_lbl <- sprintf("%s (%.1f%%)", pct_df$method, pct_df$pct)
  names(pct_lbl) <- pct_df$method
  legend_labels <- c(
    pct_lbl["High-PBD"], pct_lbl["sPCA"], pct_lbl["limma-DE"], "None"
  )
  names(legend_labels) <- c("High-PBD", "sPCA", "limma-DE", "None")
  
  # --- Volcano plot ---
  fig_volcano <- ggplot(res_df, aes(x = logFC, y = neglogP, color = method_label)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = method_colors, labels = legend_labels) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "",
      x = "logFC (Bulk vs. PB)",
      y = "-log10(p-value)",
      color = "Method"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "right"
    )
  
  # Highlight top genes if requested
  if (highlight_genes) {
    top10_highPBD <- res_df %>%
      filter(PBD_flag == "High-PBD") %>%
      arrange(desc(PBD)) %>%
      head(20)
    
    top_genes <- top10_highPBD 
    fig_volcano <- fig_volcano +
      geom_text_repel(
        data = top_genes,
        aes(label = Gene),
        size = 3.2,
        # fontface = "bold",
        color = "black"
      )
  }
  
  return(list(
    volcano_plot = fig_volcano,
    summary_df = res_df,
    pct_summary = pct_df,
    cut_pbd = vp_cut,
    top_genes = top10_highPBD
  ))
}




# ---------3. Downstream Biology Interpretation ------------------------------ 

# ======== GO Enrichment Analysis ========

run_gsva_enrichment <- function(expr_mat, genes_of_interest, metaInfo, category = "C5", subset_prefix = "GOBP", 
                                minSize = 15, maxSize = 500, fdr_cutoff = 0.05, top_n = 10,
                                title_prefix = "Bulk vs Pseudo-bulk", heatmap_option = c("ComplexHeatmap", "pheatmap")) {
  
  heatmap_option <- match.arg(heatmap_option)
  # ---- 1. Subset expression matrix ----
  expr_sub <- expr_mat[rownames(expr_mat) %in% genes_of_interest, , drop = FALSE]
  storage.mode(expr_sub) <- "numeric"
  
  # ---- 2. Gene sets ----
  msig <- msigdbr(species = "Homo sapiens", category = category)
  gene_sets <- split(msig$gene_symbol, msig$gs_name)
  prefix <- sub("^([^_]+)_.*", "\\1", names(gene_sets))
  gene_sets <- gene_sets[prefix == subset_prefix]
  
  # ---- 3. Run GSVA ----
  params <- gsvaParam(exprData = expr_sub, geneSets = gene_sets,
                      minSize = minSize, maxSize = maxSize, kcdf = "Gaussian")
  gsva_scores <- gsva(params, verbose = FALSE)
  
  # ---- 4. Metadata & grouping ----
  platform <- metaInfo$platform[match(colnames(expr_sub), metaInfo$sample_id)]
  
  # ---- 5. Pathway-level tests ----
  gsva_res <- apply(gsva_scores, 1, function(x) t.test(x ~ platform)$p.value)
  gsva_res <- data.frame(
    pathway = rownames(gsva_scores), pval = gsva_res, padj = p.adjust(gsva_res, method = "BH"),
    stringsAsFactors = FALSE
  )
  
  mean_diff <- rowMeans(gsva_scores[, platform == "bulk"]) - rowMeans(gsva_scores[, platform == "pb"])
  gsva_res$mean_diff <- mean_diff[gsva_res$pathway]
  
  # ---- 6. Clean names ----
  clean_pathway_names <- function(names_vec, max_len = 50) {
    clean <- gsub("^(GOMF_|GOBP_|HP_)", "", names_vec)
    clean <- gsub("_", " ", clean)
    clean <- trimws(clean)
    clean <- ifelse(nchar(clean) > max_len, paste0(substr(clean, 1, max_len - 3), "..."), clean)
    clean
  }
  
  clean_pathway_names <- function(names_vec, max_len = 60) {
    clean <- gsub("^(GOMF_|GOBP_|GOCC_|HP_)", "", names_vec)
    clean <- gsub("_", " ", clean)
    clean <- trimws(clean)
    
    # Convert to sentence case (first letter uppercase, rest lowercase)
    clean <- paste0(toupper(substr(clean, 1, 1)), 
                    tolower(substr(clean, 2, nchar(clean))))
    
    # Force uppercase for specific biological terms (case-insensitive)
    terms_to_upper <- c("dna", "rna", "mrna", "trna", "rrna", "atp", "gtp", 
                        "gdp", "nad", "fad", "camp", "cgmp", "g1", "g2")
    
    for (term in terms_to_upper) {
      pattern <- paste0("\\b", term, "\\b")
      clean <- gsub(pattern, toupper(term), clean, ignore.case = TRUE)
    }
    
    # Special handling for mixed case like "miRNA" (mi in lowercase, RNA in uppercase)
    clean <- gsub("\\bmirna\\b", "miRNA", clean, ignore.case = TRUE)
    clean <- gsub("\\bsirna\\b", "siRNA", clean, ignore.case = TRUE)
    clean <- gsub("\\bncrna\\b", "ncRNA", clean, ignore.case = TRUE)
    clean <- gsub("mrna", "mRNA", clean)
    clean <- gsub("trna", "tRNA", clean)
    clean <- gsub("rrna", "rRNA", clean)
    clean <- gsub("RRNA", "rRNA", clean)
    
    
    # Truncate if too long
    ifelse(nchar(clean) > max_len, 
           paste0(substr(clean, 1, max_len - 3), "..."), 
           clean)
  }
  
  gsva_res$pathway_clean <- clean_pathway_names(gsva_res$pathway)
  
  # ---- 7. Top enriched terms ----
  top_bulk <- gsva_res %>%
    filter(mean_diff > 0, padj < fdr_cutoff) %>%
    arrange(desc(mean_diff)) %>%
    slice_head(n = top_n) %>%
    mutate(group = "Bulk enriched")
  
  top_pb <- gsva_res %>%
    filter(mean_diff < 0, padj < fdr_cutoff) %>%
    arrange(mean_diff) %>%
    slice_head(n = top_n) %>%
    mutate(group = "PB enriched")
  
  top_terms <- bind_rows(top_bulk, top_pb)
  
  
  # ---- 9. Barplot ----
  p_barplot <- ggplot(top_terms, aes(x = reorder(pathway_clean, mean_diff), y = mean_diff, fill = group)) +
    geom_bar(stat = "identity") + coord_flip() +
    facet_wrap(~group, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Bulk enriched" = "#1f78b4", "PB enriched"   = "#33a02c")) +
    labs(title = " ", x = "Pathway / GO Process", y = "Enrichment score (NES)", fill = "Group") +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16)
    )
  
  # ---- 11. Heatmap ----
  var_pathways <- unique(top_terms$pathway)
  gsva_sub <- gsva_scores[var_pathways, , drop = FALSE]
  rownames(gsva_sub) <- gsva_res$pathway_clean[match(rownames(gsva_sub), gsva_res$pathway)]
  
  
  # font size for row names
  rn_fontsize <- 10
  if (heatmap_option == "ComplexHeatmap") {
    # compute required text width
    rn_width <- ComplexHeatmap::max_text_width(rownames(gsva_sub),
                                               gp = gpar(fontsize = rn_fontsize)
    )
    
    ha <- HeatmapAnnotation(Platform = platform,
                            col = list(Platform = c("bulk" = "#1f78b4", "pb" = "#33a02c"))
    )
    col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
    p_heatmap <- Heatmap(gsva_sub,
                         name = "GSVA score",
                         top_annotation = ha,
                         col = col_fun,
                         show_column_names = FALSE,
                         row_names_side = "left",
                         row_names_gp = gpar(fontsize = rn_fontsize),
                         row_names_max_width = rn_width + unit(2, "mm"),
                         cluster_rows = FALSE,
                         cluster_columns = TRUE,
                         # column_title = paste0("Top GSVA pathways (", title_prefix, ")"),
                         heatmap_legend_param = list(title = "GSVA score",direction = "horizontal")
    )
    
    draw(p_heatmap,
         heatmap_legend_side = "top",
         annotation_legend_side = "top",
         padding = unit(c(5, 5, 5, 5), "mm"))
    
  } else {
    # pheatmap
    annotation_col <- data.frame(Platform = platform)
    rownames(annotation_col) <- colnames(gsva_sub)
    
    ann_colors <- list(Platform = c("bulk" = "#1f78b4", "pb" = "#33a02c"))
    # enlarge left margin
    op <- par(no.readonly = TRUE)
    par(mar = c(5, 12, 4, 2))
    p_heatmap <- pheatmap(gsva_sub,
                          annotation_col = annotation_col,
                          annotation_colors = ann_colors,
                          show_colnames = FALSE,
                          scale = "row",
                          fontsize_row = rn_fontsize,
                          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                          # main = paste0("Top GSVA pathways (", title_prefix, ")"),
                          legend = TRUE,
                          silent = TRUE
    )
    par(op)
  }
  
  # ---- 12. Return outputs ----
  list(gsva_scores = gsva_scores,
       gsva_res = gsva_res,
       top_terms = top_terms,
       barplot = p_barplot,
       heatmap = p_heatmap
  )
}


# Improved version
run_gsva_enrichment2 <- function(expr_mat, genes_of_interest, metaInfo, category = "C5", subset_prefix = "GOBP", 
                                 minSize = 15, maxSize = 500, fdr_cutoff = 0.05, top_n = 10, top_cat.n = 15,
                                 title_prefix = "Platform-Divergent Genes", 
                                 heatmap_option = c("ComplexHeatmap", "pheatmap"),
                                 plot_width = 12, plot_height = 8, path) {
  
  heatmap_option <- match.arg(heatmap_option)
  
  # ---- 1. Subset expression matrix ----
  cat("Subsetting expression matrix...\n")
  expr_sub <- expr_mat[rownames(expr_mat) %in% genes_of_interest, , drop = FALSE]
  storage.mode(expr_sub) <- "numeric"
  cat(sprintf("   %d genes × %d samples\n", nrow(expr_sub), ncol(expr_sub)))
  
  # ---- 2. Gene sets ----
  cat("Loading gene sets...\n")
  msig <- msigdbr(species = "Homo sapiens", category = category) 
  gene_sets <- split(msig$gene_symbol, msig$gs_name)
  prefix <- sub("^([^_]+)_.*", "\\1", names(gene_sets))
  gene_sets <- gene_sets[prefix == subset_prefix]
  cat(sprintf("   %d gene sets loaded\n", length(gene_sets)))
  
  # ---- 3. Run GSVA ----
  cat("Running GSVA...\n")
  params <- gsvaParam(exprData = expr_sub, geneSets = gene_sets,
                      minSize = minSize, maxSize = maxSize, kcdf = "Gaussian")
  gsva_scores <- gsva(params, verbose = FALSE)
  
  # ---- 4. Metadata & grouping ----
  platform <- metaInfo$platform[match(colnames(expr_sub), metaInfo$sample_id)]
  
  # ---- 5. Pathway-level tests ----
  cat("Performing pathway-level tests...\n")
  gsva_res <- apply(gsva_scores, 1, function(x) t.test(x ~ platform)$p.value)
  gsva_res <- data.frame(
    pathway = rownames(gsva_scores), 
    pval = gsva_res, 
    padj = p.adjust(gsva_res, method = "BH"), 
    stringsAsFactors = FALSE
  )
  
  mean_diff <- rowMeans(gsva_scores[, platform == "bulk"]) - 
    rowMeans(gsva_scores[, platform == "pb"])
  gsva_res$mean_diff <- mean_diff[gsva_res$pathway]
  
  
  
  clean_pathway_names <- function(names_vec, max_len = 60) {
    clean <- gsub("^(GOMF_|GOBP_|GOCC_|HP_)", "", names_vec)
    clean <- gsub("_", " ", clean)
    clean <- trimws(clean)
    
    # Convert to sentence case (first letter uppercase, rest lowercase)
    clean <- paste0(toupper(substr(clean, 1, 1)), 
                    tolower(substr(clean, 2, nchar(clean))))
    
    # Force uppercase for specific biological terms (case-insensitive)
    terms_to_upper <- c("dna", "rna", "mrna", "trna", "rrna", "atp", "gtp", 
                        "gdp", "nad", "fad", "camp", "cgmp", "g1", "g2")
    
    for (term in terms_to_upper) {
      pattern <- paste0("\\b", term, "\\b")
      clean <- gsub(pattern, toupper(term), clean, ignore.case = TRUE)
    }
    
    # Special handling for mixed case like "miRNA" (mi in lowercase, RNA in uppercase)
    clean <- gsub("\\bmirna\\b", "miRNA", clean, ignore.case = TRUE)
    clean <- gsub("\\bsirna\\b", "siRNA", clean, ignore.case = TRUE)
    clean <- gsub("\\bncrna\\b", "ncRNA", clean, ignore.case = TRUE)
    clean <- gsub("mrna", "mRNA", clean)
    clean <- gsub("trna", "tRNA", clean)
    clean <- gsub("rrna", "rRNA", clean)
    clean <- gsub("RRNA", "rRNA", clean)
    
    
    # Truncate if too long
    ifelse(nchar(clean) > max_len, 
           paste0(substr(clean, 1, max_len - 3), "..."), 
           clean)
  }
  
  gsva_res$pathway_clean <- clean_pathway_names(gsva_res$pathway)
  
  
  # ---- 7. Categorize GO terms ----
  categorize_go_terms <- function(terms) {
    categories <- case_when(
      # 1. REGULATION TERMS (highest priority)
      str_detect(toupper(terms), "^GOBP_(POSITIVE|NEGATIVE)_REGULATION_OF_") ~ "Regulation (Positive/Negative)",
      str_detect(toupper(terms), "^GOBP_REGULATION_OF_") ~ "Regulation",
      
      # 2. SIGNALING PATHWAYS (specific before general)
      str_detect(toupper(terms), "(?i)GTPASE_|RAS_|RHO_|SMALL_GTPASE") ~ "Signaling Pathways (GTPase)",
      str_detect(toupper(terms), "(?i)AUTOPHOSPHORYLATION|PHOSPHORYLATION|KINASE_ACTIVITY|PHOSPHATASE") ~ "Kinase/Phosphorylation Signaling",
      str_detect(toupper(terms), "(?i)_SIGNALING_PATHWAY$|_SIGNAL_TRANSDUCTION$|SIGNAL_CASCADE$") ~ "Signaling Pathways",
      str_detect(toupper(terms), "(?i)MAPK|NF_KAPPAB|WNT|NOTCH|PI3K|TOR|JAK_STAT|SMAD|HEDGEHOG|SMOOTHENED|TGF_BETA|EGFR|VEGF|INSULIN_RECEPTOR|ERBB|BMP|TUMOR_NECROSIS_FACTOR|INTERFERON|CALCIUM_MEDIATED") ~ "Signaling Pathways",
      
      # 3. RESPONSE & STRESS
      str_detect(toupper(terms), "GOBP_RESPONSE_TO_") ~ "Response to Stimuli",
      str_detect(toupper(terms), "(?i)STRESS_RESPONSE|DNA_DAMAGE_RESPONSE|UNFOLDED_PROTEIN_RESPONSE|OXIDATIVE_STRESS|HEAT_SHOCK") ~ "Stress Response",
      
      # 4. RNA & GENE EXPRESSION
      str_detect(toupper(terms), "(?i)NCRNA_|MIRNA_|SIRNA_|REGULATORY_NCRNA") ~ "RNA Processing & ncRNA",
      str_detect(toupper(terms), "(?i)RNA_MODIFICATION|RNA_METHYLATION|RNA_PROCESSING|RNA_SPLICING|SPLICEOSOME|SNRNP_") ~ "RNA Processing",
      str_detect(toupper(terms), "(?i)MRNA_METABOLIC|MRNA_CATABOLIC|MRNA_PROCESSING|MRNA_EXPORT|MRNA_SPLICING") ~ "mRNA Metabolism",
      str_detect(toupper(terms), "(?i)TRANSCRIPTION_|RNA_POLYMERASE|PREINITIATION_COMPLEX|GENE_EXPRESSION") ~ "Transcription & Gene Expression",
      str_detect(toupper(terms), "(?i)EPIGENETIC|CHROMATIN_REMODELING|HETEROCHROMATIN|DNA_METHYLATION|HISTONE_MODIFICATION") ~ "Epigenetic Regulation",
      
      # 5. TRANSLATION & PROTEIN SYNTHESIS
      str_detect(toupper(terms), "(?i)TRANSLATION|TRANSLATIONAL_|RIBOSOM|RIBONUCLEOPROTEIN|TRNA_|RRNA_") ~ "Translation & Protein Synthesis",
      
      # 6. PROTEIN PROCESSING & MODIFICATION
      str_detect(toupper(terms), "(?i)GLYCOSYLATION|O_LINKED|N_LINKED|PROTEIN_MODIFICATION") ~ "Protein Modification",
      str_detect(toupper(terms), "(?i)UBIQUITINATION|SUMOYLATION|ACETYLATION|METHYLATION") ~ "Protein Post-translational Modification",
      str_detect(toupper(terms), "(?i)PROTEOLYSIS|PROTEASOME|PROTEIN_FOLDING|PROTEIN_MATURATION|CHAPERONE") ~ "Protein Processing & Folding",
      str_detect(toupper(terms), "(?i)PROTEIN_CATABOLIC|PROTEASOMAL|UBIQUITIN_DEPENDENT") ~ "Protein Degradation",
      
      # 7. TRANSPORT & TRAFFICKING
      str_detect(toupper(terms), "(?i)PROTEIN_LOCALIZATION|PROTEIN_TARGETING|PROTEIN_EXPORT|PROTEIN_IMPORT") ~ "Protein Localization & Targeting",
      str_detect(toupper(terms), "(?i)ION_TRANSPORT|MEMBRANE_TRANSPORT|CHANNEL_ACTIVITY|TRANSPORTER_ACTIVITY|SODIUM_ION|POTASSIUM_ION|CALCIUM_ION") ~ "Ion & Membrane Transport",
      str_detect(toupper(terms), "(?i)_TRANSPORT$|_TRANSPORT_|TRANSMEMBRANE_TRANSPORT|VESICLE_|ENDOCYTOSIS|EXOCYTOSIS|SECRETION|TARGETING|LOCALIZATION|BUDDING|FUSION|DOCKING") ~ "Transport & Trafficking",
      
      # 8. METABOLIC PROCESSES
      str_detect(toupper(terms), "(?i)MITOCHONDR|MITOPHAGY|MITOCHONDRIAL_ELECTRON|MITOCHONDRIAL_GENE|MITOCHONDRIAL_TRANSLATION") ~ "Mitochondrial Processes",
      str_detect(toupper(terms), "(?i)OXIDATIVE_PHOSPHORYLATION|ELECTRON_TRANSPORT|ATP_SYNTHESIS|PROTON_MOTIVE|AEROBIC_RESPIRATION|ATP_BIOSYNTHETIC") ~ "Energy Metabolism",
      str_detect(toupper(terms), "(?i)NUCLEOTIDE_|PURINE|PYRIMIDINE|RIBOSE_PHOSPHATE") ~ "Nucleotide Metabolism",
      str_detect(toupper(terms), "(?i)CARBOHYDRATE_METABOL|GLUCOSE_METABOL|GLYCOLYSIS|GLUCONEOGENESIS") ~ "Carbohydrate Metabolism",
      str_detect(toupper(terms), "(?i)LIPID_METABOL|FATTY_ACID|PHOSPHOLIPID|GLYCEROLIPID|STEROID|CHOLESTEROL|SPHINGOLIPID") ~ "Lipid Metabolism",
      str_detect(toupper(terms), "(?i)AMINO_ACID_METABOL|PROTEIN_METABOL|PEPTIDE_METABOL") ~ "Amino Acid/Protein Metabolism",
      str_detect(toupper(terms), "(?i)_METABOLIC_PROCESS$|_BIOSYNTHETIC_PROCESS$|_CATABOLIC_PROCESS$") ~ "Metabolic Processes",
      
      # 9. CELL CYCLE & DIVISION
      str_detect(toupper(terms), "(?i)CELL_CYCLE|MITOTIC|MEIOTIC|CYTOKINESIS|CHROMOSOME_SEGREGATION|SISTER_CHROMATID|SPINDLE|CHECKPOINT|NUCLEAR_DIVISION|CENTROSOME|TELOMERE|DNA_REPLICATION") ~ "Cell Cycle & Division",
      
      # 10. DEVELOPMENT & DIFFERENTIATION (specific to general)
      str_detect(toupper(terms), "(?i)SENSORY_PERCEPTION|SENSORY_SYSTEM|SENSORY_ORGAN|CHEMICAL_STIMULUS|MECHANICAL_STIMULUS|LIGHT_STIMULUS") ~ "Sensory & Perception",
      str_detect(toupper(terms), "(?i)EYE_MORPHOGENESIS|CAMERA_TYPE_EYE|RETINA|LENS|OPTIC|VISUAL") ~ "Visual System Development",
      str_detect(toupper(terms), "(?i)MYOTUBE|MYOBLAST|MUSCLE_CELL|STRIATED_MUSCLE|SKELETAL_MUSCLE|CARDIAC_MUSCLE|SMOOTH_MUSCLE") ~ "Muscle Development",
      str_detect(toupper(terms), "(?i)SKELETAL_SYSTEM|BONE|CARTILAGE|OSSIFICATION|OSTEO|CHONDRO") ~ "Skeletal System Development",
      str_detect(toupper(terms), "(?i)NEUROGENESIS|GLIOGENESIS|MYOGENESIS|OSTEOGENESIS|CHONDROGENESIS|ANGIOGENESIS|VASCULOGENESIS|HEMATOPOIESIS|HEMOPOIESIS") ~ "Cell Type Differentiation",
      str_detect(toupper(terms), "(?i)EMBRYONIC|FETAL|IN_UTERO") ~ "Embryonic Development",
      str_detect(toupper(terms), "(?i)TISSUE_DEVELOPMENT|ORGAN_MORPHOGENESIS|TUBE_DEVELOPMENT|BRANCHING_MORPHOGENESIS|PATTERN_SPECIFICATION|CELL_FATE") ~ "Tissue & Organ Development",
      str_detect(toupper(terms), "(?i)ORGAN_GROWTH|TISSUE_GROWTH|DEVELOPMENTAL_GROWTH") ~ "Growth & Organogenesis",
      str_detect(toupper(terms), "(?i)_DEVELOPMENT$|_DIFFERENTIATION$|MORPHOGENESIS$") ~ "Development & Differentiation",
      
      # 11. NERVOUS SYSTEM
      str_detect(toupper(terms), "(?i)DENDRITE|AXON|AXONOGENESIS|NEURITE") ~ "Neuronal Structure & Development",
      str_detect(toupper(terms), "(?i)POSTSYNAPSE|PRESYNAPSE|SYNAPTIC_VESICLE|SYNAPTIC_TRANSMISSION|SYNAPTIC_PLASTICITY") ~ "Synaptic Function",
      str_detect(toupper(terms), "(?i)NEURON|NEURO|SYNAPSE|SYNAPTIC|MYELIN|GLIAL|BRAIN|CORTEX|CEREBRAL|FOREBRAIN|HINDBRAIN|NEUROTRANSMITTER") ~ "Nervous System",
      
      # 12. IMMUNE SYSTEM
      str_detect(toupper(terms), "(?i)IMMUNE|IMMUNITY|ANTIGEN|CYTOKINE|INTERLEUKIN|INTERFERON|LYMPHOCYTE|LEUKOCYTE|T_CELL|B_CELL|MYELOID|MACROPHAGE|DEFENSE_RESPONSE|ANTIMICROBIAL|PHAGOCYTOSIS|INFLAMMATORY_RESPONSE|HUMORAL_IMMUNE|IMMUNOGLOBULIN") ~ "Immune System & Defense",
      
      # 13. CELL STRUCTURE & ADHESION
      str_detect(toupper(terms), "(?i)CYTOSKELETON|ACTIN|MICROTUBULE|FILAMENT") ~ "Cytoskeleton Organization",
      str_detect(toupper(terms), "(?i)MOTILITY|MIGRATION|CHEMOTAXIS|LOCOMOTION|MOTILE|CILIUM|FLAGELLUM|LAMELLIPODIUM|SPREADING|AMEBOIDAL") ~ "Cell Motility & Migration",
      str_detect(toupper(terms), "(?i)ADHESION|JUNCTION|TIGHT_JUNCTION|ADHERENS|DESMOSOME|FOCAL_ADHESION|CELL_MATRIX|CELL_CELL|HOMOPHILIC|INTEGRIN") ~ "Cell Adhesion & Junctions",
      str_detect(toupper(terms), "(?i)CONTRACTION|MUSCLE_CONTRACTION") ~ "Muscle Contraction",
      
      # 14. HOMEOSTASIS & MAINTENANCE
      str_detect(toupper(terms), "(?i)HOMEOSTASIS|MAINTENANCE|REGENERATION") ~ "Homeostasis & Maintenance",
      
      # 15. CELL DEATH & SENESCENCE
      str_detect(toupper(terms), "(?i)APOPTOSIS|CELL_DEATH|PROGRAMMED_CELL_DEATH") ~ "Apoptosis & Cell Death",
      str_detect(toupper(terms), "(?i)SENESCENCE|AGING|CELLULAR_SENESCENCE") ~ "Cell Senescence & Aging",
      
      # 16. DNA REPAIR & MAINTENANCE
      str_detect(toupper(terms), "(?i)DNA_REPAIR|BASE_EXCISION|NUCLEOTIDE_EXCISION|DOUBLE_STRAND_BREAK|RECOMBINATIONAL_REPAIR") ~ "DNA Repair & Maintenance",
      
      # 17. REPRODUCTION
      str_detect(toupper(terms), "(?i)REPRODUCTION|GAMETE|SPERM|OVARY|TESTIS|SEXUAL|FERTILIZATION|SPERMATID|OOCYTE|MEIOSIS") ~ "Reproduction & Gametogenesis",
      
      # 18. ORGANELLE ORGANIZATION
      str_detect(toupper(terms), "(?i)ORGANELLE_ORGANIZATION|ORGANELLE_ASSEMBLY|MEMBRANE_ORGANIZATION|NUCLEUS_ORGANIZATION|ENDOMEMBRANE_SYSTEM|GOLGI_ORGANIZATION|ENDOPLASMIC_RETICULUM_ORGANIZATION") ~ "Organelle Organization",
      
      # 19. PROTEIN COMPLEX ASSEMBLY
      str_detect(toupper(terms), "(?i)COMPLEX_ASSEMBLY|RESPIRATORY_CHAIN_COMPLEX|NADH_DEHYDROGENASE_COMPLEX") ~ "Protein Complex Assembly",
      
      # 20. BEHAVIOR & COGNITION
      str_detect(toupper(terms), "(?i)BEHAVIOR|LEARNING|MEMORY|COGNITION|RHYTHMIC|CIRCADIAN|ASSOCIATIVE_LEARNING") ~ "Behavior & Cognition",
      
      # 21. AUTOPHAGY & LYSOSOMAL
      str_detect(toupper(terms), "(?i)AUTOPHAGY|MACROAUTOPHAGY|AUTOPHAGOSOME|LYSOSOMAL|VACUOLE|PHAGOSOME") ~ "Autophagy & Lysosomal Processes",
      
      # 22. VIRAL PROCESSES
      str_detect(toupper(terms), "(?i)VIRAL_|VIRUS|HOST_INTERACTION|SYMBIONT") ~ "Viral Processes & Host Interactions",
      
      # Default
      TRUE ~ "Other Biological Processes"
    )
    
    return(categories)
  }
  
  gsva_res$category <- categorize_go_terms(gsva_res$pathway)
  
  # ---- 8. Top enriched terms ----
  cat("Selecting top enriched terms...\n")
  top_bulk <- gsva_res %>%
    filter(mean_diff > 0, padj < fdr_cutoff) %>%
    arrange(desc(mean_diff)) %>%
    slice_head(n = top_n) %>%
    mutate(group = "Bulk enriched")
  
  top_pb <- gsva_res %>%
    filter(mean_diff < 0, padj < fdr_cutoff) %>%
    arrange(mean_diff) %>%
    slice_head(n = top_n) %>%
    mutate(group = "PB enriched")
  
  top_terms <- bind_rows(top_bulk, top_pb)
  
  # ---- 9. Enhanced Barplot ----
  cat("Creating enhanced visualizations...\n")
  p_barplot <- ggplot(top_terms, aes(x = reorder(pathway_clean, mean_diff), y = mean_diff, fill = group)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.85) + 
    coord_flip() +
    facet_grid(group ~ ., scales = "free_y", space = "free") +
    scale_fill_manual(values = c("Bulk enriched" = "#2E86AB", 
                                 "PB enriched" = "#A23B72")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = sprintf("%.2f", mean_diff), 
                  hjust = ifelse(mean_diff > 0, -0.2, 1.2)),
              size = 3.5, color = "black") +
    labs(title = title_prefix,
         subtitle = paste("FDR <", fdr_cutoff, "| Top", top_n, "terms per group"),
         x = "Pathway/ GO Process", 
         y = "Enrichment Score (NES)",
         fill = "Enrichment") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold", color = "white"),
      strip.background = element_rect(fill = "gray30"),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40", size = 12),
      panel.spacing = unit(.2, "lines")
    )
  
  
  # ---- 10. Volcano Plot ----
  volcano_data <- gsva_res %>%
    mutate(
      significance = case_when(
        padj < 0.001 & abs(mean_diff) > 0.5 ~ "High",
        padj < 0.05 & abs(mean_diff) > 0.3 ~ "Medium",
        TRUE ~ "NotSig"
      ),
      label = ifelse(padj < 0.001 & abs(mean_diff) > 0.8, pathway_clean, "")
    )
  
  p_volcano <- ggplot(volcano_data, aes(x = mean_diff, y = -log10(padj))) +
    geom_point(aes(color = significance, size = significance), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = filter(volcano_data, label != ""),
      aes(label = label),
      size = 3.5,
      max.overlaps = 10,
      box.padding = 0.5,
      segment.color = "gray50"
    ) +
    scale_color_manual(
      values = c("High" = "#D62728", "Medium" = "#FF7F0E", "NotSig" = "gray50"),
      name = "Significance"
    ) +
    scale_size_manual(
      values = c("High" = 3, "Medium" = 2, "NotSig" = 1),
      guide = "none"
    ) +
    labs(
      title = " ",
      subtitle = paste(nrow(gsva_res), "GO terms"),
      x = "Mean Difference (Bulk - Pseudo-bulk)",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      panel.grid.major = element_line(color = "gray90")
    )
  
  
  # ---- 8. Top enriched terms ----
  cat("Selecting top enriched terms...\n")
  top_bulk <- gsva_res %>%
    filter(mean_diff > 0, padj < fdr_cutoff) %>%
    arrange(desc(mean_diff)) %>%
    slice_head(n = top_cat.n) %>%
    mutate(group = "Bulk enriched")
  
  top_pb <- gsva_res %>%
    filter(mean_diff < 0, padj < fdr_cutoff) %>%
    arrange(mean_diff) %>%
    slice_head(n = top_cat.n) %>%
    mutate(group = "PB enriched")
  
  top_terms1 <- bind_rows(top_bulk, top_pb)
  
  bubble_data <- top_terms1 %>%
    filter(padj < fdr_cutoff) %>%
    group_by(category, sign = ifelse(mean_diff > 0, "Bulk", "Pseudo-bulk")) %>%
    slice_max(order_by = abs(mean_diff), n = 2) %>%
    ungroup()
  
  if(nrow(bubble_data) > 0) {
    
    # Get unique categories
    categories <- unique(bubble_data$category)
    
    # Create a named color vector (optional - if you want colored strips)
    library(viridis)
    strip_colors <- viridis::viridis(length(categories))
    names(strip_colors) <- categories
    
    p_bubble <- ggplot(bubble_data, 
                       aes(x = factor(sign, levels = c("Bulk", "Pseudo-bulk")), 
                           y = reorder(pathway_clean, mean_diff))) +
      geom_point(aes(size = abs(mean_diff), 
                     color = -log10(padj)), 
                 alpha = 0.8) +
      scale_color_gradientn(
        colors = c("blue", "cyan", "green", "yellow", "orange", "red"),
        name = "-log10(padj)"
      ) +
      scale_size_continuous(
        range = c(3, 10),
        name = "NES"
      ) +
      # Use switch = "y" to move strips to left side
      facet_grid(
        category ~ ., 
        scales = "free_y", 
        space = "free",
        switch = "y",  # This moves strips to left
        labeller = labeller(category = ~.x)  # Keep original category names
      ) +
      labs(
        title = " ",
        subtitle = " ",
        x = " ",
        y = "Pathway/GO Biological Process"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        # Use strip.text.y.left for left-positioned strips
        strip.text.y.left = element_text(
          angle = 0, 
          size = 10, 
          face = "bold", 
          margin = margin(l = 5, r = 5),
          hjust = 0.5  # Center text horizontally
        ),
        strip.background = element_rect(
          fill = "#e8f5e9", 
          color = "black", 
          size = 0.01
        ),
        strip.placement = "outside",  # Ensure strips are outside plot area
        panel.grid.major = element_line(color = "gray90"),
        legend.position = "right"
      )
    
  } else {
    p_bubble <- NULL
  }
  
  
  # ---- 12. Category Distribution Plot ----
  category_summary <- top_terms1 %>%
    filter(padj < fdr_cutoff) %>%
    mutate(direction = ifelse(mean_diff > 0, "Bulk", "Pseudo-bulk")) %>%
    group_by(category, direction) %>%
    summarise(
      n_terms = n(),
      mean_effect = mean(abs(mean_diff)),
      .groups = 'drop'
    )
  
  if(nrow(category_summary) > 0) {
    p_category <- ggplot(category_summary, 
                         aes(x = reorder(category, n_terms), 
                             y = n_terms, fill = direction)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
      coord_flip() +
      scale_fill_manual(
        values = c("Bulk" = "#2E86AB", "Pseudo-bulk" = "#A23B72"),
        name = "Enriched in"
      ) +
      geom_text(
        aes(label = n_terms, group = direction),
        position = position_dodge(width = 0.9),
        hjust = -0.2, size = 3.5
      ) +
      labs(
        title = "Distribution of Enriched GO Terms by Category",
        x = "Biological Category",
        y = "Number of Enriched Terms"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        axis.text.y = element_text(size = 11)
      )
  } else {
    p_category <- NULL
  }
  
  # ---- 13. Enhanced Heatmap ----
  var_pathways <- unique(top_terms$pathway)
  gsva_sub <- gsva_scores[var_pathways, , drop = FALSE]
  rownames(gsva_sub) <- gsva_res$pathway_clean[match(rownames(gsva_sub), gsva_res$pathway)]
  
  # Order samples by platform
  platform_order <- order(platform)
  gsva_sub <- gsva_sub[, platform_order]
  platform_ordered <- platform[platform_order]
  
  if(heatmap_option == "ComplexHeatmap") {
    # Create annotations
    ha <- HeatmapAnnotation(
      Platform = platform_ordered,
      col = list(Platform = c("bulk" = "blue", "pb" = "red")),
      annotation_name_side = "right",
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        Platform = list(title = "Platform", 
                        labels = c("bulk" = "Bulk", "pb" = "Pseudo-bulk"))
      )
    )
    
    
    row_categories <- gsva_res$category[match(rownames(gsva_sub), gsva_res$pathway_clean)]
    unique_cats <- unique(row_categories)
    n_cats <- length(unique_cats)
    
    # Generate distinct colors with randomcoloR
    cat_colors <- setNames(
      randomcoloR::distinctColorPalette(n_cats),
      unique_cats
    )
    
    ra <- rowAnnotation(
      Category = row_categories,
      col = list(Category = cat_colors),
      show_annotation_name = FALSE,
      annotation_name_side = "top"
    )
    
    # Color scale
    col_fun <- colorRamp2(
      breaks = c(min(gsva_sub), 0, max(gsva_sub)),
      colors = c("#2E86AB", "white", "#A23B72")
    )
    
    # Create heatmap
    p_heatmap <- Heatmap(
      gsva_sub,
      name = "GSVA\nScore",
      top_annotation = ha,
      left_annotation = ra,
      col = col_fun,
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10),
      row_names_max_width = max_text_width(
        rownames(gsva_sub),
        gp = gpar(fontsize = 12)
      ) * 1.2,  # Add 20% padding
      column_split = platform_ordered,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      cluster_row_slices = FALSE,
      show_row_dend = FALSE,
      row_title = NULL,
      column_title = title_prefix,
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      heatmap_legend_param = list(
        title = "GSVA Score",
        direction = "horizontal",
        title_position = "topcenter"
      )
    )
    
    # Draw heatmap
    draw(p_heatmap,
         heatmap_legend_side = "top",
         annotation_legend_side = "top",
         padding = unit(c(2, 2, 2, 2), "cm"))
    
  } else {
    # pheatmap version
    annotation_col <- data.frame(
      Platform = platform_ordered,
      row.names = colnames(gsva_sub)
    )
    
    annotation_row <- data.frame(
      Category = row_categories,
      row.names = rownames(gsva_sub)
    )
    
    ann_colors <- list(
      Platform = c("bulk" = "#2E86AB", "pb" = "#A23B72"),
      Category = cat_colors
    )
    
    p_heatmap <- pheatmap(
      gsva_sub,
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      show_colnames = FALSE,
      show_rownames = TRUE,
      scale = "row",
      fontsize_row = 9,
      fontsize_col = 8,
      color = colorRampPalette(c("#2E86AB", "white", "#A23B72"))(100),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      border_color = NA,
      main = paste(title_prefix, "- Pathway Enrichment"),
      legend = TRUE,
      silent = FALSE
    )
  }
  
  # ---- 14. Summary Statistics ----
  cat("\n=== ENRICHMENT ANALYSIS SUMMARY ===\n")
  cat(sprintf("Total GO terms tested: %d\n", nrow(gsva_res)))
  cat(sprintf("Significant terms (FDR < %.3f): %d\n", 
              fdr_cutoff, sum(gsva_res$padj < fdr_cutoff)))
  cat(sprintf("Bulk-enriched terms: %d\n", sum(gsva_res$padj < fdr_cutoff & gsva_res$mean_diff > 0)))
  cat(sprintf("Pseudo-bulk-enriched terms: %d\n", sum(gsva_res$padj < fdr_cutoff & gsva_res$mean_diff < 0)))
  cat(sprintf("Most significant FDR: %.2e\n", min(gsva_res$padj, na.rm = TRUE)))
  cat(sprintf("Largest effect size: %.3f\n", max(abs(gsva_res$mean_diff), na.rm = TRUE)))
  
  # ---- 15. Return outputs ----
  results <- list(
    gsva_scores = gsva_scores,
    gsva_results = gsva_res,
    top_terms = top_terms,
    barplot = p_barplot,
    volcano_plot = p_volcano,
    bubble_plot = p_bubble,
    category_plot = p_category,
    heatmap = p_heatmap,
    summary = list(
      total_tested = nrow(gsva_res),
      significant = sum(gsva_res$padj < fdr_cutoff),
      bulk_enriched = sum(gsva_res$padj < fdr_cutoff & gsva_res$mean_diff > 0),
      pb_enriched = sum(gsva_res$padj < fdr_cutoff & gsva_res$mean_diff < 0),
      min_fdr = min(gsva_res$padj, na.rm = TRUE),
      max_effect = max(abs(gsva_res$mean_diff), na.rm = TRUE)
    )
  )
  
  return(results)
}




# -------- 4. UTILITY AND HELPER FUNCTIONS -----------------------------------

# ======== QC Functions ========

is_outlier <- function(SeuratObject, metric, nmads = 5) {
  # Ensure the metric is available in meta.data
  if (!metric %in% colnames(SeuratObject@meta.data)) {
    stop(paste("Metric", metric, "not found in meta.data"))
  }
  # Extract values
  M <- SeuratObject@meta.data[[metric]]
  
  # Robust thresholds
  med <- median(M, na.rm = TRUE)
  mad_val <- mad(M, na.rm = TRUE)
  lower <- med - nmads * mad_val
  upper <- med + nmads * mad_val
  
  # Flag outliers
  outliers <- (M < lower) | (M > upper)
  names(outliers) <- colnames(SeuratObject)
  
  return(outliers)
}

#  QC filter
qc_filter <- function(SeuratObject, nmads = 5, metrics = c("nFeature_RNA", "nCount_RNA", "percent.mito")) {
  outlier_flags <- lapply(metrics, function(m) is_outlier(SeuratObject, m, nmads))
  names(outlier_flags) <- metrics
  
  # Combine: flag cell as outlier if TRUE for any metric
  combined_outliers <- Reduce("|", outlier_flags)
  
  message("Outliers detected per metric:")
  print(sapply(outlier_flags, sum))
  message("Total outliers to remove: ", sum(combined_outliers))
  
  # Subset Seurat object to retain only high-quality cells
  filtered <- SeuratObject[, !combined_outliers]
  
  return(list(seurat_orig = SeuratObject,  
              seurat_filt = filtered))
}



# ======== Plot Saving Function ========

save_plot <- function(plot_obj, filename, path = ".", width = 10, height = 8, dpi = 300) {
  # ensure output directory exists
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  outfile <- file.path(path, filename)
  if (inherits(plot_obj, "Heatmap") || inherits(plot_obj, "HeatmapList")) {
    # ComplexHeatmap objects
    ext <- tools::file_ext(outfile)
    if (ext == "") {
      outfile <- paste0(outfile, ".png")
      ext <- "png"
    }
    
    if (ext %in% c("png", "tiff", "jpeg", "jpg")) {
      do.call(ext, list(outfile, width = width, height = height, units = "in", res = dpi))
      draw(plot_obj)
      dev.off()
    } else if (ext == "pdf") {
      pdf(outfile, width = width, height = height)
      draw(plot_obj)
      dev.off()
    } else {
      stop("Unsupported file type for ComplexHeatmap. Use png/tiff/jpeg/pdf.")
    }
    
  } else {
    # ggplot2/patchwork objects
    ggsave(outfile, plot = plot_obj, width = width, height = height, dpi = dpi)
  }
  
  message("Saved: ", outfile)
}



# ======== Cell Type Distribution Visualization ========

plot_cell_dist <- function(SeuratObject, group_var = "celltype_major",
                           main_title = "Cell Type Distribution",
                           col = "steelblue",
                           cex_names = 0.8,    # x-axis labels
                           cex_axis = 0.8,     # y-axis labels
                           ylim = NULL,
                           return_plot = TRUE) {  # NEW option
  if (!group_var %in% colnames(SeuratObject@meta.data)) {
    stop(paste("Grouping variable", group_var, "not found in meta.data"))
  }
  
  # counts
  counts <- table(SeuratObject@meta.data[[group_var]])
  df <- data.frame(CellType = names(counts), Count = as.numeric(counts))
  print(counts)
  
  # base R version (still works, side-effect only)
  old_mar <- par("mar")
  par(mar = c(8, old_mar[2], old_mar[3], old_mar[4]))
  if (is.null(ylim)) {
    ylim <- c(0, max(counts) * 1.1)
  }
  barplot(counts,
          main = main_title,
          ylab = "Number of Cells",
          col = col,
          las = 2,
          cex.names = cex_names,
          cex.axis = cex_axis, 
          ylim = ylim)
  par(mar = old_mar)
  
  # ggplot2 object (returned if requested)
  library(ggplot2)
  p <- ggplot(df, aes(x = CellType, y = Count, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = cex_names*12),
          axis.text.y = element_text(size = cex_axis*12),
          legend.position = "none") +
    labs(title = main_title, x = NULL, y = "Number of Cells") +
    scale_fill_manual(values = rep(col, length(unique(df$CellType))))
  
  if (return_plot) {
    return(p)   # ggplot object
  } else {
    invisible(NULL)
  }
}




# ======== Bimodal Threshold Determination ========

find_peak_threshold <- function(scores, bandwidth = .05, select = "last") {
  # Kernel density estimation
  dens <- density(scores, bw = bandwidth)
  
  # CORRECT: Find local maxima (PEAKS) where derivative changes + to -
  # diff(sign(diff(y))) == -2 indicates a local maximum
  peaks <- which(diff(sign(diff(dens$y))) == -2) + 1
  
  cat("Peaks detected at indices:", peaks, "\n")
  cat("Peak values (x):", round(dens$x[peaks], 3), "\n")
  cat("Peak densities (y):", round(dens$y[peaks], 3), "\n")
  
  if (length(peaks) >= 2) {
    # Bimodal distribution detected
    if (select == "last") {
      # Use the LAST peak (second mode in bimodal distribution)
      threshold <- dens$x[peaks[length(peaks)]]
      cat("Bimodal distribution: using last peak as threshold\n")
      threshold_l <- threshold - 1 * sd(scores[scores>threshold]) # lower limit
      threshold_h <- threshold + 1 * sd(scores[scores>threshold]) # high limit
    } else if (select == "mean") {
      # Use mean of peaks
      threshold <- mean(dens$x[peaks])
      cat("Bimodal distribution: using mean of peaks as threshold\n")
      threshold_l <- threshold - 1 * sd(scores[scores>threshold]) # lower limit
      threshold_h <- threshold + 1 * sd(scores[scores>threshold]) # high limit
    }
  } else if (length(peaks) == 1) {
    # Unimodal distribution
    threshold <- quantile(scores, 0.90)
    cat("Unimodal distribution: using 90th percentile\n")
  } else {
    # No clear peaks
    threshold <- quantile(scores, 0.90)
    cat("No peaks detected: using 90th percentile\n")
  }
  
  return(list(
    threshold = threshold,
    threshold_low = threshold_l,
    threshold_high = threshold_h,
    peaks = peaks,
    peak_x = dens$x[peaks],
    peak_y = dens$y[peaks],
    density = dens
  ))
}




# ======== Peak Threshold Visualization ========

visualize_peak_detection <- function(scores, result, title = "Peak Detection", xname ="Score") {
  
  df <- data.frame(
    x = result$density$x,
    y = result$density$y
  )
  
  τ      <- result$threshold
  τ_low  <- result$threshold_low
  τ_high <- result$threshold_high
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(color = "black", size = 1.1)
  
  # --------------------------------------------------------------
  # 1) Stable region (x < τ_low)
  # --------------------------------------------------------------
  if (!is.na(τ_low)) {
    p <- p + geom_area(
      data = df[df$x < τ_low, ],
      aes(x = x, y = y),
      fill = "green", alpha = 0.20
    )
  }
  
  # --------------------------------------------------------------
  # 2) Transition region (τ_low ≤ x < τ_high)
  # --------------------------------------------------------------
  if (!is.na(τ_low) & !is.na(τ_high)) {
    p <- p + geom_area(
      data = df[df$x >= τ_low & df$x < τ_high, ],
      aes(x = x, y = y),
      fill = "orange", alpha = 0.20
    )
  }
  
  # --------------------------------------------------------------
  # 3) Divergent region (x ≥ τ_high or fallback region)
  # --------------------------------------------------------------
  if (!is.na(τ_high)) {
    p <- p + geom_area(
      data = df[df$x >= τ_high, ],
      aes(x = x, y = y),
      fill = "red", alpha = 0.20
    )
  } else {
    p <- p + geom_area(
      data = df[df$x >= τ, ],
      aes(x = x, y = y),
      fill = "red", alpha = 0.20
    )
  }
  
  # --------------------------------------------------------------
  # Stable Proxy & Divergent labels
  # --------------------------------------------------------------
  p <- p +
    annotate("text",
             x = mean(df$x[df$x < τ]),
             y = max(df$y) * 0.7,
             label = "Stable Proxy",
             color = "darkgreen", size = 4) +
    
    annotate("text",
             x = mean(df$x[df$x >= τ]),
             y = max(df$y) * 0.7,
             label = "Divergent",
             color = "darkred", size = 4)
  
  # --------------------------------------------------------------
  # Threshold lines
  # --------------------------------------------------------------
  p <- p +
    geom_vline(xintercept = τ,
               color = "red", linetype = "dashed", size = 1.3) +
    annotate("text",
             x = τ,
             y = max(df$y) * 0.95,
             label = paste0("τ = ", round(τ, 3)),
             color = "red", hjust = -0.1, fontface = "bold")
  
  if (!is.na(τ_low)) {
    p <- p +
      geom_vline(xintercept = τ_low,
                 color = "darkgreen", linetype = "dotted", size = 1) 
    # + annotate("text",
    #            x = τ_low,
    #            y = max(df$y) * 0.85,
    #            label = paste0("τ_low = ", round(τ_low, 3)),
    #            color = "darkgreen", hjust = 1.1, size = 4)
  }
  
  if (!is.na(τ_high)) {
    p <- p +
      geom_vline(xintercept = τ_high,
                 color = "darkred", linetype = "dotted", size = 1) 
    # + annotate("text",
    #            x = τ_high,
    #            y = max(df$y) * 0.85,
    #            label = paste0("τ_high = ", round(τ_high, 3)),
    #            color = "darkred", hjust = -0.1, size = 4)
  }
  
  # --------------------------------------------------------------
  # Labels + Theme
  # --------------------------------------------------------------
  p <- p +
    labs(
      title = title,
      x = xname,
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title  = element_text(face = "bold", size = 14),
      axis.title  = element_text(size = 16),
      axis.text   = element_text(size = 14)
    )
  
  return(p)
}





# ======== Marker Gene Expression Visualization ========

plot_marker_expression <- function(logCPM_all, metaInfo, marker_genes, 
                                   title = "Marker gene expression") {

  # Check inputs
  if (!all(marker_genes %in% rownames(logCPM_all))) {
    warning("Some marker genes not found in expression matrix; ignoring missing ones.")
    marker_genes <- marker_genes[marker_genes %in% rownames(logCPM_all)]
  }
  
  # Prepare expression data frame
  expr_df <- data.frame(
    gene     = rep(rownames(logCPM_all), ncol(logCPM_all)),
    expr     = as.vector(logCPM_all),
    sample   = rep(colnames(logCPM_all), each = nrow(logCPM_all)),
    platform = rep(metaInfo$platform, each = nrow(logCPM_all))
  )
  
  # Filter only marker genes
  expr_markers <- expr_df %>% filter(gene %in% marker_genes)
  
  # Plot
  p <- ggplot(expr_markers, aes(x = platform, y = expr, fill = platform)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.11, outlier.size = 0.8) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("bulk" = "#1f78b4", "pb" = "#33a02c")) +
    labs(title = title, x = "", y = "logCPM") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y   = element_text(size = 14),
      axis.title.x  = element_text(size = 16),
      axis.title.y  = element_text(size = 16),
      plot.title    = element_text(size = 16, hjust = 0.5),
      panel.grid    = element_blank(),
      legend.position = "none"
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    )
  return(p)
}




