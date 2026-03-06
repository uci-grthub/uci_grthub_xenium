#!/usr/bin/env Rscript
# Heatmaps of top DEGs and curated pathway gene sets.
#
# Outputs (src/final_figures/output/figures/):
#   heatmap_logFC_heatvsctrl.pdf      – logFC matrix, heat vs ctrl (needs no h5ad)
#   heatmap_logFC_agedvsadult.pdf     – logFC matrix, aged vs adult
#   heatmap_expression_curated.pdf    – mean expression of curated pathway genes
#                                        (requires 00_extract_h5ad_data.py output)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(RColorBrewer)
  library(grid)
})

BASE    <- "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  <- file.path(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR <- file.path(BASE, "src/final_figures/output")
OUTDIR  <- file.path(BASE, "src/final_figures/output/figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

N_CLUSTERS <- 17   # Leiden clusters 0–16
TOP_N      <- 6    # top genes per cluster to keep in heatmap

# ── Helper: load all cluster files for one pseudobulk comparison ──────────────
load_pb_comp <- function(comparison) {
  map_dfr(0:(N_CLUSTERS - 1), function(cl) {
    f <- file.path(DEG_PB, paste0("DEG_pb_", comparison, "_cluster", cl, ".csv"))
    if (!file.exists(f)) return(NULL)
    read.csv(f, row.names = 1) %>%
      mutate(gene = rownames(.), cluster = as.character(cl))
  })
}

# ── Helper: build gene × cluster logFC matrix ─────────────────────────────────
build_logfc_matrix <- function(df, p_thresh = 0.05, fc_thresh = 0.25,
                                top_n = TOP_N) {
  # Select genes significant in ≥1 cluster, take top_n per cluster by |logFC|
  sig <- df %>%
    filter(adj.P.Val < p_thresh, abs(logFC) > fc_thresh)

  top_genes <- sig %>%
    group_by(cluster) %>%
    slice_max(order_by = abs(logFC), n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(gene) %>%
    pull(gene)

  if (length(top_genes) == 0) {
    message("No significant genes found."); return(NULL)
  }

  mat <- df %>%
    filter(gene %in% top_genes) %>%
    select(gene, cluster, logFC) %>%
    pivot_wider(names_from = cluster, values_from = logFC, values_fill = 0) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  # Order columns numerically
  col_order <- as.character(sort(as.integer(colnames(mat))))
  mat[, col_order, drop = FALSE]
}

# ── Helper: draw and save heatmap ─────────────────────────────────────────────
save_heatmap <- function(mat, title, outfile,
                         col_annot_df = NULL,
                         row_annot_df = NULL,
                         row_split    = NULL) {
  if (is.null(mat)) { message("Skipping (null matrix): ", outfile); return() }

  fc_lim <- max(abs(mat), na.rm = TRUE)
  fc_lim <- min(fc_lim, 3)
  col_fun <- colorRamp2(
    c(-fc_lim, 0, fc_lim),
    c("#4575B4", "white", "#D73027")
  )

  col_ha <- NULL
  if (!is.null(col_annot_df)) {
    age_col  <- c(adult = "#2166AC", aged = "#B2182B")
    cond_col <- c(control = "#66C2A5", heatstress = "#FC8D62")
    ha_list  <- list()
    if ("age" %in% names(col_annot_df))
      ha_list[["age"]] <- col_annot_df$age
    if ("condition" %in% names(col_annot_df))
      ha_list[["condition"]] <- col_annot_df$condition

    col_ha <- HeatmapAnnotation(
      df     = col_annot_df,
      col    = list(age = age_col, condition = cond_col),
      show_annotation_name = TRUE,
      annotation_name_gp   = gpar(fontsize = 8)
    )
  }

  row_ha <- NULL
  if (!is.null(row_annot_df)) {
    pathway_colors <- setNames(
      colorRampPalette(brewer.pal(8, "Set2"))(length(unique(row_annot_df$pathway))),
      unique(row_annot_df$pathway)
    )
    row_ha <- rowAnnotation(
      pathway = row_annot_df$pathway,
      col     = list(pathway = pathway_colors),
      show_annotation_name = FALSE,
      annotation_width     = unit(3, "mm")
    )
  }

  ht <- Heatmap(
    mat,
    name                 = "logFC",
    col                  = col_fun,
    column_title         = title,
    column_title_gp      = gpar(fontsize = 11, fontface = "bold"),
    cluster_rows         = TRUE,
    cluster_columns      = FALSE,
    row_names_gp         = gpar(fontsize = 7),
    column_names_gp      = gpar(fontsize = 8),
    top_annotation       = col_ha,
    left_annotation      = row_ha,
    row_split            = row_split,
    row_gap              = unit(1, "mm"),
    border               = TRUE,
    heatmap_legend_param = list(
      title    = "log2 FC",
      at       = c(-fc_lim, 0, fc_lim),
      labels   = c(sprintf("-%.1f", fc_lim), "0", sprintf("+%.1f", fc_lim)),
      direction = "horizontal"
    )
  )

  h  <- max(5, nrow(mat) * 0.22 + 2)
  w  <- max(6, ncol(mat) * 0.5  + 2)
  pdf(outfile, width = w, height = h)
  draw(ht, heatmap_legend_side = "bottom")
  invisible(dev.off())
  message("Saved: ", outfile)
}

# ══ Heatmap 1: pseudobulk heat vs ctrl logFC ══════════════════════════════════
message("\n── Heatmap: heat vs ctrl ──")
df_heat <- load_pb_comp("comp1_heatvsctrl")
mat_heat <- build_logfc_matrix(df_heat)

col_anno <- data.frame(
  condition = rep("heatstress", ncol(mat_heat)),
  row.names = colnames(mat_heat)
)

save_heatmap(
  mat_heat,
  title   = "Top DEGs – Heat Stress vs Control (pseudobulk, all clusters)",
  outfile = file.path(OUTDIR, "heatmap_logFC_heatvsctrl.pdf"),
  col_annot_df = NULL
)

# ══ Heatmap 2: pseudobulk aged vs adult logFC ═════════════════════════════════
message("\n── Heatmap: aged vs adult ──")
df_age  <- load_pb_comp("comp1_agedvsadult")
mat_age <- build_logfc_matrix(df_age)

save_heatmap(
  mat_age,
  title   = "Top DEGs – Aged vs Adult (pseudobulk, all clusters)",
  outfile = file.path(OUTDIR, "heatmap_logFC_agedvsadult.pdf")
)

# ══ Heatmap 3: curated gene expression (requires 00_extract_h5ad_data.py) ════
expr_file    <- file.path(DATADIR, "cluster_mean_expr.csv")
curated_file <- file.path(DATADIR, "curated_gene_sets.csv")

if (!file.exists(expr_file)) {
  message("\nSkipping curated expression heatmap – run 00_extract_h5ad_data.py first.")
} else {
  message("\n── Heatmap: curated pathway gene expression ──")
  expr_long <- read_csv(expr_file, show_col_types = FALSE)
  curated   <- read_csv(curated_file, show_col_types = FALSE)

  curated_genes <- unique(curated$gene)
  gene_cols <- intersect(curated_genes, colnames(expr_long))

  if (length(gene_cols) == 0) {
    message("None of the curated genes found in expression matrix.")
  } else {
    # Build column label: Cluster_age_condition
    expr_long <- expr_long %>%
      mutate(col_label = paste0("C", leiden_harmony, "_",
                                str_sub(age, 1, 2), "_",
                                str_sub(condition, 1, 4)))

    mat_expr <- expr_long %>%
      select(col_label, all_of(gene_cols)) %>%
      tibble::column_to_rownames("col_label") %>%
      t() %>%
      as.matrix()

    # Z-score per gene
    mat_z <- t(scale(t(mat_expr)))
    mat_z[is.nan(mat_z)] <- 0

    # Column annotation
    col_meta <- expr_long %>%
      select(col_label, age, condition) %>%
      tibble::column_to_rownames("col_label")
    col_meta <- col_meta[colnames(mat_z), , drop = FALSE]

    # Row annotation
    row_annot <- curated %>%
      filter(gene %in% rownames(mat_z)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      tibble::column_to_rownames("gene")
    row_annot <- row_annot[rownames(mat_z), , drop = FALSE]
    row_split <- row_annot$pathway

    z_lim   <- 2.5
    col_fun <- colorRamp2(c(-z_lim, 0, z_lim),
                          c("#4575B4", "#FFFFBF", "#D73027"))

    age_col  <- c(adult = "#2166AC", aged = "#B2182B")
    cond_col <- c(control = "#66C2A5", heatstress = "#FC8D62")

    top_ha <- HeatmapAnnotation(
      age       = col_meta$age,
      condition = col_meta$condition,
      col       = list(age = age_col, condition = cond_col),
      annotation_name_gp = gpar(fontsize = 8)
    )

    pw_colors <- setNames(
      colorRampPalette(brewer.pal(8, "Set2"))(length(unique(row_split))),
      unique(row_split)
    )
    left_ha <- rowAnnotation(
      pathway = row_split,
      col     = list(pathway = pw_colors),
      show_annotation_name = FALSE,
      annotation_width     = unit(3, "mm")
    )

    ht_expr <- Heatmap(
      mat_z,
      name                 = "z-score",
      col                  = col_fun,
      column_title         = "Curated Pathway Gene Sets – Mean Expression",
      column_title_gp      = gpar(fontsize = 11, fontface = "bold"),
      cluster_rows         = FALSE,
      cluster_columns      = FALSE,
      row_split            = row_split,
      row_gap              = unit(2, "mm"),
      row_names_gp         = gpar(fontsize = 7),
      column_names_gp      = gpar(fontsize = 6),
      top_annotation       = top_ha,
      left_annotation      = left_ha,
      border               = TRUE,
      heatmap_legend_param = list(
        title     = "z-score",
        direction = "horizontal"
      )
    )

    h  <- max(8, length(gene_cols) * 0.22 + 3)
    w  <- max(10, ncol(mat_z) * 0.18 + 3)
    pdf(file.path(OUTDIR, "heatmap_expression_curated.pdf"),
        width = w, height = h)
    draw(ht_expr, heatmap_legend_side = "bottom")
    invisible(dev.off())
    message("Saved: heatmap_expression_curated.pdf")
  }
}
