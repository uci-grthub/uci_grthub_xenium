#!/usr/bin/env Rscript
# Volcano plots for pseudobulk and cluster-wise (Wilcoxon) DEGs.
#
# Outputs (src/final_figures/output/figures/):
#   volcano_pb_<comparison>.pdf   – pseudobulk, one panel per Leiden cluster
#   volcano_std_<file>.pdf        – cluster-wise Wilcoxon, one panel per cluster

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(readxl)
  library(stringr)
  library(purrr)
  library(patchwork)
})

BASE   <- "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB <- file.path(BASE, "results/to_client/all/degs/pseudobulk")
DEG_ST <- file.path(BASE, "results/to_client/all/degs/standard")
OUTDIR <- file.path(BASE, "src/final_figures/output/figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ── Theme ──────────────────────────────────────────────────────────────────────
theme_volcano <- function() {
  theme_bw(base_size = 8) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.25, color = "grey90"),
      strip.background  = element_rect(fill = "grey20"),
      strip.text        = element_text(color = "white", size = 7),
      legend.position   = "none",
      plot.title        = element_text(size = 8, face = "bold"),
      axis.title        = element_text(size = 7)
    )
}

COL_UP <- "#C0392B"
COL_DN <- "#2980B9"
COL_NS <- "grey75"

# ── Core volcano helper ────────────────────────────────────────────────────────
make_volcano <- function(df, title,
                         fc_col    = "logFC",
                         pval_col  = "adj.P.Val",
                         gene_col  = "gene",
                         fc_thresh = 0.5,
                         p_thresh  = 0.05,
                         top_n     = 6) {
  df <- df %>%
    mutate(
      fc   = as.numeric(.data[[fc_col]]),
      pval = as.numeric(.data[[pval_col]]),
      nlp  = -log10(pmax(pval, 1e-300)),
      sig  = case_when(
        pval < p_thresh & fc >  fc_thresh ~ "Up",
        pval < p_thresh & fc < -fc_thresh ~ "Down",
        TRUE                              ~ "NS"
      )
    )

  label_df <- df %>%
    filter(sig != "NS") %>%
    slice_max(order_by = nlp, n = top_n)

  ggplot(df, aes(fc, nlp, color = sig)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = -log10(p_thresh),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_label_repel(
      data        = label_df,
      aes(label   = .data[[gene_col]]),
      size        = 2.0,
      max.overlaps = 15,
      box.padding = 0.25,
      label.size  = 0.15,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(Up = COL_UP, Down = COL_DN, NS = COL_NS)) +
    labs(title = title, x = "log2 FC", y = expression(-log[10](adj.p))) +
    theme_volcano()
}

# ── Pseudobulk multi-panel volcano ────────────────────────────────────────────
pb_volcano_pdf <- function(comparison, outfile) {
  pattern  <- paste0("DEG_pb_", comparison, "_cluster")
  files    <- sort(list.files(DEG_PB, pattern = pattern, full.names = TRUE))
  if (length(files) == 0) {
    message("No files for: ", comparison); return(invisible(NULL))
  }

  all_data <- map_dfr(files, function(f) {
    cl_num <- str_extract(basename(f), "(?<=cluster)\\d+$") %>% as.integer()
    df <- read.csv(f, row.names = 1) %>%
      mutate(gene = rownames(.), cluster_n = cl_num,
             cluster = paste0("Cluster ", cl_num))
    df
  }) %>%
    mutate(cluster = factor(cluster,
                            levels = paste0("Cluster ", sort(unique(cluster_n)))))

  plots <- all_data %>%
    split(.$cluster) %>%
    imap(function(d, cl) make_volcano(d, cl, top_n = 5))

  ncols    <- min(6, length(plots))
  nrows    <- ceiling(length(plots) / ncols)
  title_str <- comparison %>%
    str_replace_all("_", " ") %>%
    str_replace("comp\\d+ ", "") %>%
    str_to_title()

  combined <- wrap_plots(plots, ncol = ncols) +
    plot_annotation(
      title = paste("Pseudobulk DEGs –", title_str),
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              hjust = 0.5))
    )

  pdf(outfile, width = 3.5 * ncols, height = 3.5 * nrows)
  print(combined)
  invisible(dev.off())
  message("Saved: ", outfile)
}

pb_comps <- c(
  "comp1_heatvsctrl",
  "comp1_agedvsadult",
  "comp2_adultonly_heatvsctrl",
  "comp2_agedonly_heatvsctrl"
)

for (comp in pb_comps) {
  outfile <- file.path(OUTDIR, paste0("volcano_pb_", comp, ".pdf"))
  pb_volcano_pdf(comp, outfile)
}

# ── Standard (Wilcoxon) cluster-wise volcano ───────────────────────────────────
std_volcano_pdf <- function(excel_file, outfile) {
  if (!file.exists(excel_file)) {
    message("Not found: ", excel_file); return(invisible(NULL))
  }
  sheets <- excel_sheets(excel_file)

  all_data <- map_dfr(sheets, function(sh) {
    df <- read_excel(excel_file, sheet = sh, col_types = "text")
    if (!"names" %in% names(df)) return(NULL)
    df %>%
      transmute(
        gene      = names,
        logFC     = as.numeric(logfoldchanges),
        adj.P.Val = as.numeric(pvals_adj),
        score     = as.numeric(scores),
        cluster   = paste0("Cluster ", sh)
      )
  }) %>%
    filter(!is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(cluster = factor(cluster,
                            levels = paste0("Cluster ", sort(as.integer(sheets)))))

  plots <- all_data %>%
    split(.$cluster) %>%
    imap(function(d, cl) make_volcano(d, cl, fc_thresh = 0.25, top_n = 5))

  ncols    <- min(6, length(plots))
  nrows    <- ceiling(length(plots) / ncols)
  title_str <- tools::file_path_sans_ext(basename(excel_file)) %>%
    str_replace_all("_ranked.*", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title()

  combined <- wrap_plots(plots, ncol = ncols) +
    plot_annotation(
      title = paste("Wilcoxon DEGs –", title_str),
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              hjust = 0.5))
    )

  pdf(outfile, width = 3.5 * ncols, height = 3.5 * nrows)
  print(combined)
  invisible(dev.off())
  message("Saved: ", outfile)
}

std_files <- c(
  "adult_heatstress_v_control_ranked_genes_per_cluster.xlsx",
  "aged_heatstress_v_control_ranked_genes_per_cluster.xlsx",
  "all_adult_v_aged_ranked_genes_per_cluster.xlsx",
  "all_control_v_heatstress_ranked_genes_per_cluster.xlsx"
)

for (f in std_files) {
  fpath   <- file.path(DEG_ST, f)
  outname <- paste0("volcano_std_",
                    tools::file_path_sans_ext(f), ".pdf")
  std_volcano_pdf(fpath, file.path(OUTDIR, outname))
}
