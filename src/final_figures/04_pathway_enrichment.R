#!/usr/bin/env Rscript
# GO Biological Process over-representation analysis on pseudobulk DEGs.
# Uses clusterProfiler + org.Mm.eg.db (no internet required for GO).
#
# Outputs (src/final_figures/output/figures/):
#   pathway_dotplot_heatvsctrl.pdf   – per-cluster GO enrichment, heat vs ctrl
#   pathway_dotplot_agedvsadult.pdf  – per-cluster GO enrichment, aged vs adult
#   pathway_barplot_combined.pdf     – top GO terms across all clusters combined
#
# Saves RDS (for downstream network script):
#   src/final_figures/output/enrich_heatvsctrl.rds
#   src/final_figures/output/enrich_agedvsadult.rds

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(RColorBrewer)
})

BASE    <- "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  <- file.path(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR <- file.path(BASE, "src/final_figures/output")
OUTDIR  <- file.path(BASE, "src/final_figures/output/figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATADIR, recursive = TRUE, showWarnings = FALSE)

N_CLUSTERS <- 17
P_GENE     <- 0.05      # gene-level significance threshold
FC_GENE    <- 0.5       # |logFC| threshold for gene selection
P_ENRICH   <- 0.05      # enrichment FDR threshold
TOP_TERMS  <- 5         # top GO terms per cluster for per-cluster dotplot
TOP_GLOBAL <- 30        # terms for combined barplot

# ── Helper: load significant DEGs for one comparison ─────────────────────────
load_sig_genes <- function(comparison, direction = "both") {
  map_dfr(0:(N_CLUSTERS - 1), function(cl) {
    f <- file.path(DEG_PB, paste0("DEG_pb_", comparison, "_cluster", cl, ".csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f, row.names = 1) %>%
      mutate(gene = rownames(.), cluster = as.character(cl))
    sig <- df %>% filter(adj.P.Val < P_GENE, abs(logFC) > FC_GENE)
    if (direction == "up")   sig <- filter(sig, logFC > 0)
    if (direction == "down") sig <- filter(sig, logFC < 0)
    sig
  })
}

# ── Helper: run enrichGO for a gene list ──────────────────────────────────────
run_enrichGO <- function(genes, universe = NULL) {
  if (length(genes) < 5) return(NULL)
  tryCatch(
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = P_ENRICH,
      qvalueCutoff  = 0.2,
      readable      = FALSE
    ),
    error = function(e) { message("enrichGO error: ", e$message); NULL }
  )
}

# ── Helper: shorten GO term descriptions ──────────────────────────────────────
shorten_term <- function(x, max_char = 45) {
  ifelse(nchar(x) > max_char, paste0(substr(x, 1, max_char - 1), "…"), x)
}

# ── Per-cluster enrichment + faceted dot plot ─────────────────────────────────
per_cluster_enrich_pdf <- function(comparison, outfile, rds_file) {
  message("\nEnrichment: ", comparison)
  sig_df <- load_sig_genes(comparison)
  if (nrow(sig_df) == 0) { message("No significant genes."); return() }

  # Universe = all measured genes in the panel
  panel_file <- file.path(DATADIR, "panel_genes.csv")
  universe   <- if (file.exists(panel_file)) read.csv(panel_file)$gene else NULL

  enrich_list <- sig_df %>%
    split(.$cluster) %>%
    imap(function(d, cl) {
      genes <- unique(d$gene)
      res   <- run_enrichGO(genes, universe)
      if (!is.null(res)) message("  Cluster ", cl, ": ", nrow(res@result), " terms")
      list(cluster = cl, result = res, genes = genes)
    })

  saveRDS(enrich_list, rds_file)
  message("Saved RDS: ", rds_file)

  # Build tidy data frame of top terms per cluster
  top_df <- map_dfr(enrich_list, function(item) {
    res <- item$result
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res@result %>%
      filter(p.adjust < P_ENRICH) %>%
      slice_min(order_by = p.adjust, n = TOP_TERMS, with_ties = FALSE) %>%
      mutate(cluster = item$cluster)
  })

  if (nrow(top_df) == 0) {
    message("No enriched terms found for any cluster.")
    return()
  }

  top_df <- top_df %>%
    mutate(
      Description = shorten_term(Description),
      GeneRatio_num = map_dbl(GeneRatio, ~ {
        parts <- str_split(.x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      cluster = factor(cluster, levels = as.character(sort(as.integer(unique(cluster)))))
    )

  p <- ggplot(top_df,
              aes(x = cluster, y = Description,
                  size = GeneRatio_num, color = -log10(p.adjust))) +
    geom_point() +
    scale_size_continuous(name = "Gene Ratio", range = c(1.5, 5)) +
    scale_color_gradientn(
      name   = "-log10(FDR)",
      colors = brewer.pal(9, "YlOrRd")[3:9]
    ) +
    labs(
      title = paste0("GO BP Enrichment – ",
                     str_replace_all(comparison, "_", " ") %>% str_to_title()),
      x = "Leiden Cluster", y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.y     = element_text(size = 7),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title      = element_text(size = 11, face = "bold")
    )

  h <- max(6, length(unique(top_df$Description)) * 0.28 + 2)
  pdf(outfile, width = max(10, N_CLUSTERS * 0.6 + 2), height = h)
  print(p)
  invisible(dev.off())
  message("Saved: ", outfile)
}

# ── Combined (all-cluster) bar plot ───────────────────────────────────────────
combined_barplot_pdf <- function(comparison1, comparison2, outfile) {
  rds1 <- file.path(DATADIR, paste0("enrich_",
                                     str_remove(comparison1, "comp\\d+_"), ".rds"))
  rds2 <- file.path(DATADIR, paste0("enrich_",
                                     str_remove(comparison2, "comp\\d+_"), ".rds"))

  load_top_terms <- function(rds_path, label) {
    if (!file.exists(rds_path)) return(NULL)
    enrich_list <- readRDS(rds_path)
    map_dfr(enrich_list, function(item) {
      res <- item$result
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      res@result %>%
        filter(p.adjust < P_ENRICH) %>%
        mutate(comparison = label, cluster = item$cluster)
    })
  }

  df1 <- load_top_terms(rds1, "Heat vs Ctrl")
  df2 <- load_top_terms(rds2, "Aged vs Adult")

  combined <- bind_rows(df1, df2)
  if (is.null(combined) || nrow(combined) == 0) {
    message("No data for combined barplot."); return()
  }

  top_terms <- combined %>%
    group_by(ID, Description, comparison) %>%
    summarise(min_padj = min(p.adjust), .groups = "drop") %>%
    group_by(comparison) %>%
    slice_min(order_by = min_padj, n = TOP_GLOBAL / 2, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(Description)

  plot_df <- combined %>%
    filter(Description %in% top_terms$Description) %>%
    group_by(Description, comparison) %>%
    summarise(
      neg_log10_padj = mean(-log10(p.adjust)),
      n_clusters     = n_distinct(cluster),
      .groups        = "drop"
    ) %>%
    mutate(Description = shorten_term(Description, 55))

  # Order by total significance
  term_order <- plot_df %>%
    group_by(Description) %>%
    summarise(score = sum(neg_log10_padj)) %>%
    arrange(desc(score)) %>%
    pull(Description)
  plot_df$Description <- factor(plot_df$Description, levels = term_order)

  p <- ggplot(plot_df,
              aes(x = neg_log10_padj, y = Description,
                  fill = comparison, size = n_clusters)) +
    geom_point(shape = 21, alpha = 0.85) +
    scale_fill_manual(
      values = c("Heat vs Ctrl" = "#FC8D62", "Aged vs Adult" = "#8DA0CB"),
      name   = "Comparison"
    ) +
    scale_size_continuous(name = "# Clusters", range = c(2, 7)) +
    labs(
      title = "Top GO BP Terms – Main Comparisons",
      x = expression(-log[10](FDR)), y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
      plot.title       = element_text(size = 12, face = "bold")
    )

  h <- max(7, nrow(top_terms) * 0.3 + 2)
  pdf(outfile, width = 10, height = h)
  print(p)
  invisible(dev.off())
  message("Saved: ", outfile)
}

# ── Run ────────────────────────────────────────────────────────────────────────
per_cluster_enrich_pdf(
  "comp1_heatvsctrl",
  file.path(OUTDIR, "pathway_dotplot_heatvsctrl.pdf"),
  file.path(DATADIR, "enrich_heatvsctrl.rds")
)

per_cluster_enrich_pdf(
  "comp1_agedvsadult",
  file.path(OUTDIR, "pathway_dotplot_agedvsadult.pdf"),
  file.path(DATADIR, "enrich_agedvsadult.rds")
)

combined_barplot_pdf(
  "comp1_heatvsctrl",
  "comp1_agedvsadult",
  file.path(OUTDIR, "pathway_barplot_combined.pdf")
)
