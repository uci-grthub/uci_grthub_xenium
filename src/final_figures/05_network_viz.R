#!/usr/bin/env Rscript
# Network visualizations for biologically relevant gene modules.
# Requires enrichment RDS files from 04_pathway_enrichment.R.
#
# Outputs (src/final_figures/output/figures/):
#   network_cnet_heatvsctrl.pdf    – concept network: top GO terms & hub genes
#   network_cnet_agedvsadult.pdf   – same, for aged vs adult comparison
#   network_coexpr_heatstress.pdf  – gene co-expression network (heat stress
#                                    genes, built from pseudobulk logFC profiles)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(tidyr)
})

BASE    <- "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
DEG_PB  <- file.path(BASE, "results/to_client/all/degs/pseudobulk")
DATADIR <- file.path(BASE, "src/final_figures/output")
OUTDIR  <- file.path(BASE, "src/final_figures/output/figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

P_GENE    <- 0.05
FC_GENE   <- 0.5
N_TERMS   <- 10     # top GO terms to display per network
N_GENES   <- 5      # max genes shown per term
CORR_CUTOFF <- 0.5  # absolute Pearson r for co-expression edges

shorten_term <- function(x, n = 40) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
}

# ── Helper: build concept (cnet) graph from enrich RDS ───────────────────────
build_cnet_graph <- function(rds_file, comparison_label,
                              n_terms = N_TERMS, n_genes = N_GENES) {
  if (!file.exists(rds_file)) {
    message("RDS not found – run 04_pathway_enrichment.R first: ", rds_file)
    return(NULL)
  }
  enrich_list <- readRDS(rds_file)

  # Collect top terms across clusters (ranked by min adjusted p across clusters)
  all_terms <- map_dfr(enrich_list, function(item) {
    res <- item$result
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res@result %>%
      filter(p.adjust < 0.05) %>%
      mutate(cluster = item$cluster)
  })

  if (nrow(all_terms) == 0) {
    message("No enriched terms found in ", rds_file); return(NULL)
  }

  top_term_ids <- all_terms %>%
    group_by(ID, Description) %>%
    summarise(best_padj = min(p.adjust), .groups = "drop") %>%
    slice_min(order_by = best_padj, n = n_terms, with_ties = FALSE) %>%
    pull(ID)

  term_gene_df <- all_terms %>%
    filter(ID %in% top_term_ids) %>%
    select(ID, Description, geneID, p.adjust) %>%
    mutate(Description = shorten_term(Description)) %>%
    group_by(ID) %>%
    slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(genes = str_split(geneID, "/")) %>%
    unnest(genes) %>%
    group_by(ID) %>%
    slice_head(n = n_genes) %>%   # limit genes per term for clarity
    ungroup()

  if (nrow(term_gene_df) == 0) return(NULL)

  # Nodes
  term_nodes <- term_gene_df %>%
    distinct(ID, Description, p.adjust) %>%
    transmute(name = Description, type = "pathway",
              size = -log10(p.adjust + 1e-10), color = p.adjust)

  gene_nodes <- term_gene_df %>%
    distinct(genes) %>%
    transmute(name = genes, type = "gene",
              size = 3, color = NA_real_)

  nodes <- bind_rows(term_nodes, gene_nodes)

  # Edges
  edges <- term_gene_df %>%
    transmute(from = Description, to = genes)

  # Build igraph
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  V(g)$type  <- V(g)$type
  V(g)$size  <- V(g)$size
  list(graph = g, label = comparison_label)
}

# ── Helper: draw and save concept network ─────────────────────────────────────
save_cnet_pdf <- function(cnet, outfile) {
  if (is.null(cnet)) return()
  g <- cnet$graph

  # Node aesthetics
  type_colors <- c(pathway = "#E74C3C", gene = "#3498DB")
  type_sizes  <- c(pathway = 8,         gene = 4)

  set.seed(42)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "grey60", linewidth = 0.4) +
    geom_node_point(
      aes(size = type, color = type, fill = type),
      shape = 21, stroke = 0.4
    ) +
    geom_node_label(
      aes(label = name, size = type),
      repel          = TRUE,
      label.padding  = unit(0.12, "lines"),
      label.size     = 0.15,
      max.overlaps   = 20,
      show.legend    = FALSE
    ) +
    scale_color_manual(values = type_colors, name = "Node type",
                       guide  = guide_legend(override.aes = list(size = 4))) +
    scale_fill_manual(values  = type_colors, name = "Node type",
                      guide   = guide_legend(override.aes = list(size = 4))) +
    scale_size_manual(values  = type_sizes, name = "Node type") +
    labs(title = paste("Concept Network –", cnet$label)) +
    theme_graph(base_family = "sans", base_size = 9) +
    theme(
      plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  pdf(outfile, width = 12, height = 10)
  print(p)
  invisible(dev.off())
  message("Saved: ", outfile)
}

# ── Co-expression network for heat stress gene module ─────────────────────────
build_coexpr_pdf <- function(comparison, outfile) {
  # Use logFC profiles across clusters as a gene × cluster expression-proxy
  files <- sort(list.files(
    DEG_PB, pattern = paste0("DEG_pb_", comparison, "_cluster"),
    full.names = TRUE
  ))
  if (length(files) == 0) {
    message("No files for co-expression: ", comparison); return()
  }

  logfc_list <- map(files, function(f) {
    cl  <- str_extract(basename(f), "(?<=cluster)\\d+")
    df  <- read.csv(f, row.names = 1) %>%
      mutate(gene = rownames(.))
    setNames(df$logFC, df$gene)
  })

  all_genes <- unique(unlist(lapply(logfc_list, names)))
  mat <- sapply(logfc_list, function(lv) {
    lv[match(all_genes, names(lv))]
  })
  rownames(mat) <- all_genes
  mat[is.na(mat)] <- 0

  # Keep genes with |mean logFC| > 0.3 (somewhat differentially expressed)
  gene_var <- apply(mat, 1, function(x) mean(abs(x)))
  keep     <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(60, sum(gene_var > 0.3)))]
  if (length(keep) < 5) { message("Too few genes for co-expression network."); return() }
  mat <- mat[keep, , drop = FALSE]

  # Pearson correlation
  cor_mat <- cor(t(mat), method = "pearson")
  diag(cor_mat) <- 0

  # Edges above threshold
  edges <- which(abs(cor_mat) > CORR_CUTOFF, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
  if (nrow(edges) == 0) {
    message("No co-expression edges at r > ", CORR_CUTOFF, " – lowering to 0.3")
    edges <- which(abs(cor_mat) > 0.3, arr.ind = TRUE)
    edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
  }
  if (nrow(edges) == 0) { message("No edges found."); return() }

  edge_df <- data.frame(
    from   = rownames(cor_mat)[edges[, 1]],
    to     = rownames(cor_mat)[edges[, 2]],
    weight = cor_mat[edges]
  )

  node_df <- data.frame(
    name     = keep,
    mean_fc  = rowMeans(mat),
    stringsAsFactors = FALSE
  )

  g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)

  # Module detection by Louvain communities
  comm  <- cluster_louvain(g)
  V(g)$module <- as.character(membership(comm))

  n_mods    <- length(unique(V(g)$module))
  mod_cols  <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(n_mods),
    unique(V(g)$module)
  )

  set.seed(123)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(alpha = abs(weight), color = weight > 0),
      linewidth = 0.5
    ) +
    geom_node_point(
      aes(fill = module, size = abs(mean_fc)),
      shape = 21, color = "grey30", stroke = 0.4
    ) +
    geom_node_label(
      aes(label = name),
      size         = 2.2,
      repel        = TRUE,
      label.padding = unit(0.1, "lines"),
      label.size   = 0.1,
      max.overlaps = 25,
      show.legend  = FALSE
    ) +
    scale_edge_color_manual(
      values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
      name   = "Correlation",
      labels = c("Positive", "Negative")
    ) +
    scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    scale_fill_manual(values = mod_cols, name = "Module") +
    scale_size_continuous(name = "|mean logFC|", range = c(2, 8)) +
    labs(
      title    = paste("Gene Co-expression Network –",
                       str_replace_all(comparison, "_", " ") %>% str_to_title()),
      subtitle = paste0("Edges: |r| > ", CORR_CUTOFF,
                        "  ·  Louvain modules  ·  n = ", length(keep), " genes")
    ) +
    theme_graph(base_family = "sans", base_size = 9) +
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0.5),
      legend.position = "right"
    )

  pdf(outfile, width = 12, height = 10)
  print(p)
  invisible(dev.off())
  message("Saved: ", outfile)
}

# ── Run ────────────────────────────────────────────────────────────────────────
cnet_heat <- build_cnet_graph(
  file.path(DATADIR, "enrich_heatvsctrl.rds"),
  "Heat Stress vs Control"
)
save_cnet_pdf(cnet_heat, file.path(OUTDIR, "network_cnet_heatvsctrl.pdf"))

cnet_age <- build_cnet_graph(
  file.path(DATADIR, "enrich_agedvsadult.rds"),
  "Aged vs Adult"
)
save_cnet_pdf(cnet_age, file.path(OUTDIR, "network_cnet_agedvsadult.pdf"))

build_coexpr_pdf(
  "comp1_heatvsctrl",
  file.path(OUTDIR, "network_coexpr_heatstress.pdf")
)
