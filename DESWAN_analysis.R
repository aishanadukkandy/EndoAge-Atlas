
# DE-SWAN style sliding-window pipeline (pseudobulk) 
#
#   v6: only mice with > 100 ECs in the organ are pseudobulked and analysed
#   v6: model as in v3 (~ 0 + group, no sex covariate); HVGs 7000
#   v6: fixed seed per organ, so the permutation p-values do not depend on organ order
#   v6: GO per crest uses only the genes of that crest

.libPaths(c("C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/Aging/renv/library/windows/R-4.4/x86_64-w64-mingw32", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(readxl)
  library(limma)
  library(edgeR)
  library(patchwork)
  library(parallel)
})

# ----------------------------
# PARAMETERS
# ----------------------------
input_dir  <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/Annotation"
output_dir <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v6"
v3_dir     <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v3"
pb_dir     <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v4/pseudobulk"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Mice with this many ECs or fewer in an organ are left out
min_cells <- 100

organs <- c("BAT", "gWAT", "Heart", "iWAT", "Kidney", "Liver", "Muscle", "Stomach")

# Age levels in order (factor order is important)
age_levels <- c("03_months", "06_months", "12_months", "16_months", "23_months")
age_map    <- c("03_months"=3, "06_months"=6, "12_months"=12, "16_months"=16, "23_months"=23)

# Window half-width (k ages on each side of center)
k <- 1

# DE thresholds (v3 script had lfc 0.1; every FDR<0.05 gene in v3 had |log2FC| > 0.4,
# so 0.25 as stated in the Methods gives the same calls)
fdr_thresh <- 0.05
lfc_thresh <- 0.25

# Plot font preference
base_family <- "Calibri"

n_hvg <- 7000

# Minimum samples per side of a window
min_ids_per_side <- 2

n_perm  <- 1000
n_cores <- 12
seed    <- 20261026

# ----------------------------
# HELPERS
# ----------------------------

pick_seurat_file <- function(folder){
  only <- list.files(folder, pattern = "onlyECs_v3\\.rds$", full.names = TRUE)
  if (length(only) == 1) return(only)
  std  <- list.files(folder, pattern = "_ECs_v3\\.rds$", full.names = TRUE)
  if (length(std) == 1)  return(std)
  stop("No unique ECs_v3 file in: ", folder)
}

# Pseudobulk (sum counts) per mouse ID on HVGs; cached so the big objects load only once
make_pseudobulk <- function(organ, n_hvg = 7000){
  cache <- file.path(pb_dir, paste0(organ, "_pseudobulk_", n_hvg, "HVG.rds"))
  if (file.exists(cache)) return(readRDS(cache))

  seu <- readRDS(pick_seurat_file(file.path(input_dir, organ)))
  DefaultAssay(seu) <- "RNA"
  seu <- FindVariableFeatures(seu, nfeatures = n_hvg, verbose = FALSE)
  vgs <- VariableFeatures(seu)

  if (length(Layers(seu[["RNA"]], search = "counts")) > 1) seu <- JoinLayers(seu)
  cnt <- GetAssayData(seu, assay = "RNA", layer = "counts")[vgs, ]

  ids <- factor(as.character(seu$ID))
  mouse_mat <- sparse.model.matrix(~ 0 + ids)
  colnames(mouse_mat) <- levels(ids)
  counts <- as.matrix(cnt %*% mouse_mat)

  meta <- seu@meta.data %>%
    dplyr::select(ID, Age_group, Sex) %>%
    distinct() %>%
    mutate(
      ID = as.character(ID),
      Age_group = factor(Age_group, levels = age_levels),
      Sex = droplevels(factor(Sex))
    )
  if (anyDuplicated(meta$ID)) stop(organ, ": a mouse ID has more than one age/sex")
  meta <- meta[match(colnames(counts), meta$ID), , drop = FALSE]
  rownames(meta) <- meta$ID

  n_cells <- as.numeric(Matrix::colSums(mouse_mat))
  meta$n_cells <- n_cells

  pb <- list(counts = counts, meta = meta)
  saveRDS(pb, cache)
  rm(seu, cnt); gc()
  pb
}

# Windows around each available center, including edges (will auto-shrink) 
make_windows <- function(ages, k){
  centers <- seq_along(ages)
  tibble::tibble(
    center_age = ages[centers],
    left_ages  = lapply(centers, function(i) {
      if (i - k < 1) ages[1:(i-1)] else ages[(i-k):(i-1)]
    }),
    right_ages = lapply(centers, function(i) {
      if (i + k > length(ages)) ages[(i+1):length(ages)] else ages[(i+1):(i+k)]
    })
  )
}

# limma-voom DE for a window: Left vs Right (v3 model, no covariates)
run_window_DE <- function(counts, meta, left_ages, right_ages){
  keep_ids <- meta$ID[meta$Age_group %in% c(left_ages, right_ages)]
  y <- counts[, keep_ids, drop=FALSE]
  m <- meta[keep_ids, , drop=FALSE]

  group <- factor(ifelse(m$Age_group %in% left_ages, "Left", "Right"), levels = c("Left", "Right"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Left","Right")

  dge <- DGEList(y)
  dge <- calcNormFactors(dge, method = "TMM")

  v   <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  cont <- makeContrasts(RightMinusLeft = Right - Left, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  tt   <- topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)

  tt %>%
    transmute(
      gene,
      log2FC  = logFC,
      AveExpr = AveExpr,
      t       = t,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val
    )
}

# Call significance and direction
call_sig <- function(df, fdr=0.05, lfc=0.25){
  df %>%
    mutate(
      sig = (adj.P.Val < fdr) & (abs(log2FC) >= lfc),
      direction = case_when(
        sig & log2FC > 0  ~ "Up_in_Right",
        sig & log2FC < 0  ~ "Up_in_Left",
        TRUE              ~ "Non_sig"
      )
    )
}

# Find crests: local maxima in #significant genes 
find_crests <- function(summary_tbl){
  n <- nrow(summary_tbl)
  crest <- rep(FALSE, n)
  for (i in seq_len(n)){
    left_ok  <- if (i==1) TRUE else summary_tbl$n_sig_total[i] >= summary_tbl$n_sig_total[i-1]
    right_ok <- if (i==n) TRUE else summary_tbl$n_sig_total[i] >= summary_tbl$n_sig_total[i+1]
    crest[i] <- left_ok & right_ok
  }
  summary_tbl$crest <- crest
  summary_tbl
}

# One permutation: shuffle age labels between mice
perm_once <- function(b, counts, meta, windows, centers_vec){
  meta_perm <- meta
  meta_perm$Age_group <- sample(meta_perm$Age_group)
  out <- setNames(rep(0, length(centers_vec)), centers_vec)
  for (i in seq_len(nrow(windows))) {
    center_age <- windows$center_age[i]
    if (!(center_age %in% centers_vec)) next
    left_ages  <- windows$left_ages[[i]]
    right_ages <- windows$right_ages[[i]]
    left_ids  <- meta_perm$ID[meta_perm$Age_group %in% left_ages]
    right_ids <- meta_perm$ID[meta_perm$Age_group %in% right_ages]
    if (length(left_ids) < min_ids_per_side || length(right_ids) < min_ids_per_side) next
    de_perm <- call_sig(run_window_DE(counts, meta_perm, left_ages, right_ages),
                        fdr = fdr_thresh, lfc = lfc_thresh)
    out[center_age] <- sum(de_perm$sig)
  }
  out
}

# ----------------------------
# MAIN PER ORGAN
# ----------------------------
cl <- makeCluster(n_cores)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(dplyr); library(limma); library(edgeR) }))
clusterExport(cl, c("run_window_DE", "call_sig", "fdr_thresh", "lfc_thresh", "min_ids_per_side"))
mice_tables <- list()

for (organ in organs){
  message("▶ DESWAN v6: ", organ)
  clusterSetRNGStream(cl, seed + match(organ, organs))

  pb <- make_pseudobulk(organ, n_hvg = n_hvg)
  meta <- pb$meta
  meta$Age_group <- factor(meta$Age_group, levels = age_levels)
  mice_tables[[organ]] <- meta %>% mutate(Organ = organ, included = n_cells > min_cells)

  # v6: drop mice with <= min_cells ECs
  meta   <- meta[meta$n_cells > min_cells, , drop = FALSE]
  counts <- pb$counts[, meta$ID, drop = FALSE]
  message("  ", sum(!mice_tables[[organ]]$included), " mice with <= ", min_cells, " ECs left out; ",
          nrow(meta), " mice analysed")

  windows <- make_windows(age_levels, k = k)
  window_results <- list()

  for (i in seq_len(nrow(windows))){
    center_age <- windows$center_age[i]
    left_ages  <- windows$left_ages[[i]]
    right_ages <- windows$right_ages[[i]]

    left_ids  <- meta$ID[meta$Age_group %in% left_ages]
    right_ids <- meta$ID[meta$Age_group %in% right_ages]
    if (length(left_ids) < min_ids_per_side || length(right_ids) < min_ids_per_side) {
      message("  Skipping center ", center_age, " (insufficient samples: L=", length(left_ids),
              ", R=", length(right_ids), ")")
      next
    }

    de <- run_window_DE(counts, meta, left_ages, right_ages)
    de <- call_sig(de, fdr = fdr_thresh, lfc = lfc_thresh)
    de$center_age <- center_age
    de$left_ages  <- paste(left_ages, collapse = "+")
    de$right_ages <- paste(right_ages, collapse = "+")
    window_results[[center_age]] <- de

    out_xlsx <- file.path(output_dir, paste0(organ, "_DESWAN_results_center_", center_age, "_v6.xlsx"))
    write.xlsx(de, out_xlsx, rowNames = FALSE)
  }

  if (length(window_results) == 0){
    warning("No window results for ", organ)
    next
  }

  summary_df <- bind_rows(window_results) %>%
    group_by(center_age) %>%
    summarise(
      n_sig_total = sum(sig),
      n_up_right  = sum(direction == "Up_in_Right"),
      n_up_left   = sum(direction == "Up_in_Left"),
      .groups = "drop"
    ) %>%
    mutate(
      center_age   = factor(center_age, levels = age_levels),
      center_month = as.numeric(age_map[as.character(center_age)])
    ) %>%
    arrange(center_age)

  summary_df <- find_crests(summary_df)
  write.xlsx(summary_df, file.path(output_dir, paste0(organ, "_DESWAN_summary_v6.xlsx")), rowNames = FALSE)

  # ---------------------------------
  # Permutation test (counts of sig genes per center)
  # ---------------------------------
  centers_vec <- as.character(summary_df$center_age)
  perm_list <- parLapply(cl, seq_len(n_perm), perm_once,
                         counts = counts, meta = meta, windows = windows, centers_vec = centers_vec)
  perm_counts <- do.call(rbind, perm_list)

  empirical_p <- sapply(centers_vec, function(age) {
    obs <- summary_df$n_sig_total[summary_df$center_age == age]
    null <- perm_counts[, age]
    (sum(null >= obs) + 1) / (n_perm + 1)
  })
  summary_df$emp_pval <- as.numeric(empirical_p)

  write.xlsx(summary_df, file.path(output_dir, paste0(organ, "_DESWAN_summary_withperm_v6.xlsx")), rowNames = FALSE)

  # ----------------------------
  # PLOTS
  # ----------------------------
  top_mat_df <- bind_rows(window_results) %>%
    mutate(center_age = factor(center_age, levels = age_levels)) %>%
    arrange(center_age)

  sig_genes <- top_mat_df %>% filter(sig) %>% pull(gene) %>% unique()

  if (length(sig_genes) > 1){
    mat <- top_mat_df %>%
      filter(gene %in% sig_genes) %>%
      dplyr::select(gene, center_age, log2FC) %>%
      tidyr::pivot_wider(names_from = center_age, values_from = log2FC) %>%
      as.data.frame()

    rownames(mat) <- mat$gene
    mat$gene <- NULL

    mat_scaled <- t(scale(t(as.matrix(mat))))

    mat_long <- as.data.frame(mat_scaled) %>%
      tibble::rownames_to_column("gene") %>%
      tidyr::pivot_longer(-gene, names_to = "center", values_to = "z") %>%
      mutate(center = factor(center, levels = age_levels))

    p_hm <- ggplot(mat_long, aes(x = center, y = gene, fill = z)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(title = paste0("DE-SWAN log2FC (Z-scaled) across centers: ", organ),
           x = "Center age", y = "Genes (sig in any center)") +
      theme_minimal(base_size = 12, base_family = base_family) +
      theme(
        plot.title = element_text(face="bold", hjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )

    ggsave(file.path(output_dir, paste0(organ, "_DESWAN_heatmap_sig_genes_withperm_v6.png")),
           p_hm, width = 7.5, height = 8, dpi = 300)
  }

  message("✔ Finished DESWAN v6 for ", organ)
}

stopCluster(cl)

# Mice per organ x age x sex after the EC filter, and the mice left out
mice_all <- bind_rows(mice_tables)
write.xlsx(
  list(
    per_mouse = mice_all %>% dplyr::select(Organ, ID, Age_group, Sex, n_cells, included) %>% arrange(Organ, Age_group, Sex, ID),
    excluded = mice_all %>% filter(!included) %>% dplyr::select(Organ, ID, Age_group, Sex, n_cells) %>% arrange(Organ, Age_group, Sex, ID),
    mice_per_group = mice_all %>% filter(included) %>% count(Organ, Age_group, Sex) %>%
      tidyr::pivot_wider(names_from = Sex, values_from = n, values_fill = 0)
  ),
  file.path(output_dir, "Mice_per_organ_age_sex_v6.xlsx"), rowNames = FALSE
)

message("All done per-organ. Results in: ", output_dir)

# ============================================
# COMBINE SUMMARIES & WAVE PLOT
# ============================================

read_summaries <- function(dir, suffix){
  files <- list.files(dir, pattern = paste0("_DESWAN_summary_withperm_", suffix, "\\.xlsx$"), full.names = TRUE)
  lapply(files, function(f) {
    df <- readxl::read_xlsx(f, col_names = TRUE)
    df$Organ <- sub(paste0("_DESWAN_summary_withperm_", suffix, "\\.xlsx$"), "", basename(f))
    df
  }) %>% bind_rows()
}

summary_all <- read_summaries(output_dir, "v6") %>%
  mutate(
    center_age   = factor(center_age, levels = age_levels),
    center_month = as.numeric(center_month)
  ) %>%
  arrange(Organ, center_age)

if (nrow(summary_all) == 0) {
  stop("Combined summary is empty—likely no windows passed the min sample threshold.")
}

p_wave <- ggplot(summary_all,
                 aes(x = center_month, y = n_sig_total, color = Organ, group = Organ)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = crest), size = 2) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), na.translate = FALSE) +
  labs(title = "DE-SWAN waves across organs",
       x = "Age (months)",
       y = "# significant genes",
       color = "Organ", shape = "Crest") +
  theme_minimal(base_size = 14, base_family = base_family) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "DE_SWAN_waves_all_organs_withperm_v6.png"),
       p_wave, width = 8, height = 6, dpi = 300)

organ_colors <- c(
  "BAT"        = "#e31a1c",
  "iWAT"       = "#ffff99",
  "gWAT"       = "#ff7f00",
  "Heart"      = "#fdb462",
  "Kidney"     = "#cab2d6",
  "Liver"      = "#8dd3c7",
  "Muscle"     = "#b2df8a",
  "Stomach"    = "#1F78B4"
)

summary_z <- summary_all %>%
  group_by(Organ) %>%
  mutate(
    n_sig_z = {
      m <- mean(n_sig_total, na.rm = TRUE)
      s <- sd(n_sig_total, na.rm = TRUE)
      if (is.na(s) || s == 0) 0 else (n_sig_total - m) / s
    }
  ) %>%
  ungroup()

p_wave_z <- ggplot(summary_z,
                   aes(x = center_month, y = n_sig_z,
                       color = Organ, group = Organ)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = crest), size = 2) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                     na.translate = FALSE) +
  scale_color_manual(values = organ_colors) +
  labs(title = "DE-SWAN waves across organs (z-scored within organ)",
       x = "Age (months)",
       y = "z-score of # significant genes",
       color = "Organ",
       shape = "Crest") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave(file.path(output_dir, "DE_SWAN_waves_all_organs_zscore_colored.png"),
       p_wave_z, width = 8, height = 6, dpi = 300)

# ========= BARPLOTS (per organ, using organ_colors) =========
library(purrr)

plot_organ_bar <- function(org) {
  df <- summary_all %>% filter(Organ == org)
  if (nrow(df) == 0) {
    warning("No rows for organ: ", org)
    return(NULL)
  }
  p_counts <- ggplot(df, aes(x = center_age, y = n_sig_total)) +
    geom_col(width = 0.7, fill = organ_colors[[org]], color = "black") +
    geom_point(data = df %>% filter(crest == TRUE),
               aes(y = n_sig_total),
               size = 3, color = "black") +
    geom_text(data = df %>% filter(crest == TRUE),
              aes(label = "crest", y = n_sig_total),
              vjust = -0.6, size = 4, family = base_family) +
    labs(
      title = paste0("DE-SWAN crest scan: ", org),
      x = "Center age",
      y = paste0("# significant genes (FDR<", fdr_thresh, ", |log2FC|≥", lfc_thresh, ")")
    ) +
    theme_minimal(base_size = 14, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  out_file <- file.path(output_dir, paste0(org, "_DESWAN_crests_barplot_withperm_v6.png"))
  ggsave(out_file, p_counts, width = 7, height = 4.5, dpi = 300)
  invisible(p_counts)
}

orgs_in_data <- intersect(names(organ_colors), unique(summary_all$Organ))
walk(orgs_in_data, plot_organ_bar)

p_facets <- ggplot(summary_all,
                   aes(x = center_age, y = n_sig_total, fill = Organ)) +
  geom_col(width = 0.7, color = "black") +
  geom_point(data = summary_all %>% filter(crest == TRUE),
             aes(y = n_sig_total), size = 2, color = "black") +
  scale_fill_manual(values = organ_colors, guide = "none") +
  labs(
    title = "DE-SWAN crests across organs",
    x = "Center age",
    y = paste0("# significant genes (FDR<", fdr_thresh, ", |log2FC|≥", lfc_thresh, ")")
  ) +
  facet_wrap(~ Organ, ncol = 4, scales = "free_y") +
  theme_minimal(base_size = 13, base_family = base_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "DE_SWAN_crests_barplots_facet_all_organs.png"),
       p_facets, width = 11, height = 8, dpi = 300)

######################
# Find the function of these genes at the crest in each organ

library(stringr)

crest_index <- summary_all %>%
  filter(crest == TRUE) %>%
  transmute(Organ, center_age = as.character(center_age))

crest_gene_rows <- list()

for (i in seq_len(nrow(crest_index))) {
  organ <- crest_index$Organ[i]
  center_age <- crest_index$center_age[i]
  de_file <- file.path(output_dir, paste0(organ, "_DESWAN_results_center_", center_age, "_v6.xlsx"))
  if (!file.exists(de_file)) next
  de <- readxl::read_xlsx(de_file) %>% call_sig(fdr = fdr_thresh, lfc = lfc_thresh)

  crest_sig <- de %>% filter(sig) %>%
    transmute(Organ = organ,
              center_age = center_age,
              gene,
              direction)

  if (nrow(crest_sig) > 0) {
    crest_gene_rows[[paste0(organ,"_",center_age)]] <- crest_sig
  }
}

summary_all_genes <- bind_rows(crest_gene_rows)
write.xlsx(summary_all_genes, file.path(output_dir, "Crest_Significant_Genes_ByOrgan_ByCenter.xlsx"))

go_ok <- requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE)
if (!go_ok) {
  message("clusterProfiler/org.Mm.eg.db not installed - skipping GO")
} else {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })

  go_dir <- file.path(output_dir, "GO_Crest_ByCenter")
  dir.create(go_dir, showWarnings = FALSE)

  # Background = all genes tested in this window
  run_go_for_center <- function(organ, crest_age) {
    de_file <- file.path(output_dir, paste0(organ, "_DESWAN_results_center_", crest_age, "_v6.xlsx"))
    if (!file.exists(de_file)) return(NULL)
    de <- readxl::read_xlsx(de_file)

    bg_map <- bitr(unique(de$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    if (nrow(bg_map) == 0) return(NULL)

    # v6: only this crest's genes (v3 compared center_age with itself, pooling all crests of the organ)
    fg_symbols <- summary_all_genes %>%
      filter(Organ == organ, center_age == crest_age) %>%
      pull(gene) %>% unique()
    if (length(fg_symbols) < 5) return(NULL)

    fg_map <- bitr(fg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    if (nrow(fg_map) == 0) return(NULL)

    ego <- enrichGO(
      gene          = unique(fg_map$ENTREZID),
      universe      = unique(bg_map$ENTREZID),
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "fdr",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.10,
      readable      = TRUE
    )
    if (is.null(ego) || nrow(ego@result) == 0) return(NULL)

    res <- ego@result %>% mutate(Organ = organ, Crest_Age = crest_age)
    openxlsx::write.xlsx(res, file.path(go_dir, paste0("GO_Crest_", organ, "_", crest_age, ".xlsx")))
    res
  }

  pairs <- summary_all_genes %>% distinct(Organ, center_age) %>% arrange(Organ, center_age)
  go_crest_list <- list()
  for (i in seq_len(nrow(pairs))) {
    out <- run_go_for_center(pairs$Organ[i], pairs$center_age[i])
    if (!is.null(out)) go_crest_list[[paste0(pairs$Organ[i], "_", pairs$center_age[i])]] <- out
  }
  go_crest_all <- bind_rows(go_crest_list)
  openxlsx::write.xlsx(go_crest_all, file.path(go_dir, "GO_Crest_AllOrgans_ByCenter.xlsx"))

  if (nrow(go_crest_all) > 0) {
    top10 <- go_crest_all %>%
      group_by(Organ, Crest_Age) %>%
      slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        Tag = paste0(Organ, "_", Crest_Age),
        neglog10FDR = -log10(p.adjust),
        Description_wrapped = stringr::str_wrap(Description, 40)
      )

    p <- ggplot(top10, aes(x = Tag, y = Description_wrapped)) +
      geom_point(aes(size = Count, color = neglog10FDR)) +
      scale_color_gradient(low = "darkblue", high = "darkred", name = "-log10(FDR)") +
      scale_size_area(max_size = 7, name = "Gene count") +
      labs(x = "Organ_Crest", y = NULL,
           title = "Top GO BP terms per Organ × Crest (crest genes only)") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(file.path(go_dir, "GO_Crest_Top10_Dotplot.png"), p, width = 10, height = 15, dpi = 300)
  }
}

message("DE-SWAN v6 finished.")
