
# DE-SWAN style sliding-window pipeline (pseudobulk) 
#
# DESWAN/DESWAN_analysis.R: top 7000 HVGs, same windows, limma-voom,
# 1000 permutations, same crest rule and plots), with these changes for the revision:
#   only mice with > 100 ECs in the organ are pseudobulked and analysed
#   sex is a covariate (~ 0 + group + sex); permutations keep each mouse's sex
#   fixed seed per organ, so the permutation p-values do not depend on organ order
#   GO per crest uses only the genes of that crest

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
output_dir <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v6sex"
base_dir   <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v6"
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

# DE thresholds 
fdr_thresh <- 0.05
lfc_thresh <- 0.25

# Plot font preference
base_family <- "Calibri"

n_hvg <- 7000

# Minimum samples per side of a window
min_ids_per_side <- 2

n_perm  <- 1000
n_cores <- 12
seed    <- 20261027

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

# Windows around each available center, including edges (will auto-shrink) - as in v3
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

# limma-voom DE for a window: Left vs Right, adjusted for sex (as v4)
run_window_DE <- function(counts, meta, left_ages, right_ages){
  keep_ids <- meta$ID[meta$Age_group %in% c(left_ages, right_ages)]
  y <- counts[, keep_ids, drop=FALSE]
  m <- meta[keep_ids, , drop=FALSE]

  group <- factor(ifelse(m$Age_group %in% left_ages, "Left", "Right"), levels = c("Left", "Right"))
  sex   <- droplevels(factor(m$Sex))

  # sex covariate (dropped only if a window holds a single sex)
  if (nlevels(sex) > 1) {
    design <- model.matrix(~ 0 + group + sex)
  } else {
    design <- model.matrix(~ 0 + group)
  }
  colnames(design)[1:2] <- c("Left","Right")
  colnames(design) <- make.names(colnames(design))

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

# Find crests: local maxima in #significant genes (same rule as v3)
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

# One permutation: shuffle age labels between mice (sex stays with the mouse)
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
  message("▶ DESWAN v6sex: ", organ)
  clusterSetRNGStream(cl, seed + match(organ, organs))

  pb <- make_pseudobulk(organ, n_hvg = n_hvg)
  meta <- pb$meta
  meta$Age_group <- factor(meta$Age_group, levels = age_levels)
  mice_tables[[organ]] <- meta %>% mutate(Organ = organ, included = n_cells > min_cells)

  # v6sex: drop mice with <= min_cells ECs
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

    out_xlsx <- file.path(output_dir, paste0(organ, "_DESWAN_results_center_", center_age, "_v6sex.xlsx"))
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
  write.xlsx(summary_df, file.path(output_dir, paste0(organ, "_DESWAN_summary_v6sex.xlsx")), rowNames = FALSE)

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

  write.xlsx(summary_df, file.path(output_dir, paste0(organ, "_DESWAN_summary_withperm_v6sex.xlsx")), rowNames = FALSE)

  