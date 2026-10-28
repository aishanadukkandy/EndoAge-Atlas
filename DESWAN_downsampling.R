
# DE-SWAN  / robustness to unequal group sizes (Reviewer 1, comment 4)
#
# mice (> 100 ECs per organ).
# Every age group is downsampled to 2 females + 2 males (all mice if fewer),
# 100 times. Each draw re-runs the DE-SWAN windows (limma-voom; v6: ~ 0 + group,
# No permutations inside the draws: the question is whether the crest pattern
# holds when every age has the same number of mice.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(readxl)
  library(limma)
  library(edgeR)
  library(parallel)
})

# ----------------------------
# PARAMETERS
# ----------------------------
v4_dir     <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs/v4"
pb_dir     <- file.path(v4_dir, "pseudobulk")
args       <- commandArgs(trailingOnly = TRUE)
version    <- if (length(args) > 0) args[1] else "v6"
use_sex    <- grepl("sex", version)
ver_dir    <- file.path("C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/DESWAN/onlyECs", version)
output_dir <- file.path(ver_dir, "downsampling")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Mice with this many ECs or fewer in an organ are left out (as v6)
min_cells <- 100

organs <- c("BAT", "gWAT", "Heart", "iWAT", "Kidney", "Liver", "Muscle", "Stomach")

age_levels <- c("03_months", "06_months", "12_months", "16_months", "23_months")
age_map    <- c("03_months"=3, "06_months"=6, "12_months"=12, "16_months"=16, "23_months"=23)

k <- 1
fdr_thresh <- 0.05
lfc_thresh <- 0.25
min_ids_per_side <- 2

n_per_sex <- 2      # mice per sex per age group in each draw
n_draws   <- 100
n_cores   <- 12

base_family <- "Calibri"

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

# ----------------------------
# HELPERS 
# ----------------------------
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

run_window_DE <- function(counts, meta, left_ages, right_ages){
  keep_ids <- meta$ID[meta$Age_group %in% c(left_ages, right_ages)]
  y <- counts[, keep_ids, drop=FALSE]
  m <- meta[keep_ids, , drop=FALSE]

  group <- factor(ifelse(m$Age_group %in% left_ages, "Left", "Right"), levels = c("Left", "Right"))
  sex   <- droplevels(factor(m$Sex))

  # sex covariate for v6sex (dropped only if a window holds a single sex)
  if (use_sex && nlevels(sex) > 1) {
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

  sum((tt$adj.P.Val < fdr_thresh) & (abs(tt$logFC) >= lfc_thresh))
}

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

# One draw: 2 females + 2 males per age group, then all windows
one_draw <- function(d, counts, meta, windows){
  ids <- meta %>%
    group_by(Age_group, Sex) %>%
    group_modify(~ .x[sample(nrow(.x), min(n_per_sex, nrow(.x))), , drop = FALSE]) %>%
    ungroup() %>%
    pull(ID)
  m <- meta[ids, , drop = FALSE]

  out <- lapply(seq_len(nrow(windows)), function(i) {
    left_ids  <- m$ID[m$Age_group %in% windows$left_ages[[i]]]
    right_ids <- m$ID[m$Age_group %in% windows$right_ages[[i]]]
    n_sig <- if (length(left_ids) < min_ids_per_side || length(right_ids) < min_ids_per_side) NA else
      run_window_DE(counts, m, windows$left_ages[[i]], windows$right_ages[[i]])
    data.frame(draw = d, center_age = windows$center_age[i], n_sig_total = n_sig,
               n_left = length(left_ids), n_right = length(right_ids),
               mice = paste(c(left_ids, right_ids), collapse = ";"))
  }) %>% bind_rows()

  # A crest needs at least one significant gene (ties of zeros are not crests)
  out <- find_crests(out)
  out$crest <- out$crest & out$n_sig_total > 0
  out$top   <- out$n_sig_total == max(out$n_sig_total) & out$n_sig_total > 0
  out
}

# ----------------------------
# MAIN
# ----------------------------
windows <- make_windows(age_levels, k = k)

cl <- makeCluster(n_cores)
clusterEvalQ(cl, suppressPackageStartupMessages({ library(dplyr); library(limma); library(edgeR) }))
clusterExport(cl, c("run_window_DE", "find_crests", "fdr_thresh", "lfc_thresh",
                    "min_ids_per_side", "n_per_sex", "use_sex"))
clusterSetRNGStream(cl, 20261024)

draws_all <- list()
for (organ in organs){
  message("▶ Downsampling (", version, "): ", organ)
  pb <- readRDS(file.path(pb_dir, paste0(organ, "_pseudobulk_7000HVG.rds")))
  meta <- pb$meta[pb$meta$n_cells > min_cells, , drop = FALSE]
  meta$Age_group <- factor(meta$Age_group, levels = age_levels)

  res <- parLapply(cl, seq_len(n_draws), one_draw,
                   counts = pb$counts, meta = meta, windows = windows)
  draws_all[[organ]] <- bind_rows(res) %>% mutate(Organ = organ)
}
stopCluster(cl)

draws_all <- bind_rows(draws_all) %>%
  mutate(center_age = factor(center_age, levels = age_levels))

# Full-data v3 results for reference
full_summary <- lapply(organs, function(organ) {
  readxl::read_xlsx(file.path(ver_dir, paste0(organ, "_DESWAN_summary_withperm_", version, ".xlsx"))) %>%
    transmute(Organ = organ, center_age = factor(center_age, levels = age_levels),
              n_sig_full = n_sig_total, crest_full = as.logical(crest), emp_pval_full = emp_pval)
}) %>% bind_rows() %>%
  mutate(sig_crest_full = crest_full & emp_pval_full < 0.05)

by_window <- draws_all %>%
  group_by(Organ, center_age) %>%
  summarise(
    mice_per_window = paste0(median(n_left), " vs ", median(n_right)),
    median_n_sig = median(n_sig_total, na.rm = TRUE),
    q25_n_sig = quantile(n_sig_total, 0.25, na.rm = TRUE),
    q75_n_sig = quantile(n_sig_total, 0.75, na.rm = TRUE),
    min_n_sig = min(n_sig_total, na.rm = TRUE),
    max_n_sig = max(n_sig_total, na.rm = TRUE),
    pct_draws_crest = 100 * mean(crest, na.rm = TRUE),
    pct_draws_top   = 100 * mean(top, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(full_summary, by = c("Organ", "center_age")) %>%
  arrange(Organ, center_age)

crest_recovery <- by_window %>%
  filter(sig_crest_full) %>%
  dplyr::select(Organ, center_age, n_sig_full, emp_pval_full, median_n_sig,
                q25_n_sig, q75_n_sig, pct_draws_crest, pct_draws_top)

# Which mice drive large counts: for every mouse in a window, median # significant genes
# in draws that include it vs draws that leave it out
mouse_meta <- lapply(organs, function(organ) {
  readRDS(file.path(pb_dir, paste0(organ, "_pseudobulk_7000HVG.rds")))$meta %>%
    transmute(Organ = organ, ID = as.character(ID), Age_group = as.character(Age_group), Sex)
}) %>% bind_rows()

mouse_influence <- draws_all %>%
  dplyr::select(Organ, center_age, draw, n_sig_total, mice) %>%
  mutate(ID = strsplit(mice, ";")) %>%
  dplyr::select(-mice) %>%
  unnest(ID) %>%
  mutate(included = TRUE) %>%
  complete(nesting(Organ, center_age, draw, n_sig_total), ID, fill = list(included = FALSE)) %>%
  semi_join(mouse_meta, by = c("Organ", "ID")) %>%
  group_by(Organ, center_age, ID) %>%
  filter(any(included)) %>%
  summarise(
    n_draws_in = sum(included),
    median_n_sig_in  = median(n_sig_total[included], na.rm = TRUE),
    median_n_sig_out = median(n_sig_total[!included], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(mouse_meta, by = c("Organ", "ID")) %>%
  mutate(ratio_in_out = (median_n_sig_in + 1) / (median_n_sig_out + 1)) %>%
  arrange(desc(ratio_in_out))

write.xlsx(
  list(crest_recovery = crest_recovery, by_window = by_window,
       mouse_influence = mouse_influence, per_draw = draws_all),
  file.path(output_dir, paste0("DESWAN_", version, "_downsampling_2F2M_100draws.xlsx")), rowNames = FALSE
)

# ----------------------------
# PLOTS
# ----------------------------
# Significant genes per window across draws; diamond = full data; label = % draws as crest
lab_df <- by_window %>%
  group_by(Organ) %>%
  mutate(y_lab = max(c(max_n_sig, n_sig_full), na.rm = TRUE) * 1.12) %>%
  ungroup()

p_ds <- ggplot(draws_all, aes(x = center_age, y = n_sig_total)) +
  geom_boxplot(aes(fill = Organ), outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.35) +
  geom_point(data = by_window, aes(y = n_sig_full), shape = 23, size = 2.6,
             fill = "white", colour = "black", stroke = 0.9) +
  geom_text(data = lab_df,
            aes(y = y_lab, label = paste0(round(pct_draws_crest), "%"),
                fontface = ifelse(sig_crest_full, "bold", "plain")),
            size = 3.2, family = base_family) +
  scale_fill_manual(values = organ_colors, guide = "none") +
  scale_x_discrete(labels = function(x) paste0(as.numeric(age_map[x]), " m")) +
  facet_wrap(~ Organ, ncol = 4, scales = "free_y") +
  labs(
    title = "DE-SWAN crests after downsampling to 2 females + 2 males per age (100 draws)",
    subtitle = paste0("Boxes: significant genes per draw. Diamond: full data (", version, "). Label: % of draws in which the window is a crest (bold = significant crest in full data)"),
    x = "Center age",
    y = paste0("# significant genes (FDR<", fdr_thresh, ", |log2FC|≥", lfc_thresh, ")")
  ) +
  theme_minimal(base_size = 12, base_family = base_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, paste0("DESWAN_", version, "_downsampling_by_window.png")), p_ds, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, paste0("DESWAN_", version, "_downsampling_by_window.pdf")), p_ds, width = 12, height = 8, device = cairo_pdf)

# Recovery of the significant crests
p_rec <- crest_recovery %>%
  mutate(Crest = paste0(Organ, " ", as.numeric(age_map[as.character(center_age)]), " m")) %>%
  ggplot(aes(x = reorder(Crest, pct_draws_crest), y = pct_draws_crest, fill = Organ)) +
  geom_col(colour = "black", width = 0.7) +
  geom_hline(yintercept = 80, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = organ_colors, guide = "none") +
  coord_flip(ylim = c(0, 100)) +
  labs(title = "Recovery of significant DE-SWAN crests after downsampling",
       x = NULL, y = "% of 100 draws in which the window is a crest") +
  theme_minimal(base_size = 12, base_family = base_family) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(output_dir, paste0("DESWAN_", version, "_crest_recovery.png")), p_rec, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, paste0("DESWAN_", version, "_crest_recovery.pdf")), p_rec, width = 7, height = 5, device = cairo_pdf)

message("Downsampling finished. Results in: ", output_dir)
