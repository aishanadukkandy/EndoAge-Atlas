
# Linear age effects in EC pseudobulk, mice with > 100 ECs


.libPaths(c("C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/Aging/renv/library/windows/R-4.4/x86_64-w64-mingw32", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
})

# ----------------------------
# PARAMETERS
# ----------------------------
input_dir <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/Annotation"
base_dir  <- "C:/Users/au730762/OneDrive - Aarhus universitet/Documents/Aging/Linear_age_v1_min100"
sp_dir    <- file.path(base_dir, "Spearman")
lm_dir    <- file.path(base_dir, "LinearRegression")
for (d in c(sp_dir, lm_dir, file.path(lm_dir, "Scaled"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

organs    <- c("BAT", "gWAT", "Heart", "iWAT", "Kidney", "Liver", "Muscle", "Stomach")
min_cells <- 100
n_perm    <- 1000
seed      <- 20261026

age_map <- c("03_months"=3, "06_months"=6, "12_months"=12, "16_months"=16, "23_months"=23)

# ----------------------------
# MAIN PER ORGAN
# ----------------------------
mice_tables <- list()

for (organ in organs){
  message("▶ Linear age v1 (> ", min_cells, " ECs): ", organ)
  set.seed(seed + match(organ, organs))

  seu_file <- list.files(file.path(input_dir, organ), pattern = "_onlyECs_v3\\.rds$", full.names = TRUE)
  if (length(seu_file) != 1) stop("No unique _onlyECs_v3.rds file for ", organ)
  seu <- readRDS(seu_file)
  seu$Age_num <- as.numeric(age_map[seu$Age_group])

  # Remove mice with <= min_cells ECs
  n_cells <- table(as.character(seu$ID))
  keep <- names(n_cells)[n_cells > min_cells]
  mice_tables[[organ]] <- data.frame(Organ = organ, ID = names(n_cells), n_cells = as.numeric(n_cells),
                                     included = names(n_cells) %in% keep)
  seu <- subset(seu, cells = colnames(seu)[as.character(seu$ID) %in% keep])
  message("  ", length(n_cells) - length(keep), " mice left out; ", length(keep), " mice analysed")

  # pseudobulk per mouse ID 
  expr <- AggregateExpression(seu, group.by = "ID", assays = "RNA", slot = "data")$RNA

  meta <- seu@meta.data %>%
    dplyr::select(ID, Age_num, Sex) %>%
    distinct() %>%
    arrange(ID)
  # Seurat sometimes prepends 'g' to IDs
  if (any(!meta$ID %in% colnames(expr))) meta$ID <- paste0("g", meta$ID)
  expr <- as.matrix(expr[, meta$ID, drop = FALSE])
  rm(seu); gc()

  n    <- nrow(meta)
  age  <- meta$Age_num
  ages_perm <- replicate(n_perm, sample(age))          # mice x permutations, shared by both tests

  # --- Spearman (as cor.test / cor with method = "spearman") ---
  rk  <- t(apply(expr, 1, rank))
  rkc <- rk - rowMeans(rk)
  rk_norm <- sqrt(rowSums(rkc^2))                      # 0 for genes with constant expression -> NA rho
  rank_cor <- function(a_mat) {
    ra <- apply(a_mat, 2, rank)
    ra <- sweep(ra, 2, colMeans(ra))
    (rkc %*% ra) / outer(rk_norm, sqrt(colSums(ra^2)))
  }
  rho <- as.numeric(rank_cor(matrix(age, ncol = 1)))
  rho[!is.finite(rho)] <- NA
  pval <- 2 * pt(-abs(rho / sqrt((1 - rho^2) / (n - 2))), df = n - 2)
  rho_perm <- rank_cor(ages_perm)
  perm_pval <- (rowSums(abs(rho_perm) >= abs(rho)) + 1) / (n_perm + 1)

  sp <- data.frame(Gene = rownames(expr), rho.rho = rho, pval = pval,
                   perm_pval = perm_pval, perm_fdr = p.adjust(perm_pval, method = "fdr")) %>%
    arrange(pval)
  write.xlsx(sp, file.path(sp_dir, paste0("Spearman_perm_", organ, "_onlyECs.xlsx")), rowNames = FALSE)

  # --- Linear regression ---
  S <- model.matrix(~ factor(meta$Sex))
  P <- S %*% solve(crossprod(S)) %*% t(S)
  Ex <- expr - expr %*% P                              # expression with sex (and intercept) removed
  age_res <- function(a_mat) a_mat - P %*% a_mat
  Ea  <- age_res(matrix(age, ncol = 1))
  ssa <- sum(Ea^2)
  beta <- as.numeric(Ex %*% Ea) / ssa
  rss  <- rowSums(Ex^2) - beta^2 * ssa
  se   <- sqrt(pmax(rss, 0) / (n - 3) / ssa)
  p_lm <- 2 * pt(-abs(beta / se), df = n - 3)

  Ea_perm   <- age_res(ages_perm)
  beta_perm <- (Ex %*% Ea_perm) / matrix(colSums(Ea_perm^2), nrow(Ex), n_perm, byrow = TRUE)
  perm_p_lm <- (rowSums(abs(beta_perm) >= abs(beta)) + 1) / (n_perm + 1)

  lmr <- data.frame(Gene = rownames(expr), beta.Estimate = beta, `se.Std. Error` = se,
                    `pval.Pr(>|t|)` = p_lm, perm_pval = perm_p_lm,
                    perm_fdr = p.adjust(perm_p_lm, "fdr"), check.names = FALSE)
  write.xlsx(lmr, file.path(lm_dir, paste0(organ, "_linear_age_adjSex.xlsx")), rowNames = FALSE)

  # beta scaled to [-1, 1] (as v1 "Scaled")
  lmr$beta_scaled <- lmr$beta.Estimate / max(abs(lmr$beta.Estimate), na.rm = TRUE)
  names(lmr)[names(lmr) == "se.Std. Error"] <- "se.Std..Error"
  write.xlsx(lmr, file.path(lm_dir, "Scaled", paste0(organ, "_linear_age_adjSex_scaled.xlsx")), rowNames = FALSE)

  message("✔ ", organ, ": ", sum(sp$perm_pval < 0.05, na.rm = TRUE), " Spearman genes and ",
          sum(lmr$perm_pval < 0.05, na.rm = TRUE), " lm genes with perm p < 0.05")
  rm(expr, rk, rkc, rho_perm, Ex, beta_perm); gc()
}

write.xlsx(bind_rows(mice_tables) %>% arrange(Organ, ID),
           file.path(base_dir, "Mice_included_v1_min100.xlsx"), rowNames = FALSE)

message("Linear age v1 (> 100 ECs) finished. Results in: ", base_dir)
