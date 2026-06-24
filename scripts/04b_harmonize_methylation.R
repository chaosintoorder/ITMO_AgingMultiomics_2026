# ============================================================================
# 04b_harmonize_methylation.R — ОБЩИЙ набор CpG для всех метиломных датасетов.
# Зачем: в 04 каждый датасет отбирал СВОИ топ-N вариабельных CpG (дисперсия на
#   сырых бетах ловит техшум, разный по датасетам) -> пересечение почти пустое.
# Решение: отбираем CpG, вариабельные ВО ВСЕХ датасетах (медианный ранг дисперсии).
#   Проход 1: читаем беты раз, полную M на диск (temp) + дисперсии (память ~10 МБ).
#   Проход 2: подменяем $M в 04_<GSE>.rds на общий набор (доли клеток/пол сохраняем).
# Память: одна матрица за раз (~3 ГБ). Запускать ПОСЛЕ 04 (5 датасетов готовы).
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })

meta    <- data.table::as.data.table(readRDS(file.path(PATHS$interim_meth, "02_meta_meth.rds")))
targets <- METH_GSE
tmpdir  <- file.path(PATHS$interim_meth, "_tmpM"); dir.create(tmpdir, showWarnings = FALSE)

read_full_M <- function(g) {                    # читает беты -> фильтрует -> полная M
  spec <- METH_BETA[[g]]; msub <- meta[gse == g]
  age  <- setNames(msub$Age, msub$gsm); gsm_ok <- msub[!is.na(Age), gsm]
  cap  <- if (!is.null(META_CAP_OVERRIDE[[g]])) META_CAP_OVERRIDE[[g]] else MAX_PER_DATASET
  keep <- cap_samples(gsm_ok, age[gsm_ok], cap)
  files <- vapply(spec$files, find_file, character(1), dir = PATHS$raw_meth)
  meta_g <- read_meth_meta(g, PATHS$raw_meth)
  if (spec$type == "series") B <- read_series_betas(files[1], keep)
  else { rd <- function(f) { cm <- map_keep_to_filecols(keep, beta_colnames(f), meta_g); read_betas_file(f, cm) }
         B <- if (spec$type == "suppl2") combine_two(rd(files[1]), rd(files[2])) else rd(files[1]) }
  storage.mode(B) <- "double"; keep <- intersect(keep, colnames(B)); B <- B[, keep, drop = FALSE]
  B <- B[!duplicated(rownames(B)), , drop = FALSE]                       # дубли (GSE87571: две матрицы)
  B <- B[!rownames(B) %in% union(sex_chr_probes_450k(), xreactive_probes_450k()), , drop = FALSE]
  B <- B[rowMeans(is.na(B)) <= NA_PROBE_MAX, , drop = FALSE]
  M <- beta_to_m(B); rm(B); gc(); M
}

# ---- проход 1: полная M на диск + дисперсии ----
varlist <- list()
for (g in targets) {
  log_msg("=== ", g, " (проход 1: читаю беты) ===")
  M <- read_full_M(g)
  saveRDS(M, file.path(tmpdir, paste0(g, ".rds")))
  varlist[[g]] <- setNames(matrixStats::rowVars(M, na.rm = TRUE), rownames(M))
  log_msg("  полная M ", nrow(M), " x ", ncol(M), " -> temp")
  rm(M); gc()
}

# ---- общий набор: топ-N по медианному рангу дисперсии (CpG, измеренные во ВСЕХ) ----
common_all <- Reduce(intersect, lapply(varlist, names))
log_msg("CpG, измеренных во ВСЕХ датасетах: ", length(common_all))
rankmat <- sapply(varlist, function(v) rank(-v[common_all]))   # ранг 1 = самый вариабельный
agg <- rowMeans(rankmat)
common <- names(sort(agg))[seq_len(min(N_VAR_CPG, length(agg)))]
log_msg("ОБЩИЙ слим-набор CpG: ", length(common))

# ---- проход 2: подменяем $M в слим-файлах ----
for (g in targets) {
  M <- readRDS(file.path(tmpdir, paste0(g, ".rds")))[common, , drop = FALSE]
  if (anyNA(M) && requireNamespace("impute", quietly = TRUE)) M <- impute::impute.knn(M, k = 10)$data
  r <- readRDS(file.path(PATHS$interim_meth, paste0("04_", g, ".rds")))
  cols <- colnames(M); r$M <- M
  r$age <- r$age[cols]; r$sex <- r$sex[match(cols, names(r$sex))]
  if (!is.null(r$cellprop)) r$cellprop <- r$cellprop[cols, , drop = FALSE]
  save_rds(r, file.path(PATHS$interim_meth, paste0("04_", g, ".rds")))
  log_msg("  ", g, ": $M -> общий набор ", nrow(M), " x ", ncol(M))
  rm(M, r); gc()
}
unlink(tmpdir, recursive = TRUE)

# ---- заодно понятный QC-график пола: столбики reported/predicted по датасетам ----
obj <- lapply(file.path(PATHS$interim_meth, paste0("04_", targets, ".rds")), readRDS)
sd  <- do.call(rbind, lapply(obj, function(x) x$sexrow))
sd$status <- ifelse(sd$predicted, "predicted", "reported")
agg2 <- as.data.frame(table(dataset = sd$gse, Sex = sd$Sex, status = sd$status))
agg2 <- agg2[agg2$Freq > 0 & !is.na(agg2$Sex), ]
p <- ggplot(agg2, aes(dataset, Freq, fill = interaction(Sex, status, sep = " / "))) +
  geom_col(position = position_stack()) + coord_flip() +
  scale_fill_manual(values = c("F / reported"="#b5495b","M / reported"="#3b6ea5",
                               "F / predicted"="#e3a3ae","M / predicted"="#9db9d6"), name = NULL) +
  labs(title = "Methylation: sex composition per dataset",
       x = NULL, y = "Samples, count") + theme_soda()
save_plot(p, file.path(PATHS$pics, "04_sex_composition_methylation"), w = 8, h = 4.5)
log_msg("04b готов: общий набор ", length(common), " CpG | QC-пол: 04_sex_composition_methylation")
