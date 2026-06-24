# ============================================================================
# scripts/04_preprocess_methylation.R
# Вход : беты (series/suppl по манифесту METH_BETA) + 02_meta_meth.rds
# Выход: per-dataset slim M-матрицы (топ-N вариабельных CpG) -> 04_<GSE>.rds
#        графики: ridgeline и violin β-значений, предсказание пола.
# Память (8 ГБ): один датасет за раз, сразу режем, чистим RAM.
# GSE55763 (9.7 ГБ) читается ДОЛГО — при нехватке RAM ставь targets <- c("GSE55763")
#   и перезапускай R между датасетами.
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })

meta    <- data.table::as.data.table(readRDS(file.path(PATHS$interim_meth, "02_meta_meth.rds")))
#targets <- METH_GSE
#targets <- setdiff(METH_GSE, "GSE55763")
targets <- c("GSE55763")

process_meth <- function(g) {
  spec <- METH_BETA[[g]]; if (is.null(spec)) { log_msg("нет манифеста ", g); return(NULL) }
  msub <- meta[gse == g]
  age  <- setNames(msub$Age, msub$gsm); sex <- setNames(msub$Sex, msub$gsm)
  gsm_ok <- msub[!is.na(Age), gsm]; if (!length(gsm_ok)) { log_msg("  нет возраста"); return(NULL) }

  cap <- if (!is.null(META_CAP_OVERRIDE[[g]])) META_CAP_OVERRIDE[[g]] else MAX_PER_DATASET
  keep <- cap_samples(gsm_ok, age[gsm_ok], cap)
  log_msg("  образцов с возрастом=", length(gsm_ok), " -> ", length(keep), " (cap=", cap, ")")

  files <- vapply(spec$files, find_file, character(1), dir = PATHS$raw_meth)
  meta_g <- read_meth_meta(g, PATHS$raw_meth)           # для сопоставления колонок (title->GSM)

  # --- ОГРОМНЫЙ файл (GSE55763): двухпроходный slim-поток, минуя полную матрицу ---
  if (isTRUE(spec$stream)) {
    sb <- read_betas_slim_stream(files[1], keep, meta_g, N_VAR_CPG)
    Bs <- sb$slim_beta; keep <- colnames(Bs); log_msg("  slim бет (поток) ", nrow(Bs), " x ", ncol(Bs))
    bs <- as.numeric(Bs[sample(nrow(Bs), min(2000, nrow(Bs))), ]); bs <- sample(bs[is.finite(bs)], min(20000, sum(is.finite(bs))))
    pr <- predict_sex_meth(if (!is.null(sb$chrY_beta)) sb$chrY_beta else Bs, sex[keep])
    meta[gse == g & gsm %in% keep, Sex := pr$sex[match(gsm, keep)]]
    meanY <- if (!is.null(sb$chrY_beta)) colMeans(sb$chrY_beta, na.rm = TRUE) else rep(NA_real_, length(keep))
    sexrow <- data.frame(gse = g, meanY = meanY, Sex = pr$sex, predicted = pr$pred)
    cellprop <- if (!is.null(sb$ref_beta)) tryCatch(deconv_blood_meth(sb$ref_beta), error = function(e) NULL) else NULL
    M <- beta_to_m(Bs); rm(Bs, sb); gc()
    if (anyNA(M) && requireNamespace("impute", quietly = TRUE)) M <- impute::impute.knn(M, k = 10)$data
    log_msg("  slim M ", nrow(M), " x ", ncol(M), " | пол предсказан у ", sum(pr$pred))
    return(list(gse = g, M = M, age = age[keep], sex = pr$sex[match(keep, names(pr$sex))],
                cellprop = cellprop, beta_sample = bs, sexrow = sexrow))
  }

  if (spec$type == "series") {
    B <- read_series_betas(files[1], keep)
  } else {
    rd <- function(f) { cm <- map_keep_to_filecols(keep, beta_colnames(f), meta_g)
      if (!length(cm)) stop("колонки бет не сопоставлены с GSM в ", basename(f), " (см. peek_gz)")
      log_msg("  читаю ", basename(f), " (", length(cm), " образцов)..."); read_betas_file(f, cm) }
    B <- if (spec$type == "suppl2") combine_two(rd(files[1]), rd(files[2])) else rd(files[1])
  }
  storage.mode(B) <- "double"
  keep <- intersect(keep, colnames(B)); B <- B[, keep, drop = FALSE]
  log_msg("  матрица бет ", nrow(B), " x ", ncol(B))

  # подвыборка β для графика (полное распределение, ДО урезаний)
  ridx <- sample(nrow(B), min(2000, nrow(B)))
  bs <- as.numeric(B[ridx, ]); bs <- sample(bs[is.finite(bs)], min(20000, sum(is.finite(bs))))

  # предсказание пола по chrY (ДО удаления половых хромосом!)
  pr <- predict_sex_meth(B, sex[colnames(B)])
  meta[gse == g & gsm %in% colnames(B), Sex := pr$sex[match(gsm, colnames(B))]]
  yp <- intersect(chrY_probes_450k(), rownames(B))
  meanY <- if (length(yp)) colMeans(B[yp, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(B))
  sexrow <- data.frame(gse = g, meanY = meanY, Sex = pr$sex, predicted = pr$pred)

  # деконволюция крови по полным бетам (EpiDISH RPC) — доли клеток как ковариаты для ComBat/трендов
  cellprop <- tryCatch(deconv_blood_meth(B), error = function(e){ log_msg("  EpiDISH не вышел: ", conditionMessage(e)); NULL })

  # фильтр зондов: половые хромосомы (+ кросс-реактивные), доля NA
  drop <- union(sex_chr_probes_450k(), xreactive_probes_450k())
  if (length(drop)) B <- B[!rownames(B) %in% drop, , drop = FALSE]
  B <- B[rowMeans(is.na(B)) <= NA_PROBE_MAX, , drop = FALSE]

  # топ-N вариабельных отбираем НА БЕТАХ, и лишь ПОТОМ beta->M на слиме:
  # так пик памяти = одна полная матрица (важно для GSE55763 даже на Colab).
  v <- matrixStats::rowVars(B, na.rm = TRUE)
  B <- B[head(order(v, decreasing = TRUE), N_VAR_CPG), , drop = FALSE]
  M <- beta_to_m(B); rm(B); gc()
  if (anyNA(M) && requireNamespace("impute", quietly = TRUE)) M <- impute::impute.knn(M, k = 10)$data
  log_msg("  slim M ", nrow(M), " x ", ncol(M), " | пол предсказан у ", sum(pr$pred),
          if (is.null(cellprop)) " | без деконволюции" else "")

  list(gse = g, M = M, age = age[colnames(M)], sex = pr$sex[match(colnames(M), names(pr$sex))],
       cellprop = cellprop, beta_sample = bs, sexrow = sexrow)
}

for (g in targets) {
  log_msg("=== ", g, " ===")
  r <- tryCatch(process_meth(g), error = function(e){ log_msg("Ошибка ", g, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) save_rds(r, file.path(PATHS$interim_meth, paste0("04_", g, ".rds")))
  rm(r); gc()
}
# обновлённые метаданные (с предсказанным полом)
save_rds(meta, file.path(PATHS$interim_meth, "02_meta_meth.rds"))

# --- графики из сохранённых slim ---
files <- list.files(PATHS$interim_meth, pattern = "^04_GSE.*\\.rds$", full.names = TRUE)
if (length(files)) {
  obj <- lapply(files, readRDS)
  nper <- setNames(sapply(obj, function(x) ncol(x$M)), sapply(obj, function(x) x$gse))
  df  <- do.call(rbind, lapply(obj, function(x) data.frame(dataset = x$gse, value = x$beta_sample)))
  xy  <- "Methylation level, \u03b2-value"
  make_ridgeline(df, METH_COL, xlab = xy, title = "Methylation by dataset (before ComBat)",
                 file = file.path(PATHS$pics, "04_meth_ridgeline_before.png"), counts = nper)
  make_violin   (df, METH_COL, ylab = xy, title = "Methylation by dataset (before ComBat)",
                 file = file.path(PATHS$pics, "04_meth_violin_before.png"), counts = nper)
  sd <- do.call(rbind, lapply(obj, function(x) x$sexrow))
  p <- ggplot(sd, aes(meanY, gse, colour = Sex, shape = predicted)) +
    geom_jitter(height = .2, alpha = .6, size = 1.3) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c("reported","predicted"), name = NULL) +
    scale_colour_manual(values = c(M = "#3b6ea5", F = "#b5495b")) +
    labs(title = "Methylation: sex assignment (mean chrY \u03b2)",
         x = "Mean chrY \u03b2-value", y = NULL) + theme_soda()
  ggsave(file.path(PATHS$pics, "04_sex_prediction_methylation.png"), p, width = 7, height = 5, dpi = 150)
  log_msg("04 готов.")
} else log_msg("04: slim-файлы не найдены.")
