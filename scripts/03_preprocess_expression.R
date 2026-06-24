# ============================================================================
# scripts/03_preprocess_expression.R
# Вход : экспрессионные серматрицы + 01_meta_expr.rds
# Выход: объединённая матрица (гены × образцы, log2 + quantile-норм) + метаданные
#        (с предсказанным полом) -> 03_expr_processed.rds
#        графики: UpSet пересечения генов, ridgeline, violin, предсказание пола.
# Память: матрицы экспрессии небольшие; держим список из 13 — это ок на 8 ГБ.
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(GEOquery); library(Biobase); library(ggplot2) })

meta <- readRDS(file.path(PATHS$interim_expr, "01_meta_expr.rds"))

read_expr_matrix <- function(g) {
  f <- list.files(PATHS$raw_expr, pattern = paste0("^",g,"_series_matrix.*\\.txt\\.gz$"), full.names = TRUE)
  if (!length(f)) return(NULL)
  es <- getGEO(filename = f[1], getGPL = FALSE)
  if (is.list(es) && !is(es, "ExpressionSet")) es <- es[[which.max(sapply(es, ncol))]]
  X <- exprs(es)
  keep <- intersect(colnames(X), meta[gse == g, gsm]); X <- X[, keep, drop = FALSE]
  if (!ncol(X) || !nrow(X)) { log_msg("  ", g, ": нет образцов с возрастом / пустая матрица"); return(NULL) }
  X[!is.finite(X)] <- NA
  X <- X[rowMeans(is.na(X)) <= NA_PROBE_MAX, , drop = FALSE]
  if (!nrow(X)) { log_msg("  ", g, ": все зонды отфильтрованы по NA"); return(NULL) }
  if (anyNA(X) && requireNamespace("impute", quietly = TRUE)) X <- impute::impute.knn(X, k = 10)$data
  if (!is_log2(X)) { X <- log2(pmax(X, 0) + 1); log_msg("  ", g, ": применён log2") }
  annotate_to_genes(X, g)
}

res <- list(); sexdf <- list()
for (g in EXPR_GSE) {
  log_msg("=== ", g, " ===")
  X <- tryCatch(read_expr_matrix(g), error = function(e){ log_msg("Ошибка ", g, ": ", conditionMessage(e)); NULL })
  if (is.null(X) || !ncol(X)) next
  msub <- meta[gse == g]; sex0 <- setNames(msub$Sex, msub$gsm)[colnames(X)]
  pr <- predict_sex_expr(X, sex0)
  meta[gse == g & gsm %in% colnames(X), Sex := pr$sex[match(gsm, colnames(X))]]
  xist <- if ("XIST" %in% rownames(X)) X["XIST", ] else rep(NA_real_, ncol(X))
  yg <- intersect(SEX_MARKERS$male, rownames(X))
  ysc <- if (length(yg)) colMeans(X[yg, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(X))
  sexdf[[g]] <- data.frame(gse = g, XIST = xist, Yscore = ysc, Sex = pr$sex, predicted = pr$pred)
  res[[g]] <- X
  log_msg("  генов=", nrow(X), " образцов=", ncol(X), " | предсказан пол у ", sum(pr$pred))
  gc()
}
if (!length(res)) stop("ни один экспрессионный датасет не обработан.")

# аннотированные в символы генов vs оставшиеся в probe ID (их нельзя пересекать с генами)
ann <- vapply(res, function(X) isTRUE(attr(X, "annotated")), logical(1))
if (any(!ann)) log_msg("НЕ аннотированы в символы (в merge не идут): ", paste(names(res)[!ann], collapse = ", "))
res_ann <- res[ann]
if (!length(res_ann)) stop("ни один датасет не аннотирован в символы — смотри логи 'колонки GPL'.")
common <- Reduce(intersect, lapply(res_ann, rownames))

# --- таблица платформ: число фич, аннотация, вклад в общие гены ---
plat <- data.table::rbindlist(lapply(names(res), function(g) data.table::data.table(
  gse = g, platform = EXPR_PLATFORM[[g]], n_samples = ncol(res[[g]]),
  n_features = nrow(res[[g]]), annotated = isTRUE(attr(res[[g]], "annotated")),
  genes_in_common = if (isTRUE(attr(res[[g]], "annotated"))) length(intersect(rownames(res[[g]]), common)) else 0L)))
data.table::fwrite(plat, file.path(PATHS$tables, "01_expr_platforms.csv"))
log_msg("  таблица: 01_expr_platforms.csv | аннотировано ", sum(ann), "/", length(res),
        " | общих генов ", length(common))

# --- UpSet пересечения генов (только аннотированные, в палитре датасетов) ---
if (requireNamespace("UpSetR", quietly = TRUE)) {
  gl <- lapply(res_ann, rownames); sets <- names(gl)
  png(file.path(PATHS$pics, "03_upset_genes.png"), width = 1200, height = 700, res = 120)
  print(UpSetR::upset(UpSetR::fromList(gl), sets = rev(sets), keep.order = TRUE,
                      sets.bar.color = unname(EXPR_COL[rev(sets)]),
                      main.bar.color = "grey30", matrix.color = "grey30",
                      nintersects = 30, order.by = "freq"))
  dev.off(); log_msg("  график: 03_upset_genes.png")
}

# --- объединение по общим генам (БЕЗ quantile-нормализации) ---
log_msg("Общих генов: ", length(common))
merged <- do.call(cbind, lapply(res_ann, function(X) X[common, , drop = FALSE]))
# кросс-датасетную quantile-нормализацию НЕ делаем: внутри-студийная норма уже есть
# (нормализовали сабмиттеры), а кросс-студийную гармонизацию выполнит ComBat (06)
# с сохранением Age/Sex/клеток. Так распределения «до» различаются — есть что выравнивать.
batch  <- meta[match(colnames(merged), gsm), gse]
meta_m <- meta[match(colnames(merged), gsm)]

save_rds(list(expr = merged, meta = meta_m, batch = batch),
         file.path(PATHS$final_expr, "03_expr_processed.rds"))

# --- ridgeline + violin (распределения экспрессии «до» ComBat) ---
set.seed(1)
df <- do.call(rbind, lapply(unique(batch), function(g) {
  v <- as.numeric(merged[, batch == g]); v <- v[is.finite(v)]
  data.frame(dataset = g, value = sample(v, min(5000, length(v))))
}))
xy <- "Expression level, log2-scaled + normalized"
nper <- table(batch)
make_ridgeline(df, EXPR_COL, xlab = xy, title = "Expression by dataset (before ComBat)",
               file = file.path(PATHS$pics, "03_expr_ridgeline_before.png"), counts = nper)
make_violin   (df, EXPR_COL, ylab = xy, title = "Expression by dataset (before ComBat)",
               file = file.path(PATHS$pics, "03_expr_violin_before.png"), counts = nper)

# --- график предсказания пола ---
sd <- do.call(rbind, sexdf)
p <- ggplot(sd, aes(XIST, Yscore, colour = Sex, shape = predicted)) +
  geom_point(alpha = .6, size = 1.5) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c("reported", "predicted"), name = NULL) +
  scale_colour_manual(values = c(M = "#3b6ea5", F = "#b5495b")) +
  labs(title = "Expression: sex assignment (XIST vs Y-genes)",
       x = "XIST expression", y = "Mean Y-gene expression") + theme_soda()
ggsave(file.path(PATHS$pics, "03_sex_prediction_expression.png"), p, width = 7, height = 5, dpi = 150)
log_msg("03 готов.")
