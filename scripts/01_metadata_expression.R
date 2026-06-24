# ============================================================================
# scripts/01_metadata_expression.R
# Вход : экспрессионные *_series_matrix.txt.gz в data/raw/expression (только шапки)
# Выход: единая таблица метаданных (gsm, gse, Age, Sex, AgeBin) -> 01_meta_expr.rds
#        + график распределения возраста по 13 датасетам.
# Пол здесь — только заявленный (M/F/NA); пропущенный пол предскажем в 03 (нужны матрицы).
# ============================================================================
source("config.R")

meta <- list()
for (g in EXPR_GSE) {
  f <- list.files(PATHS$raw_expr, pattern = paste0("^",g,"_series_matrix.*\\.txt\\.gz$"), full.names = TRUE)
  if (!length(f)) { log_msg("ПРОПУСК ", g, ": нет серматрицы"); next }
  m  <- read_geo_meta(f[1])
  dt <- data.table::data.table(gsm = m$gsm, gse = g, Age = extract_age(m), Sex = extract_sex(m))
  meta[[g]] <- dt
  log_msg(sprintf("%s: n=%d | возраст %s | M=%d F=%d NA_sex=%d", g, nrow(dt),
          paste(suppressWarnings(round(range(dt$Age, na.rm=TRUE))), collapse="-"),
          sum(dt$Sex=="M",na.rm=TRUE), sum(dt$Sex=="F",na.rm=TRUE), sum(is.na(dt$Sex))))
}
meta <- data.table::rbindlist(meta)

# возраст — целевая переменная: образцы без Age удаляем (и пометим для матриц)
n0 <- nrow(meta); meta <- meta[!is.na(Age)]
log_msg("Удалено образцов без возраста: ", n0 - nrow(meta), " из ", n0)
meta[, AgeBin := age_bin(Age)]                    # только для визуализации (PCA), не ковариата

save_rds(meta, file.path(PATHS$interim_expr, "01_meta_expr.rds"))
data.table::fwrite(meta, file.path(PATHS$tables, "01_meta_expr.csv"))

# график: распределение возраста по датасетам
df <- data.frame(dataset = meta$gse, value = meta$Age)
make_age_ridgeline(df, EXPR_COL, xlab = "Age, years",
                   title = "Expression cohorts: age distribution",
                   file = file.path(PATHS$pics, "01_age_distribution_expression.png"))
log_msg("01 готов.")
