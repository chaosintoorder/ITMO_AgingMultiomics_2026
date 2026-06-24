# ============================================================================
# scripts/02_metadata_methylation.R
# Вход : метаданные метилирования (из *_series_matrix.txt.gz; для GSE87571 — xlsx;
#        для GSE55763 — серматрица 49 КБ или сеть). Беты тут НЕ грузим.
# Выход: единая таблица (gsm, gse, Age, Sex, AgeBin) -> 02_meta_meth.rds + график.
# ВАЖНО: положи в data/raw/methylation маленькие серматрицы GSE87571 и GSE55763
#        (28.5/49.1 КБ) — иначе скрипт попробует xlsx/сеть как запасной вариант.
# ============================================================================
source("config.R")

meta <- list()
for (g in METH_GSE) {
  m <- tryCatch(read_meth_meta(g, PATHS$raw_meth),
                error = function(e) { log_msg("ПРОПУСК ", g, ": ", conditionMessage(e)); NULL })
  if (is.null(m)) next
  age <- extract_age(m); sex <- extract_sex(m)
  if (g %in% AGE_IN_MONTHS) { age <- age / 12; log_msg("  ", g, ": возраст из месяцев в годы") }
  dt <- data.table::data.table(gsm = m$gsm, gse = g, Age = age, Sex = sex)
  meta[[g]] <- dt
  log_msg(sprintf("%s: n=%d | возраст %s | M=%d F=%d NA_sex=%d", g, nrow(dt),
          paste(suppressWarnings(round(range(dt$Age, na.rm=TRUE))), collapse="-"),
          sum(dt$Sex=="M",na.rm=TRUE), sum(dt$Sex=="F",na.rm=TRUE), sum(is.na(dt$Sex))))
}
meta <- data.table::rbindlist(meta)

n0 <- nrow(meta); meta <- meta[!is.na(Age)]
log_msg("Удалено образцов без возраста: ", n0 - nrow(meta), " из ", n0)
meta[, AgeBin := age_bin(Age)]

save_rds(meta, file.path(PATHS$interim_meth, "02_meta_meth.rds"))
data.table::fwrite(meta, file.path(PATHS$tables, "02_meta_meth.csv"))

df <- data.frame(dataset = meta$gse, value = meta$Age)
make_age_ridgeline(df, METH_COL, xlab = "Age, years",
                   title = "Methylation cohorts: age distribution",
                   file = file.path(PATHS$pics, "02_age_distribution_methylation.png"))
log_msg("02 готов.")
