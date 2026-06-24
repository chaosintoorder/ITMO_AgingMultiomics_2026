# ============================================================================
# R/utils.R — общие функции и оптимизация памяти. Источается из config.R.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# ---- оптимизация сессии (вызывается в начале каждого скрипта через config) -
setup_session <- function() {
  invisible(gc(full = TRUE))
  data.table::setDTthreads(parallel::detectCores())    # все ядра для fread
  mem_report("сессия инициализирована")
}
mem_report <- function(tag = "") {
  g <- gc(); used <- sum(g[, 2])
  log_msg("RAM в R: ", round(used), " МБ занято", if (nzchar(tag)) paste0(" | ", tag) else "")
}

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

theme_soda <- function(base = 14) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(face = "bold", size = base + 2),
                   axis.title  = ggplot2::element_text(size = base),
                   axis.text   = ggplot2::element_text(size = base - 2, colour = "grey20"),
                   legend.text = ggplot2::element_text(size = base - 2),
                   strip.text  = ggplot2::element_text(size = base - 1, face = "bold"))
}
# сохранить график СРАЗУ в png и svg (вектор для слайдов) — требование комиссии
save_plot <- function(p, path_noext, w = 8, h = 5, dpi = 200) {
  ggplot2::ggsave(paste0(path_noext, ".png"), p, width = w, height = h, dpi = dpi)
  ok <- requireNamespace("svglite", quietly = TRUE)
  ggplot2::ggsave(paste0(path_noext, ".svg"), p, width = w, height = h,
                  device = if (ok) svglite::svglite else "svg")
  log_msg("  график: ", path_noext, ".png + .svg")
}
# аккуратная подпись p-значения (без технического 0.0e+00)
fmt_pval <- function(p) {
  ifelse(is.na(p), "p = NA",
  ifelse(p < 2.2e-16, "p < 2.2e-16",
  ifelse(p < 1e-3, sprintf("p = %.1e", p), sprintf("p = %.3f", p))))
}

geo_matrix_url <- function(gse) {
  stub <- sub("\\d{3}$", "nnn", gse)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
          stub, gse, gse)
}

# ---- метаданные из локального series_matrix (НЕ грузит матрицу данных) -----
read_geo_meta <- function(gz_path) {
  con <- gzfile(gz_path, "rt"); on.exit(close(con))
  rows <- list()
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (length(l) == 0L) break
    l <- gsub("\r", "", l, fixed = TRUE)
    if (startsWith(l, "!series_matrix_table_begin")) break
    if (startsWith(l, "!Sample_")) {
      p <- strsplit(l, "\t", fixed = TRUE)[[1]]
      rows[[length(rows)+1L]] <- list(key = sub("^!","",p[1]), vals = gsub('^"|"$',"",p[-1]))
    }
  }
  acc <- unlist(lapply(rows, function(r) if (r$key=="Sample_geo_accession") r$vals))
  if (is.null(acc)) stop("нет Sample_geo_accession: ", basename(gz_path))
  out <- data.frame(gsm = acc, stringsAsFactors = FALSE)
  addf <- function(o, key, nm) { r <- Filter(function(x) x$key==key, rows)
    if (length(r) && length(r[[1]]$vals)==nrow(o)) o[[nm]] <- r[[1]]$vals; o }
  out <- addf(out, "Sample_title", "title")
  out <- addf(out, "Sample_source_name_ch1", "source")
  for (r in Filter(function(x) x$key=="Sample_characteristics_ch1", rows)) {
    v <- r$vals; if (length(v)!=nrow(out)) next
    has <- grepl(":", v)
    nm <- if (any(has)) trimws(sub(":.*$","",v[which(has)[1]])) else "char"
    nm <- make.names(ifelse(is.na(nm)||nm=="","char",nm))
    nm <- make.unique(c(names(out), nm))[length(names(out))+1L]
    out[[nm]] <- ifelse(has, trimws(sub("^[^:]*:\\s*","",v)), v)
  }
  out
}

# парсинг возраста из свободного текста: "age 11", "age: 11", "11 years/yo"
parse_age_text <- function(txt) {
  a <- suppressWarnings(as.numeric(sub(".*age[^0-9]{0,4}([0-9]{1,3}(?:\\.[0-9]+)?).*", "\\1",
                                       txt, ignore.case = TRUE, perl = TRUE)))
  need <- is.na(a)
  a[need] <- suppressWarnings(as.numeric(sub(".*?([0-9]{1,3}(?:\\.[0-9]+)?)\\s*(?:years|yrs|yo|y/o|year).*",
                                             "\\1", txt[need], ignore.case = TRUE, perl = TRUE)))
  a
}
extract_age <- function(meta) {
  cols <- setdiff(names(meta), "gsm"); age <- rep(NA_real_, nrow(meta))
  for (cn in grep("age", cols, ignore.case = TRUE, value = TRUE)) {           # 1) колонка с "age"
    v <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", meta[[cn]]))); age[is.na(age)] <- v[is.na(age)]
  }
  if (anyNA(age)) {                                                           # 2) свободный текст
    txt <- do.call(paste, c(lapply(cols, function(cn) as.character(meta[[cn]])), sep = " || "))
    a <- parse_age_text(txt); age[is.na(age)] <- a[is.na(age)]
  }
  age
}
extract_sex <- function(meta) {
  cols <- setdiff(names(meta), "gsm"); sex <- rep(NA_character_, nrow(meta))
  for (cn in grep("sex|gender", cols, ignore.case = TRUE, value = TRUE)) {    # 1) колонка sex/gender
    s <- tolower(trimws(as.character(meta[[cn]])))
    sex[is.na(sex) & (startsWith(s,"f") | grepl("female", s))] <- "F"
    sex[is.na(sex) & (startsWith(s,"m") | grepl("\\bmale\\b", s)) & !grepl("female", s)] <- "M"
  }
  if (anyNA(sex)) {                                                           # 2) свободный текст
    txt <- tolower(do.call(paste, c(lapply(cols, function(cn) as.character(meta[[cn]])), sep = " || ")))
    sex[is.na(sex) & grepl("female", txt)] <- "F"
    sex[is.na(sex) & grepl("\\bmale\\b", txt) & !grepl("female", txt)] <- "M"
  }
  sex
}
age_bin <- function(age) cut(age, breaks = AGE_BIN_BREAKS, right = FALSE,
                             labels = paste0("[", head(AGE_BIN_BREAKS,-1), ",",
                                             tail(AGE_BIN_BREAKS,-1), ")"))
window_of <- function(age, w = AGE_WINDOWS) {
  out <- rep(NA_character_, length(age))
  for (nm in names(w)) out[!is.na(age) & age>=w[[nm]][1] & age<=w[[nm]][2]] <- nm
  out
}
save_rds <- function(obj, path) { saveRDS(obj, path)
  log_msg("saved: ", path, " (", round(file.size(path)/1e6,1), " MB)") }
