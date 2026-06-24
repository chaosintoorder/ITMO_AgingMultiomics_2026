# ============================================================================
# R/utils_io.R — тяжёлый ввод-вывод, омикс-хелперы и общие графики.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2) })

# ---------- ПОИСК / ЧТЕНИЕ БОЛЬШИХ МАТРИЦ ----------------------------------
find_file <- function(fname, dir) {
  p <- file.path(dir, fname); if (file.exists(p)) return(p)
  h <- list.files(dir, pattern = paste0("^", gsub("([.()|])","\\\\\\1",fname), "$"),
                  full.names = TRUE, recursive = TRUE)
  if (length(h)) h[1] else stop("файл не найден: ", fname, " (", dir, ")")
}
peek_gz <- function(gz, n = 3) {
  con <- gzfile(gz, "rt"); on.exit(close(con))
  ln <- gsub("\r","", readLines(con, n = n), fixed = TRUE)
  sep <- if (grepl("\t", ln[1])) "\t" else ","
  h <- gsub('"',"", strsplit(ln[1], sep, fixed = TRUE)[[1]])
  cat(sprintf("  %s | колонок=%d | первые имена: %s\n", basename(gz), length(h),
              paste(head(h, 6), collapse=" | ")))
  invisible(list(ncol = length(h), header = h))
}
cap_samples <- function(gsm, age, cap) {
  if (length(gsm) <= cap) return(gsm)
  b <- cut(age, breaks = seq(0,110,10), include.lowest = TRUE)
  idx <- unlist(tapply(seq_along(gsm), b, function(ix){
    k <- max(1L, round(cap*length(ix)/length(gsm))); if (length(ix)<=k) ix else sample(ix,k)
  }), use.names = FALSE)
  gsm[sort(idx)]
}
beta_colnames <- function(f) names(data.table::fread(f, nrows = 0))

# сопоставить нужные GSM с реальными колонками файла.
# Если колонки не GSM — мапим через ЛЮБОЕ текстовое поле метаданных (title/source/characteristics).
map_keep_to_filecols <- function(keep_gsm, file_cols, meta) {
  # 1) колонки уже GSM
  if (any(startsWith(file_cols, "GSM"))) { s <- intersect(keep_gsm, file_cols); return(setNames(s, s)) }
  # данные-колонки: всё, кроме ID-колонки (1-я) и Detection Pval
  is_pv <- grepl("detection|p[-_ ]?val", file_cols, ignore.case = TRUE)
  data_cols <- file_cols[-1][!is_pv[-1]]
  # 2) маппинг по любому текстовому полю метаданных (title/source/characteristics)
  lut <- character(0)
  for (cn in setdiff(names(meta), "gsm")) {
    v <- as.character(meta[[cn]]); ok <- !is.na(v) & nzchar(v)
    lut <- c(lut, setNames(meta$gsm[ok], v[ok]))
  }
  lut <- lut[!duplicated(names(lut))]
  g <- lut[data_cols]; hit <- !is.na(g)
  if (any(hit)) { gg <- g[hit]; nm <- data_cols[hit]; ok <- gg %in% keep_gsm
    return(setNames(unname(gg[ok]), nm[ok])) }
  # 3) позиционный маппинг: порядок колонок файла = порядок образцов в метаданных
  if (length(data_cols) == nrow(meta)) {
    log_msg("  ВНИМАНИЕ: колонки бет не подписаны GSM — позиционный маппинг (тот же порядок образцов)")
    ok <- meta$gsm %in% keep_gsm
    return(setNames(meta$gsm[ok], data_cols[ok]))
  }
  character(0)
}
read_series_betas <- function(gz, keep_gsm = NULL) {
  dt <- data.table::fread(gz, skip = "ID_REF", header = TRUE, showProgress = FALSE)
  idc <- names(dt)[1]; dt <- dt[!grepl("^!", get(idc))]
  if (!is.null(keep_gsm)) dt <- dt[, c(idc, intersect(keep_gsm, names(dt))), with = FALSE]
  m <- as.matrix(dt[, -1L, with = FALSE]); rownames(m) <- dt[[idc]]; m
}
read_betas_file <- function(f, cols_map) {
  hdr <- beta_colnames(f); idc <- hdr[1]
  dt <- data.table::fread(f, select = intersect(c(idc, names(cols_map)), hdr), showProgress = FALSE)
  m <- as.matrix(dt[, -1L, with = FALSE]); rownames(m) <- dt[[idc]]
  colnames(m) <- cols_map[colnames(m)]; m
}
# ПОТОКОВОЕ чтение огромного gz (GSE55763, ~24 ГБ распакованных): без отображения в память.
# Память ~ выходная матрица (зонды × выбранные образцы). Медленно, но влезает в 8 ГБ.
read_betas_stream <- function(file, keep_gsm, meta, chunk = 20000, n_alloc = 6e5) {
  con <- gzfile(file, "rt"); on.exit(close(con))
  hdr <- gsub('"',"", strsplit(gsub("\r","",readLines(con, 1L)), "\t", fixed = TRUE)[[1]])
  cm  <- map_keep_to_filecols(keep_gsm, hdr, meta)
  if (!length(cm)) stop("колонки бет не сопоставлены с GSM (stream) — глянь peek_gz()")
  ci  <- match(names(cm), hdr); new_names <- unname(cm[hdr[ci]])
  log_msg("  поток: ", length(ci), " колонок, читаю чанками по ", chunk, " строк...")
  B <- matrix(NA_real_, n_alloc, length(ci)); rid <- character(n_alloc); r <- 0L
  repeat {
    ln <- readLines(con, n = chunk); if (!length(ln)) break
    pr <- strsplit(gsub("\r","",ln), "\t", fixed = TRUE)
    pr <- pr[lengths(pr) >= max(ci)]; k <- length(pr); if (!k) next
    sub <- do.call(rbind, lapply(pr, `[`, ci))
    B[(r+1L):(r+k), ] <- suppressWarnings(as.numeric(sub))
    rid[(r+1L):(r+k)]  <- vapply(pr, `[`, character(1), 1L)
    r <- r + k
    if (r %% 100000L < chunk) log_msg("    ...", r, " строк")
  }
  B <- B[seq_len(r), , drop = FALSE]; rownames(B) <- rid[seq_len(r)]; colnames(B) <- new_names
  B
}
combine_two <- function(a, b) {
  if (length(intersect(colnames(a), colnames(b))) > 0.5*min(ncol(a),ncol(b))) {
    cc <- intersect(colnames(a), colnames(b)); rbind(a[,cc,drop=FALSE], b[,cc,drop=FALSE])
  } else { rr <- intersect(rownames(a), rownames(b)); cbind(a[rr,,drop=FALSE], b[rr,,drop=FALSE]) }
}
beta_to_m <- function(b) { b[b<=0] <- 1e-6; b[b>=1] <- 1-1e-6; log2(b/(1-b)) }
m_to_beta <- function(m) 2^m/(1+2^m)

# ---------- метаданные метилирования (series / xlsx / сеть) ----------------
read_meth_meta <- function(gse, dir) {
  sm <- list.files(dir, pattern = paste0("^",gse,"_series_matrix.*\\.txt\\.gz$"), full.names = TRUE)
  if (length(sm)) return(read_geo_meta(sm[1]))
  xl <- list.files(dir, pattern = paste0("^",gse,".*\\.xlsx$"), full.names = TRUE)
  if (length(xl) && requireNamespace("readxl", quietly = TRUE)) {
    df <- as.data.frame(readxl::read_excel(xl[1])); names(df)[1] <- "gsm"; return(df)
  }
  es <- tryCatch(GEOquery::getGEO(gse, GSEMatrix=TRUE, getGPL=FALSE, destdir=dir)[[1]],
                 error = function(e) NULL)
  if (!is.null(es)) { m <- Biobase::pData(es); m$gsm <- rownames(m); return(m) }
  stop("нет источника метаданных для ", gse, " — положи его *_series_matrix.txt.gz (даже маленький) в ", dir)
}

# ---------- аннотация 450k (с приатачиванием пакета — иначе minfi падает) ---
.anno_cache <- new.env()
.anno450k <- function() {
  if (!is.null(.anno_cache$an)) return(.anno_cache$an)
  ns <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  ok <- requireNamespace("minfi", quietly = TRUE) && requireNamespace(ns, quietly = TRUE)
  if (!ok) return(NULL)
  suppressPackageStartupMessages(require(ns, character.only = TRUE))   # ATTACH: нужно minfi::getAnnotation
  .anno_cache$an <- minfi::getAnnotation(get(ns)); .anno_cache$an
}
sex_chr_probes_450k <- function() { an <- .anno450k(); if (is.null(an)) character(0)
  else rownames(an)[an$chr %in% c("chrX","chrY")] }
chrY_probes_450k <- function() { an <- .anno450k(); if (is.null(an)) character(0)
  else rownames(an)[an$chr == "chrY"] }
xreactive_probes_450k <- function() {
  if (requireNamespace("maxprobes", quietly=TRUE))
    unlist(maxprobes::xreactive_probes(array_type="450K"), use.names=FALSE) else character(0)
}

# ---------- ЭКСПРЕССИЯ: лог, аннотация, схлопывание, пол -------------------
is_log2 <- function(X) {              # эвристика NCBI GEO2R
  qx <- as.numeric(quantile(X, c(0,.25,.5,.75,.99,1.0), na.rm = TRUE))
  !( (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) )
}
collapse_to_gene <- function(X, sym) {
  ok <- !is.na(sym) & !sym %in% c("","---","NA"); X <- X[ok,,drop=FALSE]; sym <- sym[ok]
  if (!nrow(X)) return(X)
  o <- order(rowMeans(X, na.rm=TRUE), decreasing=TRUE); X <- X[o,,drop=FALSE]; sym <- sym[o]
  d <- !duplicated(sym); X <- X[d,,drop=FALSE]; rownames(X) <- sym[d]; X
}
annotate_to_genes <- function(X, gse) {
  gpl <- EXPR_PLATFORM[[gse]]; pkg <- BIOC_ANNO[[gpl]]; sym <- NULL
  if (!is.null(pkg) && !is.na(pkg) && requireNamespace(pkg, quietly = TRUE)) {     # 1) BioC-пакет
    sym <- tryCatch(AnnotationDbi::mapIds(getExportedValue(pkg, pkg), keys = rownames(X),
                    column = "SYMBOL", keytype = "PROBEID", multiVals = "first"),
                    error = function(e){ log_msg("  ", gse, ": mapIds fail: ", conditionMessage(e)); NULL })
  }
  if (is.null(sym)) {                                                              # 2) таблица GPL из GEO
    gp <- tryCatch(GEOquery::getGEO(gpl, destdir = PATHS$raw_expr), error = function(e) NULL)
    if (!is.null(gp)) {
      tb <- GEOquery::Table(gp)
      log_msg("  ", gse, " (", gpl, ") колонки GPL: ", paste(names(tb), collapse = ", "))
      sc <- grep("gene[_ .]?symbol|^symbol$|ilmn_gene|^gene$", names(tb), ignore.case = TRUE, value = TRUE)
      if (length(sc)) {
        s <- as.character(tb[[sc[1]]]); names(s) <- as.character(tb[[1]])
        sym <- trimws(sub("\\s*//.*$","", s[rownames(X)]))
      } else {                                                                     # 3) RefSeq/GenBank -> символ
        rc <- grep("refseq|genbank|gb_acc|gb_list|accession", names(tb), ignore.case = TRUE, value = TRUE)
        if (length(rc) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
          acc <- sub("\\..*$","", as.character(tb[[rc[1]]])); names(acc) <- as.character(tb[[1]])
          m <- tryCatch(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = acc[rownames(X)],
                        column = "SYMBOL", keytype = "REFSEQ", multiVals = "first"),
                        error = function(e) NULL)
          if (!is.null(m)) sym <- unname(m)
        }
      }
    }
  }
  if (is.null(sym)) { log_msg("  ", gse, ": символы генов не получены — оставляю probe ID")
    attr(X, "annotated") <- FALSE; return(X) }
  out <- collapse_to_gene(X, as.character(sym))
  if (!nrow(out)) { log_msg("  ", gse, ": аннотация дала 0 генов — оставляю probe ID")
    attr(X, "annotated") <- FALSE; return(X) }
  attr(out, "annotated") <- TRUE; out
}
predict_sex_expr <- function(X, sex) {
  g <- rownames(X)
  f <- if ("XIST" %in% g) X["XIST",] else rep(NA_real_, ncol(X))
  yg <- intersect(SEX_MARKERS$male, g)
  m <- if (length(yg)) colMeans(X[yg,,drop=FALSE], na.rm=TRUE) else rep(NA_real_, ncol(X))
  if (all(is.na(f))||all(is.na(m))) return(list(sex = sex, pred = rep(FALSE, ncol(X))))
  pred <- ifelse(scale(m)[,1] >= scale(f)[,1], "M", "F")
  out <- sex; miss <- is.na(out); out[miss] <- pred[miss]
  list(sex = out, pred = miss)
}
predict_sex_meth <- function(B, sex) {            # средняя бета chrY + 2 кластера
  yp <- intersect(chrY_probes_450k(), rownames(B))
  if (length(yp) < 5) return(list(sex = sex, pred = rep(FALSE, ncol(B))))
  ym <- colMeans(B[yp,,drop=FALSE], na.rm=TRUE)
  km <- stats::kmeans(ym, centers = 2)
  male_cl <- which.max(tapply(ym, km$cluster, mean))
  pred <- ifelse(km$cluster == male_cl, "M", "F")
  out <- sex; miss <- is.na(out); out[miss] <- pred[miss]
  list(sex = out, pred = miss)
}

# ---------- ОБЩИЕ ГРАФИКИ --------------------------------------------------
.darker <- function(cols) vapply(cols, function(c) grDevices::adjustcolor(c, red.f=.8, green.f=.8, blue.f=.8),
                                 character(1))
# подписи "(n=...) GSE": counts — это число ОБРАЗЦОВ на датасет (не число точек графика!)
.labels <- function(lv, df, counts) {
  cnt <- if (is.null(counts)) as.integer(table(factor(df$dataset, levels = lv))[lv]) else as.integer(counts[lv])
  setNames(sprintf("(n=%d) %s", cnt, lv), lv)
}
make_ridgeline <- function(df, pal, xlab, title, file, counts = NULL) {
  stopifnot(requireNamespace("ggridges", quietly = TRUE))
  lv <- intersect(names(pal), unique(df$dataset)); lab <- .labels(lv, df, counts)
  df$lab <- factor(lab[as.character(df$dataset)], levels = lab[lv])
  p <- ggplot(df, aes(value, lab, fill = dataset, colour = dataset)) +
    ggridges::geom_density_ridges(alpha = .55, quantile_lines = TRUE, quantiles = .5,
                                  linewidth = .4, scale = 1.1) +
    scale_fill_manual(values = pal) + scale_colour_manual(values = .darker(pal)) +
    labs(title = title, x = xlab, y = NULL) + theme_soda() + theme(legend.position = "none")
  save_plot(p, sub("\\.png$", "", file), w = 8, h = 6)
}
make_violin <- function(df, pal, ylab, title, file, counts = NULL) {
  lv <- intersect(names(pal), unique(df$dataset)); lab <- .labels(lv, df, counts)
  df$lab <- factor(lab[as.character(df$dataset)], levels = lab[lv])
  p <- ggplot(df, aes(lab, value, fill = dataset, colour = dataset)) +
    geom_violin(alpha = .55, linewidth = .4, scale = "width") +
    stat_summary(fun = median, geom = "crossbar", width = .5, linewidth = .3, colour = "grey20") +
    scale_fill_manual(values = pal) + scale_colour_manual(values = .darker(pal)) +
    labs(title = title, x = NULL, y = ylab) + theme_soda() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, sub("\\.png$", "", file), w = 9, h = 5)
}
# распределение ВОЗРАСТА: ступенчатый (гистограммный) ridge + суммарная строка "All"
make_age_ridgeline <- function(df, pal, xlab, title, file, bins = 30) {
  stopifnot(requireNamespace("ggridges", quietly = TRUE))
  lv  <- intersect(names(pal), unique(df$dataset))
  df2 <- rbind(df[, c("dataset","value")], data.frame(dataset = "All", value = df$value))
  pal2 <- c(pal[lv], All = "#555555")
  cnt  <- sapply(c("All", lv), function(d) sum(df2$dataset == d))
  lab  <- setNames(sprintf("(n=%d) %s", cnt, c("All", lv)), c("All", lv))
  df2$lab <- factor(lab[df2$dataset], levels = lab[c("All", lv)])   # All — внизу
  p <- ggplot(df2, aes(value, lab, fill = dataset, colour = dataset)) +
    ggridges::geom_density_ridges(stat = "binline", bins = bins, scale = 0.95,
                                  alpha = .6, linewidth = .35, draw_baseline = FALSE) +
    scale_fill_manual(values = pal2) + scale_colour_manual(values = .darker(pal2)) +
    labs(title = title, x = xlab, y = NULL) + theme_soda() + theme(legend.position = "none")
  ggsave(file, p, width = 8, height = 6.5, dpi = 150); log_msg("  график: ", file)
}

# ---------- ДЕКОНВОЛЮЦИЯ СОСТАВА КРОВИ -------------------------------------
# метилирование: EpiDISH RPC, референс DHS-blood (B/NK/CD4T/CD8T/Mono/Neutro/Eosino)
deconv_blood_meth <- function(B) {
  if (!requireNamespace("EpiDISH", quietly = TRUE)) return(NULL)
  e <- new.env(); utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
  EpiDISH::epidish(beta.m = B, ref.m = e$centDHSbloodDMC.m, method = "RPC")$estF  # образцы × клетки
}
# экспрессия: маркерные скоры (стиль MCPcounter). Те же 7 типов, что и EpiDISH-референс
# (B/CD4T/CD8T/NK/Mono/Neutro/Eosino) — чтобы популяции совпадали с метилированием.
EXPR_CELL_MARKERS <- list(
  B      = c("CD19","MS4A1","CD79A","CD79B"),
  CD4T   = c("CD4","IL7R","CD40LG","CD3D"),
  CD8T   = c("CD8A","CD8B","GZMK"),
  NK     = c("NCAM1","KLRD1","NKG7","GNLY"),
  Mono   = c("CD14","LYZ","CSF1R","FCN1"),
  Neutro = c("FCGR3B","CSF3R","S100A8","S100A9"),
  Eosino = c("CCR3","IL5RA","PRG2","CLC","SIGLEC8"))
deconv_blood_expr <- function(X) {                # X: гены × образцы (log2)
  z <- t(scale(t(X)))                             # z-скор по образцам на ген
  m <- sapply(EXPR_CELL_MARKERS, function(mk) {
    g <- intersect(mk, rownames(z)); if (!length(g)) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
  rownames(m) <- colnames(X); m[, colSums(is.na(m)) < nrow(m), drop = FALSE]
}
# хитмэп клеточного состава: строки = типы клеток, столбцы = ВОЗРАСТНЫЕ ДЕКАДЫ
# (медиана доли), фасет по датасету. Без тысяч образцов-полосок -> читаемые блоки.
cellcomp_heatmap <- function(cp, meta, palette, title, value_label, file_noext) {
  age <- meta$Age[match(rownames(cp), meta$gsm)]; gse <- meta$gse[match(rownames(cp), meta$gsm)]
  bin <- cut(age, breaks = AGE_BIN_BREAKS, include.lowest = TRUE)
  keep <- !is.na(bin) & !is.na(gse)
  df <- do.call(rbind, lapply(colnames(cp), function(k)
    data.frame(cell = k, gse = gse[keep], bin = bin[keep], value = cp[keep, k])))
  ag <- aggregate(value ~ cell + gse + bin, df, median)               # медиана по (тип, датасет, декада)
  ag$cell <- factor(ag$cell, levels = rev(colnames(cp)))
  ag$bin  <- factor(ag$bin, levels = levels(bin))
  p <- ggplot(ag, aes(bin, cell, fill = value)) + geom_tile(colour = "white", linewidth = .3) +
    facet_grid(~ gse, scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours = palette, name = value_label) +
    labs(title = title, x = "Age, decade", y = NULL) +
    theme_soda() + theme(panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8),
                         legend.position = "right")
  save_plot(p, file_noext, w = 12, h = 4.6)
}
# собрать матрицу ковариат под порядок колонок матрицы; NA-колонки выкинуть, прочие NA -> среднее
build_covars <- function(cols, meta, cellprop = NULL) {
  i <- match(cols, meta$gsm)
  d <- data.frame(Age = as.numeric(meta$Age[i]), Sex = factor(meta$Sex[i]))
  if (anyNA(d$Sex)) { lv <- names(which.max(table(d$Sex))); d$Sex[is.na(d$Sex)] <- lv; d$Sex <- factor(d$Sex) }
  if (!is.null(cellprop)) {
    cp <- cellprop[match(cols, rownames(cellprop)), , drop = FALSE]
    cp <- cp[, colSums(is.na(cp)) < nrow(cp), drop = FALSE]
    for (j in seq_len(ncol(cp))) cp[is.na(cp[,j]), j] <- mean(cp[,j], na.rm = TRUE)
    # композиционные доли (EpiDISH) суммируются в 1 -> коллинеарны с интерсептом:
    # выкидываем одну колонку-референс, иначе ComBat ругается "confounded with batch"
    if (ncol(cp) > 2 && stats::median(rowSums(cp), na.rm = TRUE) > 0.95 &&
                        stats::median(rowSums(cp), na.rm = TRUE) < 1.05)
      cp <- cp[, -which.min(colMeans(cp)), drop = FALSE]
    d <- cbind(d, cp)
  }
  if (nlevels(d$Sex) < 2) d$Sex <- NULL          # один пол в выборке -> ковариата не нужна
  d$Age[is.na(d$Age)] <- stats::median(d$Age, na.rm = TRUE)
  d
}

# ---------- ComBat + PCA-диагностика batch --------------------------------
run_combat <- function(mat, batch, covdf) {
  batch <- factor(batch); ref <- names(which.max(table(batch)))   # крупнейший датасет = опорный
  mm <- stats::model.matrix(~ ., data = covdf)
  log_msg("  ComBat: ", nrow(mat), " фич × ", ncol(mat), " образцов | батчей ", nlevels(batch),
          " | опорный ", ref)
  sva::ComBat(dat = mat, batch = batch, mod = mm, ref.batch = ref, par.prior = TRUE)
}
pca_batch_plot <- function(mat, batch, title, file_noext, age = NULL, nfeat = N_PCA_FEATURES) {
  v <- matrixStats::rowVars(mat, na.rm = TRUE)
  idx <- head(order(v, decreasing = TRUE), min(nfeat, nrow(mat)))
  pc <- stats::prcomp(t(mat[idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  d  <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], batch = factor(batch))
  pb <- ggplot(d, aes(PC1, PC2, colour = batch)) + geom_point(alpha = .5, size = 1.1) +
    labs(title = title, subtitle = "coloured by dataset (batch)",
         x = paste0("PC1, ", ve[1], "%"), y = paste0("PC2, ", ve[2], "%")) +
    theme_soda() + theme(legend.position = "right")
  save_plot(pb, paste0(file_noext, "_batch"), w = 7.8, h = 5.6)
  if (!is.null(age)) {                                  # биологический градиент (должен сохраниться)
    d$age <- as.numeric(age)
    pa <- ggplot(d, aes(PC1, PC2, colour = age)) + geom_point(alpha = .6, size = 1.1) +
      scale_colour_viridis_c(option = "C", name = "Age, years") +
      labs(title = title, subtitle = "coloured by age (biology, preserved)",
           x = paste0("PC1, ", ve[1], "%"), y = paste0("PC2, ", ve[2], "%")) + theme_soda()
    save_plot(pa, paste0(file_noext, "_age"), w = 7.8, h = 5.6)
  }
  if (requireNamespace("uwot", quietly = TRUE)) {        # UMAP (нелинейная проекция), если есть пакет
    um <- uwot::umap(pc$x[, 1:min(20, ncol(pc$x))], n_neighbors = 15, min_dist = 0.3)
    du <- data.frame(U1 = um[,1], U2 = um[,2], batch = factor(batch))
    pu <- ggplot(du, aes(U1, U2, colour = batch)) + geom_point(alpha = .5, size = 1.1) +
      labs(title = title, subtitle = "UMAP, coloured by dataset",
           x = "UMAP-1", y = "UMAP-2") + theme_soda() + theme(legend.position = "right")
    save_plot(pu, paste0(file_noext, "_umap"), w = 7.8, h = 5.6)
  }
}

# ---------- ДВУХПРОХОДНОЕ потоковое слим-чтение огромного gz (GSE55763) -----
# Проход 1: дисперсия каждого CpG (храним только числа) + строки chrY и референса EpiDISH.
# Проход 2: читаем ТОЛЬКО топ-N вариабельных CpG. Пик памяти ~ slim-матрица (а не вся).
read_betas_slim_stream <- function(file, keep_gsm, meta, n_top, chunk = 8000) {
  chrY <- chrY_probes_450k(); drop_set <- union(sex_chr_probes_450k(), xreactive_probes_450k())
  refc <- tryCatch({ e <- new.env(); utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
                     rownames(e$centDHSbloodDMC.m) }, error = function(e) character(0))
  con <- gzfile(file, "rt")
  hdr <- gsub('"',"", strsplit(gsub("\r","",readLines(con, 1L)), "\t", fixed = TRUE)[[1]])
  cm  <- map_keep_to_filecols(keep_gsm, hdr, meta)
  if (!length(cm)) { close(con); stop("колонки бет не сопоставлены (slim stream) — глянь peek_gz()") }
  ci <- match(names(cm), hdr); new_names <- unname(cm[hdr[ci]])
  parse <- function(pr) { suppressWarnings(matrix(as.numeric(do.call(rbind, lapply(pr, `[`, ci))),
                                                  nrow = length(pr))) }
  # ---- проход 1 ----
  log_msg("  поток 1/2: дисперсии CpG, ", length(ci), " образцов...")
  vn <- 6e5; vid <- character(vn); vv <- numeric(vn); k <- 0L; chrYr <- list(); refr <- list(); seen <- 0L
  repeat {
    ln <- readLines(con, n = chunk); if (!length(ln)) break
    pr <- strsplit(gsub("\r","",ln), "\t", fixed = TRUE); pr <- pr[lengths(pr) >= max(ci)]; if (!length(pr)) next
    ids <- vapply(pr, `[`, character(1), 1L); mat <- parse(pr); rownames(mat) <- ids
    iy <- ids %in% chrY; if (any(iy)) chrYr[[length(chrYr)+1L]] <- mat[iy, , drop = FALSE]
    ir <- ids %in% refc; if (any(ir)) refr[[length(refr)+1L]]  <- mat[ir, , drop = FALSE]
    kv <- !(ids %in% drop_set)
    if (any(kv)) { mm <- mat[kv, , drop = FALSE]; vr <- matrixStats::rowVars(mm, na.rm = TRUE)
      nk <- length(vr); vid[(k+1L):(k+nk)] <- rownames(mm); vv[(k+1L):(k+nk)] <- vr; k <- k + nk }
    seen <- seen + length(pr); if (seen %% 100000L < chunk) log_msg("    ...", seen, " CpG")
  }
  close(con)
  top <- vid[seq_len(k)][head(order(vv[seq_len(k)], decreasing = TRUE), n_top)]
  # ---- проход 2 ----
  log_msg("  поток 2/2: сбор ", length(top), " вариабельных CpG...")
  con <- gzfile(file, "rt"); readLines(con, 1L)
  slim <- matrix(NA_real_, length(top), length(ci), dimnames = list(top, new_names))
  pos <- setNames(seq_along(top), top)
  repeat {
    ln <- readLines(con, n = chunk); if (!length(ln)) break
    pr <- strsplit(gsub("\r","",ln), "\t", fixed = TRUE); pr <- pr[lengths(pr) >= max(ci)]; if (!length(pr)) next
    ids <- vapply(pr, `[`, character(1), 1L); hit <- ids %in% top; if (!any(hit)) next
    slim[pos[ids[hit]], ] <- parse(pr[hit])
  }
  close(con)
  cb <- function(L) if (length(L)) { m <- do.call(rbind, L); colnames(m) <- new_names; m } else NULL
  list(slim_beta = slim, chrY_beta = cb(chrYr), ref_beta = cb(refr))
}
