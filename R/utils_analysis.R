# ============================================================================
# R/utils_analysis.R — общие функции анализа трендов (для обеих омик).
#   Источается из config.R. Держит дорожки экспрессии/метилирования на одних
#   и тех же методах, но СКРИПТЫ запускаются раздельно (разные корзины).
# ============================================================================

# ---- (A) limma: ассоциация фичи с НЕПРЕРЫВНЫМ возрастом + BH-FDR -----------
# Ковариаты — доли клеток (внутри пола, поэтому Sex не нужен). Знак logFC = направление с возрастом.
age_limma <- function(mat, age, cellprop = NULL) {
  if (!is.null(cellprop)) {
    cp <- cellprop[match(colnames(mat), rownames(cellprop)), , drop = FALSE]
    cp <- cp[, colSums(is.na(cp)) < nrow(cp), drop = FALSE]
    for (j in seq_len(ncol(cp))) cp[is.na(cp[,j]), j] <- mean(cp[,j], na.rm = TRUE)
    # композиционные доли (EpiDISH) суммируются в 1 -> коллинеарны -> дизайн вырожден.
    # Выкидываем одну долю-референс (иначе limma: "coefficients not estimable").
    if (ncol(cp) > 2 && stats::median(rowSums(cp), na.rm = TRUE) > 0.95 &&
                        stats::median(rowSums(cp), na.rm = TRUE) < 1.05)
      cp <- cp[, -which.min(colMeans(cp)), drop = FALSE]
    design <- stats::model.matrix(~ age + cp)
  } else design <- stats::model.matrix(~ age)
  fit <- limma::eBayes(limma::lmFit(mat, design))
  tt <- limma::topTable(fit, coef = "age", number = Inf, sort.by = "none")
  tt$feature <- rownames(mat); tt$dir <- ifelse(tt$logFC > 0, "up", "down")
  tt[order(tt$adj.P.Val), ]
}

# ---- (B) GAM s(age): скорость изменения + точки разрыва -> критические возрасты ----
critical_ages <- function(mat, age, feats, omics, sx, n_breaks = 120) {
  grid <- AGE_GRID[AGE_GRID >= min(age, na.rm=TRUE) & AGE_GRID <= max(age, na.rm=TRUE)]
  eps <- 0.5; D <- matrix(NA_real_, length(grid), length(feats)); brks <- c(); fb <- list()
  for (i in seq_along(feats)) {
    y <- as.numeric(mat[feats[i], ])
    gm <- tryCatch(mgcv::gam(y ~ s(age)), error = function(e) NULL); if (is.null(gm)) next
    D[, i] <- abs((predict(gm, data.frame(age = grid + eps)) -
                   predict(gm, data.frame(age = grid - eps))) / (2*eps))
    if (i <= n_breaks && requireNamespace("segmented", quietly = TRUE)) {
      sg <- tryCatch(segmented::segmented(stats::lm(y ~ age), seg.Z = ~ age), error = function(e) NULL)
      if (!is.null(sg) && !is.null(sg$psi)) { b <- as.numeric(sg$psi[, "Est."])
        brks <- c(brks, b); fb[[feats[i]]] <- b[1] }   # точка перехода гена (для «transition genes»)
    }
  }
  rate <- rowMeans(D, na.rm = TRUE); rate <- rate / max(rate, na.rm = TRUE)
  peaks <- grid[which(diff(sign(diff(c(-Inf, rate, -Inf)))) == -2)]
  p <- ggplot2::ggplot(data.frame(age = grid, rate = rate), ggplot2::aes(age, rate)) +
    ggplot2::geom_line(colour = "#3b6ea5", linewidth = 1) +
    ggplot2::geom_vline(xintercept = peaks, linetype = "dashed", colour = "#b5495b", alpha = .7) +
    ggplot2::labs(title = paste0(omics, " (", sx, "): rate of age-related change"),
                  x = "Age, years", y = "Mean |GAM derivative|, scaled") + theme_soda()
  save_plot(p, file.path(PATHS$pics, sprintf("08_critical_ages_%s_%s", omics, sx)), w = 8, h = 4.2)
  if (length(brks)) {
    brks <- brks[brks >= min(grid) & brks <= max(grid)]
    pb <- ggplot2::ggplot(data.frame(bp = brks), ggplot2::aes(bp)) +
      ggplot2::geom_histogram(binwidth = 5, fill = "#7eb3a8", colour = "white") +
      ggplot2::labs(title = paste0(omics, " (", sx, "): breakpoint distribution"),
                    x = "Breakpoint age, years", y = "Features, count") + theme_soda()
    save_plot(pb, file.path(PATHS$pics, sprintf("08_breakpoints_%s_%s", omics, sx)), w = 8, h = 3.8)
  }
  list(peaks = peaks, breakpoints = brks, rate = data.frame(age = grid, rate = rate),
       feat_breaks = if (length(fb)) data.frame(feature = names(fb), breakpoint = unlist(fb)) else NULL)
}

# ---- (C) Mfuzz: мягкая кластеризация траекторий «фича vs возраст» ----------
# Число кластеров выбираем ПО ДАННЫМ (метод локтя на профилях), не фиксируем.
choose_k <- function(prof, omics = "x", sx = "x", kmax = 10) {
  kmax <- min(kmax, nrow(prof) - 1, 10); ks <- 2:kmax
  wss <- sapply(ks, function(k) sum(stats::kmeans(prof, centers = k, nstart = 5, iter.max = 50)$withinss))
  # «колено» = точка кривой (k, wss), максимально удалённая от прямой (первая->последняя)
  x <- (ks - min(ks)) / (max(ks) - min(ks)); y <- (wss - min(wss)) / (max(wss) - min(wss))
  d <- abs(y - (1 - x))                                  # расстояние до диагонали (норм.)
  k <- ks[which.max(d)]
  pe <- ggplot2::ggplot(data.frame(k = ks, wss = wss), ggplot2::aes(k, wss)) +
    ggplot2::geom_line(colour = "#7799c4") + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = k, linetype = "dashed", colour = "#b5495b") +
    ggplot2::scale_x_continuous(breaks = ks) +
    ggplot2::labs(title = paste0(omics, " (", sx, "): cluster-number elbow (k=", k, ")"),
                  x = "Number of clusters k", y = "Within-cluster sum of squares") + theme_soda()
  save_plot(pe, file.path(PATHS$pics, sprintf("08_elbow_%s_%s", omics, sx)), w = 6, h = 4)
  max(2, k)
}
trend_clusters <- function(mat, age, feats, omics, sx, k = NULL) {
  if (length(feats) < 15) return(NULL)
  if (!requireNamespace("Mfuzz", quietly = TRUE) || !requireNamespace("e1071", quietly = TRUE)) return(NULL)
  suppressPackageStartupMessages({ require(e1071); require(Mfuzz) })   # mfuzz зовёт cmeans() из e1071 без префикса
  bins <- cut(age, breaks = AGE_BIN_BREAKS, include.lowest = TRUE)
  prof <- t(apply(mat[feats, , drop = FALSE], 1, function(v) tapply(v, bins, mean, na.rm = TRUE)))
  prof <- prof[, colSums(is.finite(prof)) > 0, drop = FALSE]
  for (j in seq_len(ncol(prof))) prof[!is.finite(prof[,j]), j] <- mean(prof[,j], na.rm = TRUE)
  if (is.null(k)) k <- choose_k(t(scale(t(prof))), omics, sx)
  log_msg("  кластеров выбрано (локоть): ", k)
  es <- Biobase::ExpressionSet(as.matrix(prof)); es <- Mfuzz::standardise(es)
  cl <- Mfuzz::mfuzz(es, c = k, m = Mfuzz::mestimate(es))
  ctr <- cl$centers
  dd <- do.call(rbind, lapply(seq_len(nrow(ctr)), function(i)
    data.frame(cluster = sprintf("cluster %d (n=%d)", i, sum(cl$cluster == i)),
               bin = seq_len(ncol(ctr)), value = ctr[i, ])))
  p <- ggplot2::ggplot(dd, ggplot2::aes(bin, value, group = cluster)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_line(colour = "#3b6ea5", linewidth = 1) + ggplot2::facet_wrap(~ cluster) +
    ggplot2::labs(title = paste0(omics, " (", sx, "): soft-clustered age trajectories (k=", k, ")"),
                  x = "Age bin, decade", y = "Standardised level") + theme_soda()
  save_plot(p, file.path(PATHS$pics, sprintf("08_trend_clusters_%s_%s", omics, sx)), w = 9, h = 6)
  data.frame(feature = rownames(prof), cluster = cl$cluster)
}

# ---- CpG -> ген (UCSC RefGene из 450k-аннотации) ---------------------------
cpg_to_gene <- function(cpgs) {
  an <- tryCatch(.anno450k(), error = function(e) NULL); if (is.null(an)) return(character(0))
  g <- unlist(strsplit(as.character(an[cpgs, "UCSC_RefGene_Name"]), ";")); unique(g[nzchar(g)])
}

# ---- Валидация трендов по GenAge (de Magalhães signature of ageing) --------
# Скачай список со страницы GenAge (over-/under-expressed с возрастом), сохрани
# CSV в data/raw/genage_signature.csv. Парсер гибкий: ищет колонку символа гена
# и колонку направления (over/under | up/down | +/-).
genage_load <- function(path = file.path(PATHS$raw_expr, "..", "genage_signature.csv")) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(data.table::fread(path), error = function(e) NULL); if (is.null(d)) return(NULL)
  nm <- tolower(names(d))
  sc <- which(nm %in% c("symbol","gene_symbol","gene","name","genesymbol"))[1]
  dc <- grep("^dir$|regulat|direction|express|over|under|trend|change", nm)[1]
  if (is.na(sc) || is.na(dc)) return(NULL)
  sym <- toupper(trimws(as.character(d[[sc]]))); val <- tolower(as.character(d[[dc]]))
  dir <- ifelse(grepl("over|up|increase|\\+", val), "up",
         ifelse(grepl("under|down|decrease|-", val), "down", NA))
  data.frame(symbol = sym, dir = dir)[!is.na(dir) & nzchar(sym), ]
}
genage_concord <- function(tt, genage, omics, sx, genes = NULL) {
  if (is.null(genage)) { log_msg("  GenAge: файл не найден — валидация пропущена (положи data/raw/genage_signature.csv)"); return(NULL) }
  feat_dir <- if (omics == "expression") setNames(tt$dir, tt$feature) else NULL
  sig <- tt[tt$adj.P.Val < FDR, ]
  if (omics == "methylation" && !is.null(genes)) {
    # gene-уровень для метилирования: знак = большинство значимых CpG гена
    return(NULL)  # для метилирования валидация направления делается в интеграции (метил↓~экспр↑)
  }
  m <- merge(data.frame(symbol = sig$feature, our = sig$dir), genage, by = "symbol")
  if (!nrow(m)) { log_msg("  GenAge: пересечений с сигнатурой нет"); return(NULL) }
  conc <- mean(m$our == m$dir); n <- nrow(m)
  pb <- stats::binom.test(sum(m$our == m$dir), n, 0.5, alternative = "greater")$p.value
  data.table::fwrite(m, file.path(PATHS$tables, sprintf("08_genage_validation_%s_%s.csv", omics, sx)))
  log_msg(sprintf("  GenAge %s/%s: совпадение направления %d/%d = %.0f%% (binom p=%.3g)",
                  omics, sx, sum(m$our == m$dir), n, 100*conc, pb))
  list(n = n, concordance = conc, p = pb, table = m)
}

# ---- обогащение Reactome (пути R-HSA, как в SODA) и GO BP ------------------
enrich_reactome <- function(genes) {
  if (!requireNamespace("ReactomePA", quietly = TRUE) || !requireNamespace("clusterProfiler", quietly = TRUE)) return(NULL)
  eid <- tryCatch(clusterProfiler::bitr(unique(genes), "SYMBOL", "ENTREZID", org.Hs.eg.db::org.Hs.eg.db)$ENTREZID,
                  error = function(e) NULL)
  if (is.null(eid) || length(eid) < 5) return(NULL)
  r <- tryCatch(ReactomePA::enrichPathway(eid, pAdjustMethod = "BH", pvalueCutoff = 0.1, readable = TRUE),
                error = function(e) NULL)
  if (is.null(r) || !nrow(as.data.frame(r))) return(NULL); as.data.frame(r)
}
enrich_go <- function(genes) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(NULL)
  r <- tryCatch(clusterProfiler::enrichGO(unique(genes), org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.1), error = function(e) NULL)
  if (is.null(r) || !nrow(as.data.frame(r))) return(NULL); as.data.frame(r)
}
# ---- карта обогащённых терминов в hallmarks of aging (López-Otín 2023) -----
HALLMARKS <- list(
  "Genomic instability"      = "DNA (repair|damage)|double.strand|mismatch|genom|mutation|recombination",
  "Telomere attrition"       = "telomer|shelterin",
  "Epigenetic alterations"   = "methylat|histone|chromatin|epigenet|acetylat|HDAC",
  "Loss of proteostasis"     = "proteasom|unfolded protein|chaperone|protein folding|ubiquitin",
  "Disabled macroautophagy"  = "autophag|lysosom|mitophag",
  "Deregulated nutrient sensing" = "insulin|mTOR|IGF|AMPK|FOXO|nutrient|growth factor",
  "Mitochondrial dysfunction"= "mitochond|oxidative phosphoryl|respiratory|ATP synth|electron transport",
  "Cellular senescence"      = "senescen|cell cycle|p53|CDKN|G1.S|RB1|cycle arrest",
  "Stem cell exhaustion"     = "stem cell|hematopoiet|differentiation|self.renewal|progenitor",
  "Altered intercellular comm." = "cytokine|chemokine|interleukin|NF.k|signaling|FLT3",
  "Chronic inflammation"     = "inflammat|interferon|innate immun|complement|immune response")
map_hallmarks <- function(terms) {
  s <- sapply(HALLMARKS, function(rx) sum(grepl(rx, terms, ignore.case = TRUE)))
  data.frame(hallmark = names(s), n_terms = as.integer(s))[order(-s), ]
}
# ---- барплот топ-обогащений (горизонтальные бары -log10 p) -----------------
barplot_enrich <- function(df, title, file_noext, n = 12, idcol = "Description", pcol = "p.adjust") {
  if (is.null(df) || !nrow(df)) return(invisible())
  d <- df[order(df[[pcol]]), ][seq_len(min(n, nrow(df))), ]
  d$term <- factor(d[[idcol]], levels = rev(d[[idcol]]))
  d$logp <- -log10(d[[pcol]])
  p <- ggplot2::ggplot(d, ggplot2::aes(logp, term)) +
    ggplot2::geom_col(fill = "#b5495b") +
    ggplot2::labs(title = title, x = "-log10 adjusted p", y = NULL) + theme_soda() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
  save_plot(p, file_noext, w = 9, h = 5)
}

# ---- комбинированные графики обогащения (bar + bubble, up/down вместе) ------
.wrap_terms <- function(x, n = 42) vapply(as.character(x), function(s) paste(strwrap(s, n), collapse = "\n"), "")
# genes_a/genes_b — два набора генов (up/down или hyper/hypo или silenced/activated);
# db = "reactome" | "go". Пишет CSV (в tables) + _bar.* + _bubble.* (в pics).
enrich_plots <- function(genes_a, genes_b, lab_a, lab_b, file_noext, title, db = "reactome", topn = 10) {
  enr <- function(g) if (db == "reactome") enrich_reactome(g) else enrich_go(g)
  fa <- enr(genes_a); fb <- enr(genes_b)
  base <- basename(file_noext)
  if (!is.null(fa)) data.table::fwrite(fa, file.path(PATHS$tables, paste0(base, "_", lab_a, ".csv")))
  if (!is.null(fb)) data.table::fwrite(fb, file.path(PATHS$tables, paste0(base, "_", lab_b, ".csv")))
  mk <- function(d, lb) { if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[order(d$p.adjust), ][seq_len(min(topn, nrow(d))), ]
    gr <- sapply(strsplit(as.character(d$GeneRatio), "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
    data.frame(term = d$Description, padj = d$p.adjust, Count = d$Count, GeneRatio = gr, direction = lb) }
  D <- rbind(mk(fa, lab_a), mk(fb, lab_b)); if (is.null(D) || !nrow(D)) return(invisible())
  D$term <- .wrap_terms(D$term)
  D <- D[order(D$direction, -D$padj), ]; D$term <- factor(D$term, levels = unique(D$term))
  D$direction <- factor(D$direction, levels = c(lab_a, lab_b))
  th <- ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                       strip.text.y = ggplot2::element_text(angle = 0, face = "bold"),
                       plot.title = ggplot2::element_text(hjust = 0))
  pb <- ggplot2::ggplot(D, ggplot2::aes(-log10(padj), term, fill = padj)) + ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "#c0392b", high = "#3b6ea5", name = "p.adjust") +
    ggplot2::facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = title, x = "-log10 adjusted p", y = NULL) + theme_soda() + th
  save_plot(pb, paste0(file_noext, "_bar"), w = 10, h = 7)
  pp <- ggplot2::ggplot(D, ggplot2::aes(GeneRatio, term, size = Count, colour = padj)) + ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#c0392b", high = "#3b6ea5", name = "p.adjust") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Count") +
    ggplot2::facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = title, x = "Gene ratio", y = NULL) + theme_soda() + th
  save_plot(pp, paste0(file_noext, "_bubble"), w = 10, h = 7)
}
