# ============================================================================
# 09_networks.R — ГЕННЫЕ СЕТИ (НЕОБЯЗАТЕЛЬНЫЙ, тяжелее остального).
#   (A) STRING PPI на «transition genes» (гены с точкой перехода) — это прямой
#       аналог SODA-фигур «genes with transitions»: сеть + хабы по связности +
#       картинка сети. Считаем по полу.
#   (B) STRING на кандидатах-драйверах (кросс-омиксные согласованные гены).
#   (C) WGCNA: модули коэкспрессии -> корреляция с возрастом -> хабы (kME).
#   Сети STRING качают данные из интернета — поднимаем таймаут.
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })
options(timeout = 1200)                                   # STRING тянет ~20 МБ, дефолтные 60с мало

# ---------- читаемая сеть через igraph: подписи топ-хабов, png+svg ----------
draw_network <- function(ints, mp, tag, label_top = 35) {
  if (!requireNamespace("igraph", quietly = TRUE) || is.null(ints) || !nrow(ints)) return(invisible())
  id2g <- setNames(mp$gene, mp$STRING_id)
  e <- data.frame(from = id2g[ints$from], to = id2g[ints$to])
  e <- e[!is.na(e$from) & !is.na(e$to) & e$from != e$to, , drop = FALSE]; if (!nrow(e)) return(invisible())
  g <- igraph::simplify(igraph::graph_from_data_frame(e, directed = FALSE))
  deg <- igraph::degree(g)
  lab <- ifelse(rank(-deg, ties.method = "first") <= label_top, igraph::V(g)$name, NA)  # подписи только хабам
  set.seed(1); lay <- igraph::layout_with_fr(g)
  draw <- function() { op <- graphics::par(mar = c(0,0,2,0)); on.exit(graphics::par(op))
    plot(g, layout = lay, vertex.size = 2 + 7 * (deg / max(deg)), vertex.color = "#7799c4",
         vertex.frame.color = "white", vertex.label = lab, vertex.label.cex = 0.85,
         vertex.label.color = "black", vertex.label.family = "sans", vertex.label.dist = 0.25,
         edge.color = "grey80", edge.width = 0.4, main = tag) }
  png(file.path(PATHS$pics, sprintf("09_string_network_%s.png", tag)), width = 2800, height = 2200, res = 220)
  draw(); grDevices::dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(PATHS$pics, sprintf("09_string_network_%s.svg", tag)), width = 13, height = 10)
    draw(); grDevices::dev.off() }
  log_msg("  сеть: 09_string_network_", tag, " (png+svg) | узлов ", igraph::vcount(g), " рёбер ", igraph::ecount(g))
}

# ---------- STRING-сеть по набору генов: хабы + картинка ----------
string_network <- function(genes, tag) {
  if (!requireNamespace("STRINGdb", quietly = TRUE) || length(genes) < 5) {
    log_msg("  STRING(", tag, ") пропущен (нет пакета/мало генов)"); return(invisible()) }
  sdb <- tryCatch(STRINGdb::STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400),
                  error = function(e) { log_msg("  STRING init: ", conditionMessage(e)); NULL })
  if (is.null(sdb)) return(invisible())
  mp <- tryCatch(sdb$map(data.frame(gene = unique(genes)), "gene", removeUnmappedRows = TRUE),
                 error = function(e) { log_msg("  STRING map(", tag, "): ", conditionMessage(e)); NULL })
  if (is.null(mp) || !nrow(mp)) return(invisible())
  ints <- sdb$get_interactions(mp$STRING_id)
  if (nrow(ints)) {
    deg <- sort(table(c(ints$from, ints$to)), decreasing = TRUE)
    hub <- merge(data.frame(STRING_id = names(deg), degree = as.integer(deg)), mp, by = "STRING_id")
    hub <- hub[order(-hub$degree), c("gene","degree")]
    data.table::fwrite(hub, file.path(PATHS$tables, sprintf("09_string_hubs_%s.csv", tag)))
    log_msg("  STRING(", tag, "): ", nrow(mp), " генов, ", nrow(ints), " связей | хабы: ",
            paste(head(hub$gene, 10), collapse = ", "))
  }
  draw_network(ints, mp, tag)                            # читаемая сеть igraph (png+svg)
}

# (A) transition genes (с точкой перехода) — по полу
for (sx in c("F","M")) {
  f <- file.path(PATHS$tables, sprintf("08_transition_genes_%s.csv", sx))
  if (file.exists(f)) {
    tg <- data.table::fread(f)
    string_network(tg$feature, paste0("transition_", sx))
  }
}
# (B) кандидаты-драйверы (оба пола вместе)
cf <- list.files(PATHS$tables, pattern = "^08_candidate_drivers_.*\\.csv$", full.names = TRUE)
cand <- unique(unlist(lapply(cf, function(f) data.table::fread(f)$gene)))
string_network(cand, "candidates")

# (C) WGCNA — модули коэкспрессии и хабы (через library, иначе конфликт cor())
tryCatch({
  if (requireNamespace("WGCNA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(WGCNA)); cor <- WGCNA::cor   # критично: иначе "unused arguments"
    ex <- readRDS(file.path(PATHS$final_expr, "06_expr_combat.rds"))
    E  <- readRDS(file.path(PATHS$final_expr, "07_trends_expression.rds"))
    meta <- data.table::as.data.table(ex$meta)
    sig <- intersect(unique(unlist(lapply(E, function(x) x$sig))), rownames(ex$expr))
    log_msg("WGCNA на ", length(sig), " возраст-ассоц. генах")
    datExpr <- t(ex$expr[sig, , drop = FALSE]); age <- as.numeric(meta$Age[match(rownames(datExpr), meta$gsm)])
    sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)
    pw  <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 8
    net <- WGCNA::blockwiseModules(datExpr, power = pw, TOMType = "signed", minModuleSize = 20,
                                   mergeCutHeight = 0.25, numericLabels = TRUE, maxBlockSize = length(sig), verbose = 0)
    ME <- WGCNA::moduleEigengenes(datExpr, net$colors)$eigengenes
    corA <- as.numeric(cor(ME, age, use = "pairwise.complete.obs"))
    modtab <- data.frame(module = colnames(ME), cor_age = round(corA, 3))[order(-abs(corA)), ]
    data.table::fwrite(modtab, file.path(PATHS$tables, "09_wgcna_modules_age.csv"))
    log_msg("  WGCNA: модулей ", nrow(modtab), " | самый возраст-связанный r=", modtab$cor_age[1])
    mt <- modtab; mt$module <- sub("^ME", "M", mt$module)
    pe <- ggplot(mt, aes(reorder(module, cor_age), cor_age, fill = cor_age > 0)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c(`TRUE` = "#b5495b", `FALSE` = "#3b6ea5"),
                        labels = c(`TRUE` = "up with age", `FALSE` = "down with age"), name = NULL) +
      labs(title = "WGCNA modules: eigengene–age correlation",
           x = "Module", y = "Correlation with age, r") + theme_soda()
    save_plot(pe, file.path(PATHS$pics, "09_wgcna_modules_age"), w = 7, h = 4.5)
    topmod <- sub("ME","", modtab$module[1]); inmod <- names(net$colors)[net$colors == as.integer(topmod)]
    kME <- as.numeric(cor(datExpr[, inmod], ME[, paste0("ME", topmod)], use = "pairwise.complete.obs"))
    hubs <- data.frame(gene = inmod, kME = round(kME, 3)); hubs <- hubs[order(-hubs$kME), ]
    data.table::fwrite(hubs, file.path(PATHS$tables, "09_wgcna_hub_genes.csv"))
    log_msg("  хабы возраст-модуля (топ): ", paste(head(hubs$gene, 10), collapse = ", "))
    cor <- stats::cor
  }
}, error = function(e) log_msg("  WGCNA пропущен: ", conditionMessage(e)))

log_msg("09 (сети) готов.")
