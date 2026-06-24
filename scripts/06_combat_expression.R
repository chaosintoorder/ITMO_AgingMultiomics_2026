# ============================================================================
# 06_combat_expression.R — кросс-студийная гармонизация экспрессии (ComBat).
#   mod = ~ Age + Sex + доли клеток (биология и состав крови СОХРАНЯЮТСЯ),
#   batch = датасет, опорный батч = крупнейший. Метод: Johnson 2007 (emp. Bayes).
# Выход: 06_expr_combat.rds + ridgeline/violin «после» + PCA «до/после».
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(sva) })

ex <- readRDS(file.path(PATHS$final_expr, "03_expr_processed.rds"))   # expr, meta, batch
cp <- readRDS(file.path(PATHS$interim_expr, "05_cellprop_expr.rds"))
X <- ex$expr; batch <- ex$batch; meta <- data.table::as.data.table(ex$meta)
log_msg("вход: ", nrow(X), " генов × ", ncol(X), " образцов | датасетов ", length(unique(batch)))

covdf <- build_covars(colnames(X), meta, cp)
log_msg("ковариаты ComBat: ", paste(names(covdf), collapse = ", "))

ageX <- as.numeric(meta$Age[match(colnames(X), meta$gsm)])
# PCA «до»
pca_batch_plot(X, batch, "Expression PCA before ComBat",
               file.path(PATHS$pics, "06_pca_expr_before"), age = ageX)

Xc <- run_combat(X, batch, covdf)

# PCA «после»
pca_batch_plot(Xc, batch, "Expression PCA after ComBat",
               file.path(PATHS$pics, "06_pca_expr_after"), age = ageX)

# распределения «после» (ridgeline + violin), counts = реальные размеры выборок
set.seed(1)
df <- do.call(rbind, lapply(unique(batch), function(g) {
  v <- as.numeric(Xc[, batch == g]); v <- v[is.finite(v)]
  data.frame(dataset = g, value = sample(v, min(5000, length(v)))) }))
xy <- "Expression level, log2-scaled + ComBat-corrected"
nper <- table(batch)
make_ridgeline(df, EXPR_COL, xlab = xy, title = "Expression by dataset (after ComBat)",
               file = file.path(PATHS$pics, "06_expr_ridgeline_after.png"), counts = nper)
make_violin(df, EXPR_COL, ylab = xy, title = "Expression by dataset (after ComBat)",
            file = file.path(PATHS$pics, "06_expr_violin_after.png"), counts = nper)

meta_out <- meta[match(colnames(Xc), gsm)]
save_rds(list(expr = Xc, meta = meta_out, batch = batch, cellprop = cp[colnames(Xc), , drop = FALSE]),
         file.path(PATHS$final_expr, "06_expr_combat.rds"))
log_msg("06 готов: data/final/expression/06_expr_combat.rds")
