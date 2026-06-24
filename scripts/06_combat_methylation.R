# ============================================================================
# 06_combat_methylation.R — гармонизация метилирования (ComBat на M-значениях).
#   1) грузим слим-датасеты (04_*.rds), снимаем дубли CpG, берём ОБЩИЕ CpG,
#   2) ComBat на M (Du 2010), mod = ~ Age + Sex + доли клеток, опорный = крупнейший,
#   3) для отображения M -> beta. PCA «до/после» — показать снижение батч-эффекта.
# Выход: 06_meth_combat.rds + ridgeline/violin «после» (beta) + PCA «до/после».
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(sva) })

mf  <- list.files(PATHS$interim_meth, pattern = "^04_GSE.*\\.rds$", full.names = TRUE)
obj <- lapply(mf, readRDS)
gse <- sapply(obj, function(x) x$gse)
# снять возможные дубли rownames внутри датасета (GSE87571: две матрицы)
Ms  <- lapply(obj, function(x) { M <- x$M; M[!duplicated(rownames(M)), , drop = FALSE] })
log_msg("слим-датасеты: ", paste0(gse, "(", sapply(Ms, ncol), ")", collapse = ", "))

common <- Reduce(intersect, lapply(Ms, rownames))
log_msg("ОБЩИХ CpG по ", length(Ms), " датасетам: ", length(common))
if (length(common) < 1000)
  log_msg("  ВНИМАНИЕ: общих CpG мало (<1000) — напиши мне число, подкрутим отбор.")

M <- do.call(cbind, lapply(Ms, function(x) x[common, , drop = FALSE]))
batch <- rep(gse, sapply(Ms, ncol))
log_msg("объединённая M: ", nrow(M), " CpG × ", ncol(M), " образцов")

meta <- data.table::as.data.table(readRDS(file.path(PATHS$interim_meth, "02_meta_meth.rds")))
cp   <- readRDS(file.path(PATHS$interim_meth, "05_cellprop_meth.rds"))
covdf <- build_covars(colnames(M), meta, cp)
log_msg("ковариаты ComBat: ", paste(names(covdf), collapse = ", "))

ageM <- as.numeric(meta$Age[match(colnames(M), meta$gsm)])
pca_batch_plot(M, batch, "Methylation PCA before ComBat (M-values)",
               file.path(PATHS$pics, "06_pca_meth_before"), age = ageM)
Mc <- run_combat(M, batch, covdf)
pca_batch_plot(Mc, batch, "Methylation PCA after ComBat (M-values)",
               file.path(PATHS$pics, "06_pca_meth_after"), age = ageM)

set.seed(1)
df <- do.call(rbind, lapply(unique(batch), function(g) {
  v <- m_to_beta(as.numeric(Mc[, batch == g])); v <- v[is.finite(v)]
  data.frame(dataset = g, value = sample(v, min(5000, length(v)))) }))
xy <- "Methylation level, \u03b2-value (ComBat-corrected)"
nper <- table(batch)
make_ridgeline(df, METH_COL, xlab = xy, title = "Methylation by dataset (after ComBat)",
               file = file.path(PATHS$pics, "06_meth_ridgeline_after.png"), counts = nper)
make_violin(df, METH_COL, ylab = xy, title = "Methylation by dataset (after ComBat)",
            file = file.path(PATHS$pics, "06_meth_violin_after.png"), counts = nper)

meta_out <- meta[match(colnames(Mc), gsm)]
save_rds(list(M = Mc, meta = meta_out, batch = batch, cellprop = cp[colnames(Mc), , drop = FALSE]),
         file.path(PATHS$final_meth, "06_meth_combat.rds"))
log_msg("06 (метилирование) готов: data/final/methylation/06_meth_combat.rds | общих CpG ", length(common))
