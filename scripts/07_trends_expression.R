# ============================================================================
# 07_trends_expression.R — тренды ЭКСПРЕССИИ, форк по полу.
#   (A) limma ~ возраст + доли клеток, BH-FDR -> возраст-ассоц. гены (знак ±)
#   (B) GAM-производная + точки разрыва -> критические возрасты
#   (C) Mfuzz -> типы трендов
#   (D) обогащение: enrichGO (BP) + enrichKEGG для up/down
#   (E) валидация GenAge: совпадает ли наш знак с over-/under-expressed с возрастом
# Выход: таблицы 07_*, графики 08_*, объект 07_trends_expression.rds (для интеграции).
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(limma); library(mgcv) })
TOP_GAM <- 200

ex <- readRDS(file.path(PATHS$final_expr, "06_expr_combat.rds"))
meta <- data.table::as.data.table(ex$meta)
age_all <- as.numeric(meta$Age[match(colnames(ex$expr), meta$gsm)])
sex_all <- meta$Sex[match(colnames(ex$expr), meta$gsm)]
genage <- genage_load()

res <- list()
for (sx in c("F","M")) {
  idx <- which(sex_all == sx & is.finite(age_all))
  if (length(idx) < MIN_SAMPLES_WIN) { log_msg(sx, ": мало образцов (", length(idx), ") — пропуск"); next }
  m <- ex$expr[, idx, drop = FALSE]; ag <- age_all[idx]
  log_msg("=== EXPRESSION / пол ", sx, " | образцов ", length(idx), " | генов ", nrow(m), " ===")

  tt <- age_limma(m, ag, ex$cellprop)
  data.table::fwrite(tt, file.path(PATHS$tables, sprintf("07_age_assoc_expression_%s.csv", sx)))
  sig <- tt$feature[tt$adj.P.Val < FDR]
  up  <- tt$feature[tt$dir=="up"   & tt$adj.P.Val < FDR]
  dn  <- tt$feature[tt$dir=="down" & tt$adj.P.Val < FDR]
  log_msg("  значимых (FDR<", FDR, "): ", length(sig), " | up=", length(up), " down=", length(dn))
  if (!length(sig)) next
  feats <- head(tt$feature[tt$adj.P.Val < FDR], TOP_GAM)

  ca <- critical_ages(m, ag, feats, "expression", sx)
  log_msg("  критические возрасты (пики скорости): ", paste(round(ca$peaks), collapse = ", "))
  trend_clusters(m, ag, feats, "expression", sx)
  gv <- genage_concord(tt, genage, "expression", sx)

  # (D) обогащение: GO BP + Reactome (R-HSA, как в SODA), up/down на одной фигуре (bar + bubble)
  tryCatch({
    enrich_plots(up, dn, "up", "down", file.path(PATHS$pics, sprintf("07_Reactome_expression_%s", sx)),
                 sprintf("Expression %s: Reactome (up vs down with age)", sx), db = "reactome")
    enrich_plots(up, dn, "up", "down", file.path(PATHS$pics, sprintf("07_GO_expression_%s", sx)),
                 sprintf("Expression %s: GO BP (up vs down with age)", sx), db = "go")
  }, error = function(e) log_msg("  обогащение EXPRESSION/", sx, " пропущено: ", conditionMessage(e)))

  res[[sx]] <- list(tt = tt, sig = sig, genes_up = up, genes_down = dn,
                    peaks = ca$peaks, rate = ca$rate, feat_breaks = ca$feat_breaks, genage = gv)
}
save_rds(res, file.path(PATHS$final_expr, "07_trends_expression.rds"))
log_msg("07 (экспрессия) готов.")
