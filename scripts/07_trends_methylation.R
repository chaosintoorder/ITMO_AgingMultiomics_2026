# ============================================================================
# 07_trends_methylation.R — тренды МЕТИЛИРОВАНИЯ, форк по полу.
#   (A) limma ~ возраст + доли клеток на M-значениях, BH-FDR -> возраст-ассоц. CpG
#       (знак: up = гиперметилирование с возрастом, down = гипо)
#   (B) GAM-производная + точки разрыва -> критические возрасты
#   (C) Mfuzz -> типы трендов
#   (D) обогащение: missMethyl::gometh (КОРРЕКЦИЯ на число CpG на ген — важно!)
#   (E) CpG -> гены (для интеграции с экспрессией)
# Выход: таблицы 07_*, графики 08_*, объект 07_trends_methylation.rds.
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(limma); library(mgcv) })
TOP_GAM <- 200

me <- readRDS(file.path(PATHS$final_meth, "06_meth_combat.rds"))
meta <- data.table::as.data.table(me$meta)
age_all <- as.numeric(meta$Age[match(colnames(me$M), meta$gsm)])
sex_all <- meta$Sex[match(colnames(me$M), meta$gsm)]

res <- list()
for (sx in c("F","M")) {
  idx <- which(sex_all == sx & is.finite(age_all))
  if (length(idx) < MIN_SAMPLES_WIN) { log_msg(sx, ": мало образцов (", length(idx), ") — пропуск"); next }
  m <- me$M[, idx, drop = FALSE]; ag <- age_all[idx]
  log_msg("=== METHYLATION / пол ", sx, " | образцов ", length(idx), " | CpG ", nrow(m), " ===")

  tt <- age_limma(m, ag, me$cellprop)
  data.table::fwrite(tt, file.path(PATHS$tables, sprintf("07_age_assoc_methylation_%s.csv", sx)))
  sig <- tt$feature[tt$adj.P.Val < FDR]
  up  <- tt$feature[tt$dir=="up"   & tt$adj.P.Val < FDR]   # гиперметилирование с возрастом
  dn  <- tt$feature[tt$dir=="down" & tt$adj.P.Val < FDR]   # гипометилирование с возрастом
  log_msg("  значимых CpG (FDR<", FDR, "): ", length(sig), " | hyper=", length(up), " hypo=", length(dn))
  if (!length(sig)) next
  feats <- head(tt$feature[tt$adj.P.Val < FDR], TOP_GAM)

  ca <- critical_ages(m, ag, feats, "methylation", sx)
  log_msg("  критические возрасты (пики скорости): ", paste(round(ca$peaks), collapse = ", "))
  trend_clusters(m, ag, feats, "methylation", sx)

  gup <- cpg_to_gene(up); gdn <- cpg_to_gene(dn)          # гены из гипер/гипо-CpG (для интеграции и путей)

  # (D) gometh — обогащение GO с поправкой на число CpG на ген (нужен пакет аннотации 450k ПОДКЛЮЧЁННЫМ)
  tryCatch({
    if (requireNamespace("missMethyl", quietly = TRUE)) {
      suppressPackageStartupMessages(require(IlluminaHumanMethylation450kanno.ilmn12.hg19))
      gm <- missMethyl::gometh(sig.cpg = sig, all.cpg = tt$feature, collection = "GO", array.type = "450K")
      gm <- gm[order(gm$P.DE), ]
      data.table::fwrite(gm[gm$FDR < 0.1, ], file.path(PATHS$tables, sprintf("07_GOmeth_methylation_%s.csv", sx)))
    }
  }, error = function(e) log_msg("  gometh METHYLATION/", sx, " пропущено: ", conditionMessage(e)))

  # (D2) графики путей по генам гипер/гипо-CpG (Reactome + GO, bar + bubble) — как у экспрессии
  tryCatch({
    enrich_plots(gup, gdn, "hyper", "hypo", file.path(PATHS$pics, sprintf("07_Reactome_methylation_%s", sx)),
                 sprintf("Methylation %s: Reactome (hyper vs hypo with age)", sx), db = "reactome")
    enrich_plots(gup, gdn, "hyper", "hypo", file.path(PATHS$pics, sprintf("07_GO_methylation_%s", sx)),
                 sprintf("Methylation %s: GO BP (hyper vs hypo with age)", sx), db = "go")
  }, error = function(e) log_msg("  пути METHYLATION/", sx, " пропущено: ", conditionMessage(e)))

  res[[sx]] <- list(tt = tt, sig_cpg = sig, genes_up = gup, genes_down = gdn,
                    peaks = ca$peaks, rate = ca$rate, feat_breaks = ca$feat_breaks)
}
save_rds(res, file.path(PATHS$final_meth, "07_trends_methylation.rds"))
log_msg("07 (метилирование) готов.")
