# ============================================================================
# 05_deconv_expression.R — состав крови для ЭКСПРЕССИИ (маркерные скоры).
# 7 типов (B/CD4T/CD8T/NK/Mono/Neutro/Eosino) — те же, что в метилировании.
# Выход: 05_cellprop_expr.rds + хитмэп по датасетам + «скор vs возраст» (CI + p GAM).
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(mgcv) })
HEAT_PAL <- c("#c60000","#df8668","#eac493","#eee8ae","#b7c0ad","#8b9fb0","#3861ae")

ex <- readRDS(file.path(PATHS$final_expr, "03_expr_processed.rds"))
cp <- deconv_blood_expr(ex$expr)                       # образцы × типы клеток
save_rds(cp, file.path(PATHS$interim_expr, "05_cellprop_expr.rds"))
log_msg("экспрессия: маркерные скоры ", nrow(cp), " × ", ncol(cp), " (", paste(colnames(cp), collapse=", "), ")")

meta <- data.table::as.data.table(ex$meta)
# --- хитмэп клеточного состава по датасетам (образцы упорядочены по возрасту) ---
cellcomp_heatmap(cp, meta, HEAT_PAL, "Expression: marker-based cell scores per dataset",
                 "Marker score\n(z-scored)", file.path(PATHS$pics, "05_cellcomp_heatmap_expression"))

# --- скор vs возраст: GAM + 95% CI, значимость тренда s(age) на фасете ---
age <- as.numeric(meta$Age[match(rownames(cp), meta$gsm)])
pv <- sapply(colnames(cp), function(k) {
  g <- tryCatch(mgcv::gam(cp[, k] ~ s(age)), error = function(e) NULL)
  if (is.null(g)) NA else summary(g)$s.pv[1] })
lab <- setNames(sprintf("%s  (%s)", colnames(cp), fmt_pval(pv)), colnames(cp))
df <- do.call(rbind, lapply(colnames(cp), function(k)
  data.frame(cell = lab[k], age = age, score = cp[, k], row.names = NULL)))
p <- ggplot(df, aes(age, score)) + geom_point(alpha = .12, size = .5) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "#3b6ea5", fill = "#3b6ea5", alpha = .2) +
  facet_wrap(~ cell, scales = "free_y") +
  labs(title = "Expression: cell scores vs age (GAM fit, 95% CI)",
       x = "Age, years", y = "Marker score, z-scored") + theme_soda()
save_plot(p, file.path(PATHS$pics, "05_cellscore_vs_age_expression"), w = 10, h = 6.5)
log_msg("05 (экспрессия) готов.")
