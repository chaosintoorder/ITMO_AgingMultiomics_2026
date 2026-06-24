# ============================================================================
# 05_deconv_methylation.R — состав крови для МЕТИЛИРОВАНИЯ (EpiDISH, посчитан в 04).
# 7 типов (B/CD4T/CD8T/NK/Mono/Neutro/Eosino). Здесь только собираем + рисуем.
# Выход: 05_cellprop_meth.rds + хитмэп по датасетам + «доля vs возраст» (CI + p GAM).
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(mgcv) })
HEAT_PAL <- c("#c60000","#df8668","#eac493","#eee8ae","#b7c0ad","#8b9fb0","#3861ae")

mf  <- list.files(PATHS$interim_meth, pattern = "^04_GSE.*\\.rds$", full.names = TRUE)
obj <- lapply(mf, readRDS)
cp  <- do.call(rbind, lapply(obj, function(x) if (!is.null(x$cellprop)) x$cellprop else NULL))
if (is.null(cp) || !nrow(cp)) stop("нет долей клеток в 04_*.rds — перезапусти 04 после обновления utils_io.")
save_rds(cp, file.path(PATHS$interim_meth, "05_cellprop_meth.rds"))
log_msg("метилирование: доли EpiDISH ", nrow(cp), " × ", ncol(cp), " (датасетов: ", length(obj), ")")

meta <- data.table::as.data.table(readRDS(file.path(PATHS$interim_meth, "02_meta_meth.rds")))
# --- хитмэп клеточного состава по датасетам ---
cellcomp_heatmap(cp, meta, HEAT_PAL, "Methylation: blood cell fractions per dataset (EpiDISH RPC)",
                 "Estimated\nfraction", file.path(PATHS$pics, "05_cellcomp_heatmap_methylation"))

# --- доля vs возраст: GAM + 95% CI, значимость тренда s(age) ---
age <- as.numeric(meta$Age[match(rownames(cp), meta$gsm)])
pv <- sapply(colnames(cp), function(k) {
  g <- tryCatch(mgcv::gam(cp[, k] ~ s(age)), error = function(e) NULL)
  if (is.null(g)) NA else summary(g)$s.pv[1] })
lab <- setNames(sprintf("%s  (%s)", colnames(cp), fmt_pval(pv)), colnames(cp))
df <- do.call(rbind, lapply(colnames(cp), function(k)
  data.frame(cell = lab[k], age = age, frac = cp[, k], row.names = NULL)))
p <- ggplot(df, aes(age, frac)) + geom_point(alpha = .12, size = .5) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "#d8694f", fill = "#d8694f", alpha = .2) +
  facet_wrap(~ cell, scales = "free_y") +
  labs(title = "Methylation: cell fractions vs age (GAM fit, 95% CI)",
       x = "Age, years", y = "Estimated fraction") + theme_soda()
save_plot(p, file.path(PATHS$pics, "05_cellfraction_vs_age_methylation"), w = 10, h = 6.5)
log_msg("05 (метилирование) готов.")
