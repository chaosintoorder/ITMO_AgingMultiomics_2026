# ============================================================================
# 08_integration.R — ЕДИНСТВЕННАЯ точка слияния омик + БИОЛОГИЯ.
#   (1) пересечение возраст-ассоц. генов (гипергеометрия) + согласованность знаков
#       (канон: экспр↓ & гиперметил = «выключение»; экспр↑ & гипометил = «активация»)
#   (2) кандидаты-драйверы = в обеих омиках + согласованы (+ флаг GenAge)
#   (3) ОБОГАЩЕНИЕ кандидатов: Reactome (R-HSA) + GO BP -> биологические выводы
#   (4) карта обогащений в HALLMARKS OF AGING (López-Otín 2023) + барплот
#   (5) «transition genes» (гены с точкой перехода) -> для STRING-сети (09)
#   (6) overlay критических возрастов экспрессии и метилирования
# ============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })

E <- readRDS(file.path(PATHS$final_expr, "07_trends_expression.rds"))
M <- readRDS(file.path(PATHS$final_meth, "07_trends_methylation.rds"))
genage <- genage_load(); ga_sym <- if (!is.null(genage)) unique(genage$symbol) else character(0)
UNIV <- 20000

hall_all <- list()
for (sx in c("F","M")) {
  e <- E[[sx]]; m <- M[[sx]]; if (is.null(e) || is.null(m)) next
  eg <- unique(c(e$genes_up, e$genes_down)); mg <- unique(c(m$genes_up, m$genes_down))
  ov <- intersect(eg, mg)
  p_hyper <- stats::phyper(length(ov)-1, length(eg), UNIV-length(eg), length(mg), lower.tail = FALSE)
  silenced  <- intersect(e$genes_down, m$genes_up)     # экспр↓ + гиперметил
  activated <- intersect(e$genes_up,   m$genes_down)   # экспр↑ + гипометил
  concord <- union(silenced, activated)

  cand <- data.frame(gene = concord,
                     pattern = ifelse(concord %in% silenced, "silenced (expr down / hyper-meth)",
                                                             "activated (expr up / hypo-meth)"),
                     in_GenAge = concord %in% ga_sym)
  cand <- cand[order(!cand$in_GenAge, cand$gene), ]
  data.table::fwrite(data.frame(gene = ov), file.path(PATHS$tables, sprintf("08_integration_overlap_%s.csv", sx)))
  data.table::fwrite(cand, file.path(PATHS$tables, sprintf("08_candidate_drivers_%s.csv", sx)))
  log_msg(sprintf("ИНТЕГРАЦИЯ %s: экспр-генов=%d, метил-генов=%d, пересечение=%d (p_hyper=%.2e)",
                  sx, length(eg), length(mg), length(ov), p_hyper))
  log_msg(sprintf("  согласованных: выключено=%d, активировано=%d | в GenAge=%d",
                  length(silenced), length(activated), sum(cand$in_GenAge)))

  # (3) обогащение кандидатов-драйверов: Reactome + GO -> биология
  genes_cand <- unique(c(concord, ov))                 # согласованные + всё пересечение
  re <- enrich_reactome(genes_cand); go <- enrich_go(genes_cand)   # для карты hallmarks (термины)
  # фигуры в стиле SODA: silenced vs activated на одной картинке (bar + bubble)
  tryCatch({
    enrich_plots(silenced, activated, "silenced", "activated",
                 file.path(PATHS$pics, sprintf("08_candidates_Reactome_%s", sx)),
                 sprintf("Integrated candidates %s: Reactome", sx), db = "reactome")
    enrich_plots(silenced, activated, "silenced", "activated",
                 file.path(PATHS$pics, sprintf("08_candidates_GO_%s", sx)),
                 sprintf("Integrated candidates %s: GO BP", sx), db = "go")
  }, error = function(e) log_msg("  обогащение кандидатов/", sx, " пропущено: ", conditionMessage(e)))

  # (4) карта в hallmarks of aging
  terms <- c(if (!is.null(re)) re$Description, if (!is.null(go)) go$Description)
  if (length(terms)) {
    hm <- map_hallmarks(terms); hm$sex <- sx; hall_all[[sx]] <- hm
    data.table::fwrite(hm, file.path(PATHS$tables, sprintf("08_hallmarks_%s.csv", sx)))
  }

  # (5) «transition genes» (гены экспрессии с точкой перехода) -> для STRING (09)
  if (!is.null(e$feat_breaks)) {
    tg <- e$feat_breaks; tg$sex <- sx
    data.table::fwrite(tg, file.path(PATHS$tables, sprintf("08_transition_genes_%s.csv", sx)))
  }

  # (6) overlay критических возрастов
  if (!is.null(e$rate) && !is.null(m$rate)) {
    ddr <- rbind(data.frame(e$rate, omics = "expression"), data.frame(m$rate, omics = "methylation"))
    p <- ggplot(ddr, aes(age, rate, colour = omics)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c(expression = "#3b6ea5", methylation = "#d8694f"), name = NULL) +
      labs(title = paste0("Critical ages overlay (", sx, ")"),
           x = "Age, years", y = "Rate of change, scaled") + theme_soda()
    save_plot(p, file.path(PATHS$pics, sprintf("08_critical_ages_overlay_%s", sx)), w = 8, h = 4.2)
  }
}

# сводный барплот hallmarks (оба пола)
if (length(hall_all)) {
  H <- do.call(rbind, hall_all); H <- H[H$n_terms > 0, ]
  p <- ggplot(H, aes(reorder(hallmark, n_terms), n_terms, fill = sex)) +
    geom_col(position = position_dodge()) + coord_flip() +
    scale_fill_manual(values = c(F = "#b5495b", M = "#3b6ea5"), name = "Sex") +
    labs(title = "Aging hallmarks enriched in integrated candidates",
         x = NULL, y = "Enriched terms mapped, count") + theme_soda()
  save_plot(p, file.path(PATHS$pics, "08_hallmarks_summary"), w = 9, h = 5)
}
log_msg("08 (интеграция) готов.")
