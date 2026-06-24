# ============================================================================
# R/install_packages.R — запусти ОДИН раз перед всем пайплайном.
#   source("R/install_packages.R")
# ============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c("data.table","R.utils","ggplot2","ggridges","mgcv","segmented",
          "matrixStats","future.apply","e1071","RColorBrewer","UpSetR","readxl")
new <- setdiff(cran, rownames(installed.packages()))
if (length(new)) install.packages(new)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c(
  "GEOquery","Biobase","limma","sva","impute",                       # ввод + норм + ComBat + импутация
  "AnnotationDbi","org.Hs.eg.db",                                    # аннотация (база)
  "hgu133plus2.db","hugene10sttranscriptcluster.db",                 # аннотации платформ экспрессии
  "hugene11sttranscriptcluster.db","huex10sttranscriptcluster.db","illuminaHumanv4.db",
  "minfi","IlluminaHumanMethylation450kanno.ilmn12.hg19","EpiDISH",  # метилирование + деконволюция
  "granulator",                                                      # деконволюция экспрессии
  "missMethyl","clusterProfiler","ReactomePA",                       # обогащение
  "WGCNA","STRINGdb","Mfuzz"                                         # сети + кластеризация трендов
)
new_b <- setdiff(bioc, rownames(installed.packages()))
if (length(new_b)) BiocManager::install(new_b, update = FALSE, ask = FALSE)

# ОПЦИОНАЛЬНО (удаление кросс-реактивных зондов 450k по Chen 2013):
#   if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
#   remotes::install_github("markgene/maxprobes")
# Без него шаг просто пропускается (фильтруем только половые хромосомы и NA).

cat("\nГотово. Проверь, что нет ошибок установки выше.\n")
