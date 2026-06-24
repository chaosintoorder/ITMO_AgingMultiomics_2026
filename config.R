# ============================================================================
# config.R — единый источник правды. Источается в начале КАЖДОГО скрипта.
# Запускай R из КОРНЯ проекта ITMO_AgingMultiomics_2026 (пути относительные).
# ============================================================================
set.seed(1234)
options(stringsAsFactors = FALSE)

# ---- 1. Пути (относительные, без абсолютных и без geo_cache) ---------------
PATHS <- list(
  raw_expr     = "data/raw/expression",
  raw_meth     = "data/raw/methylation",
  interim_expr = "data/interim/expression",
  interim_meth = "data/interim/methylation",
  final_expr   = "data/final/expression",
  final_meth   = "data/final/methylation",
  pics         = "output/pics",
  tables       = "output/tables",
  logs         = "logs"
)
invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- 2. Датасеты ----------------------------------------------------------
# Экспрессия: 9 датасетов на платформах с официальной BioC-аннотацией.
EXPR_GSE <- c("GSE14642","GSE30453","GSE47353","GSE49058","GSE56033",
              "GSE56045","GSE56580","GSE65907","GSE68759")
METH_GSE <- c("GSE36054","GSE87571","GSE40279","GSE51032","GSE64495")

# ---- 3. Палитры -----------------------------------------------------------
EXPR_COL <- setNames(c("#b280a8","#998abd","#7799c4","#69a8bb","#7eb3a8",
                       "#8fba9b","#a2c28f","#b8ca84","#cdd37c"),
                     EXPR_GSE)
METH_COL <- setNames(c("#c4849a","#da7a8c","#e17a82","#e48178","#d8ad6a"),
                     METH_GSE)

# ---- 4. Платформы экспрессии (GSE -> GPL) и аннотация (GPL -> BioC-пакет) --
EXPR_PLATFORM <- c(GSE14642="GPL570", GSE30453="GPL5175", GSE47353="GPL6244",
                   GSE49058="GPL11532", GSE56033="GPL6244", GSE56045="GPL10558",
                   GSE56580="GPL10558", GSE65907="GPL10558", GSE68759="GPL11532")
BIOC_ANNO <- c(GPL570  = "hgu133plus2.db",
               GPL6244 = "hugene10sttranscriptcluster.db",
               GPL11532= "hugene11sttranscriptcluster.db",
               GPL5175 = "huex10sttranscriptcluster.db",
               GPL10558= "illuminaHumanv4.db")

# ---- 5. Метилирование: источники бет и метаданных -------------------------
# беты:  series — в *_series_matrix.txt.gz | suppl1 — один suppl | suppl2 — две части
METH_BETA <- list(
  GSE36054 = list(type="series", files="GSE36054_series_matrix.txt.gz"),
  GSE64495 = list(type="series", files="GSE64495_series_matrix.txt.gz"),
  GSE40279 = list(type="series", files="GSE40279_series_matrix.txt.gz"),
  GSE51032 = list(type="series", files="GSE51032_series_matrix.txt.gz"),
  GSE87571 = list(type="suppl2", files=c("GSE87571_matrix1of2.txt.gz","GSE87571_matrix2of2.txt.gz")),
  GSE55763 = list(type="suppl1", files="GSE55763_normalized_betas.txt.gz", stream=TRUE)
)
# метаданные берём из *_series_matrix.txt.gz. Для GSE87571 и GSE55763 положи в
# data/raw/methylation их МАЛЕНЬКИЕ серматрицы (28.5/49.1 КБ) — там метаданные.
# Если их нет: GSE87571 -> читаем xlsx; GSE55763 -> тянем метаданные сетью (см. utils_io).
AGE_IN_MONTHS <- c("GSE36054")          # возраст в месяцах -> /12

# ---- 6. Возрастные окна и бины --------------------------------------------
AGE_WINDOWS    <- list(teen=c(14,16), young=c(20,25), mid=c(37,40), older=c(58,62))
AGE_GRID       <- seq(5, 90, by = 1)
AGE_BIN_BREAKS <- seq(0, 120, by = 10)  # [0,10),[10,20),... только для визуализации (PCA)

# ---- 7. Пол: маркеры ------------------------------------------------------
SEX_MARKERS <- list(female = "XIST",
                    male   = c("RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY"))

# ---- 8. ПАМЯТЬ / ФИЛЬТРАЦИЯ (под 8 ГБ) ------------------------------------
N_VAR_CPG       <- 30000
N_VAR_GENE      <- Inf
MAX_PER_DATASET <- 3000                         # реально ограничивает только большие когорты
META_CAP_OVERRIDE <- list()                     # GSE55763 берём ПОЛНОСТЬЮ (2711) — двухпроходный slim-поток локально
NA_PROBE_MAX    <- 0.10
NA_KNN_MAX      <- 0.05

# ---- 9. Статистика / ковариаты / ComBat -----------------------------------
FDR <- 0.05; SMOOTH_K <- 5; MIN_SAMPLES_WIN <- 150
COMBAT_COVARS <- c("Age","Sex")                 # биологию сохраняем в mod (+ доли клеток)
N_PCA_FEATURES <- 3000                          # сабсэмпл фич для PCA-диагностики batch

# ---- 10. Утилиты + оптимизация памяти на старте ---------------------------
source(file.path("R","utils.R"))
source(file.path("R","utils_io.R"))
source(file.path("R","utils_analysis.R"))
setup_session()      # gc + потоки data.table + отчёт о памяти (в начале каждого скрипта)
