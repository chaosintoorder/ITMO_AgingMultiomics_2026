library(readxl)
library(data.table)

# Путь к скачанному файлу (поправь, если лежит в другом месте)
excel_path <- "data/raw/genage_signature_raw.xlsx"

# 1. Читаем 2-й лист (overexpressed). Пропускаем 1 строчку "Table S3...", так как она сдвигает заголовки
up_genes <- as.data.table(read_excel(excel_path, sheet = 2, skip = 1))
up_genes[, dir := "up"] # Добавляем маркер направления

# 2. Читаем 3-й лист (underexpressed). Тоже пропускаем 1 строчку "Table S7..."
down_genes <- as.data.table(read_excel(excel_path, sheet = 3, skip = 1))
down_genes[, dir := "down"] # Добавляем маркер направления

# 3. Объединяем таблицы по строкам
genage_combined <- rbind(up_genes, down_genes)

# Чистим имена колонок (на всякий случай, убираем пробелы)
setnames(genage_combined, old = colnames(genage_combined), new = trimws(colnames(genage_combined)))

# 4. Сохраняем в папку, откуда её будет забирать функция genage_load()
# Обычно такие константы прописаны в PATHS$raw или PATHS$interim внутри config.R
# Например, сохраним в корень или в специальную подпапку:
output_dir <- "output/tables" # Измени на актуальный PATHS из твоего config.R
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fwrite(genage_combined, file.path(output_dir, "genage_signature.csv"))

cat("Готово. Объединено строк:", nrow(genage_combined), "\n")