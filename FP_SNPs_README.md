# Комманды для преподготовки файла 

#!/bin/bash

# Задаем имена файлов
input_file="FP_SNPs.txt"
output_file="FP_SNPs_10k_GB38_twoAllelsFormat.tsv"

# 1. Удаляем столбец GB37_position (3-й столбец) и сохраняем промежуточный файл
awk -F'\t' 'BEGIN {OFS="\t"} { $3=""; print $0 }' "$input_file" | sed 's/\t\t/\t/' > temp1.tsv

# 2. Меняем порядок колонок и сохраняем результат
awk -F'\t' 'BEGIN {OFS="\t"} { print $2, $4, $1, $5, $6 }' temp1.tsv > temp2.tsv

# 3. Переименовываем колонки, а также добавляем "chr" к значениям chromosome и "rs" к значениям rs#
awk -F'\t' 'BEGIN {OFS="\t"; print "CHROM", "POS", "ID", "allele1", "allele2"} 
     NR>1 { $1="chr"$1; $3="rs"$3; print $1, $2, $3, $4, $5 }' temp2.tsv > temp3.tsv

# 4. Удаляем строки, где CHROM имеет значение "X"
awk -F'\t' 'BEGIN {OFS="\t"} $1 != "chrX" {print $0}' temp3.tsv > "$output_file"

# 5. Удаляем промежуточные файлы
rm temp1.tsv temp2.tsv temp3.tsv

echo "Файл успешно создан: $output_file"
