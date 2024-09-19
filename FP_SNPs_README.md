# Комманды для преподготовки файла 

#!/bin/bash

\# Задаем имена файлов
input_file="FP_SNPs.txt"
output_file="FP_SNPs_10k_GB38_twoAllelsFormat.tsv"

\# 1. Удаляем столбец GB37_position (3-й столбец) и сохраняем промежуточный файл
awk -F'\t' 'BEGIN {OFS="\t"} { $3=""; print $0 }' "$input_file" | sed 's/\t\t/\t/' > temp1.tsv

\# 2. Меняем порядок колонок и сохраняем результат
awk -F'\t' 'BEGIN {OFS="\t"} { print $2, $3, $1, $4, $5 }' temp1.tsv > temp2.tsv

\# 3. Переименовываем колонки, а также добавляем "chr" к значениям chromosome и "rs" к значениям rs#
awk -F'\t' 'BEGIN {OFS="\t"; print "CHROM", "POS", "ID", "allele1", "allele2"} 
     NR>1 { $1="chr"$1; $3="rs"$3; print $1, $2, $3, $4, $5 }' temp2.tsv > temp3.tsv

\# 4. Удаляем строки, где CHROM имеет значение "23" (Х-хромосома)
awk -F'\t' 'BEGIN {OFS="\t"} $1 != "chr23" {print $0}' temp3.tsv > "$output_file"

\# 5. Удаляем промежуточные файлы
rm temp1.tsv temp2.tsv temp3.tsv

echo "Файл успешно создан: $output_file"

# Скрипт для определения референсного аллеля (не удовлетворяет всем требованиям задания)
import os
import pysam
import csv

\# Директория с файлами FASTA для каждой хромосомы
ref_genome_dir = '~/bioinfo-docker/Genom'

\# Путь к входному файлу
input_file = '../FP_SNPs_10k_GB38_twoAllelsFormat.tsv'

\# Путь к выходному файлу
output_file = 'FP_SNPs_10k_GB38_with_REF_ALT.tsv'

\# Функция для получения референсного аллеля по хромосоме и позиции
def get_ref_allele(chrom, pos):
    \# Пути к файлам с хромосомами
    fasta_file = os.path.join(ref_genome_dir, f'{chrom}.fa')
    
    \# Проверяем существование файла
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Файл {fasta_file} не найден!")
    
    \# Открываем файл с соответствующей хромосомой
    ref_genome = pysam.FastaFile(fasta_file)
    
    \# Получаем аллель с позиции 
    \# В pysam 0-индексация
    ref_allele = ref_genome.fetch(f'{chrom}', pos - 1, pos)
    return ref_allele

\# Открываем входной файл и создаём выходной файл с добавленными колонками
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    
    \# Определяем поля для нового файла (добавляем REF и ALT)
    fieldnames = reader.fieldnames + ['REF', 'ALT']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    
    \# Записываем заголовки
    writer.writeheader()
    
    \# Обрабатываем каждую строку входного файла
    for row in reader:
        chrom = row['CHROM']
        pos = int(row['POS'])
        allele1 = row['allele1']
        allele2 = row['allele2']
        
        \# Получаем референсный аллель из файла референсного генома
        try:
            ref_allele = get_ref_allele(chrom, pos)
        except FileNotFoundError as e:
            print(e)
            continue
        
        \# Определяем, какой из аллелей является референсным
        if ref_allele == allele1:
            row['REF'] = 'allele1'
            row['ALT'] = 'allele2'
        elif ref_allele == allele2:
            row['REF'] = 'allele2'
            row['ALT'] = 'allele1'
        else:
            \# Если референсный аллель не совпал ни с одним, оставляем пустыми
            row['REF'] = 'unknown'
            row['ALT'] = 'unknown'
        
        \# Записываем строку в выходной файл
        writer.writerow(row)

print(f"Обработка завершена. Результат записан в файл {output_file}")



