import os
import pysam
import csv
import argparse
import logging
import time

# Настройка логирования в файл
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("process.log"),  # Указываем файл для логов
        logging.StreamHandler()  # Дополнительно выводим логи в консоль
    ]
)

# Функция для отметок времени
def timestamp():
    return time.strftime('%Y-%m-%d %H:%M:%S')

# Путь к директории с файлами референсного генома
REF_GENOME_DIR = '/ref/GRCh38.d1.vd1_mainChr/sepChrs/'

# Функция для получения референсного аллеля по хромосоме и позиции
def get_ref_allele(chrom, pos):
    fasta_file = os.path.join(REF_GENOME_DIR, f'{chrom}.fa')
    
    # Проверяем существование файла
    if not os.path.exists(fasta_file):
        logging.error(f"Файл {fasta_file} не найден!")
        return None
    
    # Открываем файл с соответствующей хромосомой
    ref_genome = pysam.FastaFile(fasta_file)
    
    # Получаем аллель с позиции (учитываем 0-индексацию в pysam)
    try:
        ref_allele = ref_genome.fetch(chrom, pos - 1, pos)
    except Exception as e:
        logging.error(f"Ошибка при извлечении референсного аллеля: {e}")
        return None
    
    return ref_allele

# Основная функция для обработки файла
def process_file(input_file, output_file):
    try:
        with open(input_file, 'r') as infile:
            # Определяем формат конца строки
            first_line = infile.readline().strip()
            if 'CHROM' not in first_line or 'POS' not in first_line:
                raise ValueError("Неверный заголовок входного файла")
            
            infile.seek(0)  # Возвращаемся в начало файла после проверки заголовка
            reader = csv.DictReader(infile, delimiter='\t')
            
            # Открываем выходной файл для записи
            with open(output_file, 'w', newline='') as outfile:
                fieldnames = [field for field in reader.fieldnames if field not in ['allele1', 'allele2']] + ['REF', 'ALT']
                writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
                
                # Записываем заголовки
                writer.writeheader()
                
                # Обрабатываем каждую строку входного файла
                for row in reader:
                    chrom = row['#CHROM']
                    pos = int(row['POS'])
                    allele1 = row['allele1']
                    allele2 = row['allele2']
                    
                    # Получаем референсный аллель
                    ref_allele = get_ref_allele(chrom, pos)
                    if ref_allele is None:
                        logging.error(f"Не удалось получить референсный аллель для {chrom}:{pos}")
                        continue
                    
                    # Определяем, какой из аллелей референсный
                    if ref_allele == allele1:
                        row['REF'] = allele1
                        row['ALT'] = allele2
                    elif ref_allele == allele2:
                        row['REF'] = allele2
                        row['ALT'] = allele1
                    else:
                        row['REF'] = 'unknown'
                        row['ALT'] = 'unknown'
                    
                    row = {key: row[key] for key in fieldnames}  # Удаление allele1, allele2                    
                    
                    writer.writerow(row)
                logging.info(f"Обработка завершена. Результат записан в файл {output_file}")
    except FileNotFoundError:
        logging.error(f"Входной файл {input_file} не найден.")
    except ValueError as ve:
        logging.error(ve)
    except Exception as e:
        logging.error(f"Произошла ошибка: {e}")

# Функция для создания аргументов командной строки
def parse_args():
    parser = argparse.ArgumentParser(description="Скрипт для обработки файла SNPs и добавления референсных аллелей.")
    parser.add_argument('-i', '--input', required=True, help="Путь к входному файлу")
    parser.add_argument('-o', '--output', required=True, help="Путь к выходному файлу")
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    # Запуск обработки файла с временной меткой
    logging.info("Запуск скрипта...")
    logging.info(f"Входной файл: {args.input}")
    logging.info(f"Выходной файл: {args.output}")
    
    process_file(args.input, args.output)
    logging.info("Скрипт завершён.")
