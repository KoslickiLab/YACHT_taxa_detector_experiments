from openpyxl import load_workbook
import csv

file = 'result.xlsx'

print('Reading ' + file)
workbook = load_workbook(filename=file)
sheet = workbook.active

name_to_taxid = {}

print('Reading ncbi_database_metadata.tsv')
with open('ncbi_database_metadata.tsv', 'r') as ref:
    reader = csv.reader(ref, delimiter='\t', lineterminator='\n')
    reader.__next__()
    for row in reader:
        name_to_taxid[row[0]] = row[1]

print('Creating genome_to_taxid.tsv')
with open('genome_to_taxid.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(['genome_id', 'taxid'])
    for key in name_to_taxid.keys():
        writer.writerow([key, name_to_taxid[key]])
