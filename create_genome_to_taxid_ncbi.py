from openpyxl import load_workbook
import csv
import sys

num_samples = snakemake.config["num_samples"]

for i in range(num_samples):
    file = snakemake.input[i]
    
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

    print('Creating ' + snakemake.output[i])
    with open(snakemake.output[i], 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['genome_id', 'taxid'])
        for key in name_to_taxid.keys():
            writer.writerow([key, name_to_taxid[key]])
