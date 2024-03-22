from openpyxl import load_workbook
from Bio import Entrez
from alive_progress import alive_bar
import csv
import time

Entrez.email = 'okb5109@psu.edu'
file = 'result.xlsx'
batch_size = 100

print('Reading ' + file)
workbook = load_workbook(filename=file)
sheet = workbook.active

output_names = [cell.value for cell in sheet['A']]
names = [cell.value.split()[0] for cell in sheet['A']]

name_to_taxid = {}

print('Reading bacteria_hierarchy.tsv')
with open('bacteria_hierarchy.tsv', 'r') as ref:
    reader = csv.reader(ref, delimiter='\t', lineterminator='\n')
    i = 0
    for row in reader:
        if row[3][:3] in ['RS_', 'GB_']:
            name_to_taxid[row[3][3:]] = row[5]

skipped = []
for i in range(1, len(names)):
    if names[i] not in name_to_taxid:
        skipped.append(names[i])
print(f'Unable to find {len(skipped)} genomes in bateria_hierarchy.tsv')

print(f'Fetching taxids for {len(skipped)} genomes')

print('Fetching UID numbers')
id_list = []
with alive_bar(len(skipped) - 1) as bar:
    for i in range(0, len(skipped), 20):
        lst = skipped[i:i + ((len(skipped) - i) if (i + 20 > len(skipped)) else 20)]
        search = Entrez.esearch(db='assembly', term=' OR '.join(lst))
        for id in Entrez.read(search)['IdList']:
            id_list.append(id)
            bar()
        time.sleep(0.3)

time.sleep(3)
print('Fetching taxids')
with alive_bar(len(id_list)) as bar:
    for i in range(0, len(id_list), batch_size):
        lst = id_list[i:i + batch_size]
        search = Entrez.epost(db='assembly', id=','.join(lst))
        result = Entrez.read(search)
        webenv, query_key = result['WebEnv'], result['QueryKey']
        search = Entrez.efetch(db='assembly', rettype='docsum', retmode='xml', webenv=webenv, query_key=query_key)
        result = Entrez.read(search)
        for j in range(len(lst)):
            name_to_taxid[skipped[i + j]] = result['DocumentSummarySet']['DocumentSummary'][((len(id_list) %  batch_size) if (i + batch_size >  len(id_list)) else batch_size) - 1 - j]['Taxid']
            bar()

taxids = []
for i in range(1, len(names)):
    try:
        taxids.append(name_to_taxid[names[i]])
    except:
        print(f'Still could not find taxid for genome {names[i]}')
        exit(-1)

print('Creating genome_to_taxid.tsv')
with open('genome_to_taxid.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(['genome_id', 'taxid'])
    for i in range(1, len(names)):
        writer.writerow([output_names[i], taxids[i - 1]])