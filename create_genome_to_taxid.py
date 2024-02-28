from openpyxl import load_workbook
from Bio import Entrez
import csv
from alive_progress import alive_bar
import time

Entrez.email = 'okb5109@psu.edu'
file = 'result.xlsx'

print('Reading ' + file)
workbook = load_workbook(filename=file)
sheet = workbook.active

output_names = [cell.value for cell in sheet['A']]
names = [cell.value.split()[0] for cell in sheet['A']]
batch_size = 500
stop = len(names)

taxids = []

print(f'Fetching taxids for {stop - 1} genomes')

print('Fetching UID numbers')
id_list = []
with alive_bar(stop - 1) as bar:
    for i in range(1, stop, 20):
        lst = names[i:i + ((stop - i) if (i + 20 > stop) else 20)]
        search = Entrez.esearch(db='assembly', term=' OR '.join(lst))
        for id in Entrez.read(search)['IdList']:
            id_list.append(id)
            bar()
        time.sleep(0.25)

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
            taxids.append(result['DocumentSummarySet']['DocumentSummary'][((len(id_list) %  batch_size) if (i + batch_size >  len(id_list)) else batch_size) - 1 - j]['Taxid'])
            bar()

print('Creating genome_to_taxid.tsv')
with open('genome_to_taxid.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(['genome_id', 'taxid'])
    with alive_bar(stop - 1) as bar:
        for i in range(1, stop):
            writer.writerow([output_names[i], taxids[i - 1]])
            bar()
