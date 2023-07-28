import math
import pandas as pd

# data source: http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html - downloaded 2023-07-27

data = pd.read_excel('./raw_graphs/Cell_marker_Human.xlsx')

cellmarker_entries = []

# TODO: should we include tissue type as marker?
for i, row in data.iterrows():
    if math.isnan(row['GeneID']):
        continue
    if not isinstance(row['cellontology_id'], str):
        continue
    entry = {}
    entry['subject_name'] = row['Symbol']
    entry['subject_id'] = row['GeneID']
    entry['subject_id_prefix'] = 'NCBIGene'
    entry['subject_category'] = 'Gene'
    if row['cell_type'] == 'Normal cell':
        entry['object_id'] = row['cellontology_id']
        entry['object_name'] = row['cell_name']
        entry['object_id_prefix'] = 'CellOntology'
    elif row['cell_type'] == 'Cancer cell':
        entry['object_id'] = row['cellontology_id']
        entry['object_name'] = row['cancer_type'] + ' ' + row['cell_name']
        entry['object_id_prefix'] = 'CellOntology'
    entry['object_category'] = 'Cell'
    entry['knowledge_source'] = row['marker_source']
    entry['primary_knowledge_source'] = 'CellMarker' 
    entry['predicate'] = 'expressed_in'
    entry['anatomical_context_qualifier'] = row['uberonongology_id']
    if not math.isnan(row['PMID']):
        entry['publications'] = 'PMID:' + str(int(row['PMID']))
    else:
        entry['publications'] = ''
    cellmarker_entries.append(entry)

cellmarker_df = pd.DataFrame(cellmarker_entries)
cellmarker_df.to_csv('./processed_graphs/cellmarker.csv', index=False)
