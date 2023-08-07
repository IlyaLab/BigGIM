import os

import igraph as ig
import pandas as pd

def load_graph(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError()
    if 'tsv' in filename:
        df = pd.read_csv(filename, sep='\t')
    else:
        df = pd.read_csv(filename)
    return df

def df_to_networkx(df, directed=False):
    """
    Converts a panda dataframe to a networkx Graph (or DiGraph), with node and edge attributes.
    """
    import networkx as nx
    create_using = nx.Graph
    if directed:
        create_using = nx.DiGraph
    df['subject_id_full'] = df['subject_id_prefix'] + '::' + df['subject_id'].astype(str)
    df['object_id_full'] = df['object_id_prefix'] + '::' + df['object_id'].astype(str)
    graph = nx.from_pandas_edgelist(df, source='subject_id_full', target='object_id_full',
            edge_attr=['predicate',
                       'Primary_Knowledge_Source',
                       'Knowledge_Source',
                       'publications'],
            create_using=create_using)
    node_attributes = {}
    for _, row in df.iterrows():
        if row['subject_id'] not in node_attributes:
            node_attributes[row['subject_id_full']] = {
                    'id': row['subject_id'],
                    'id_prefix': row['subject_id_prefix'],
                    'name': row['subject_name'],
                    'category': row['subject_category']}
        if row['object_id'] not in node_attributes:
            node_attributes[row['object_id_full']] = {
                    'id': row['object_id'],
                    'id_prefix': row['object_id_prefix'],
                    'name': row['object_name'],
                    'category': row['object_category']}
    nx.set_node_attributes(graph, node_attributes)
    return graph

def df_to_graph(df, directed=False):
    """
    Converts a panda dataframe to an igraph.Graph (or DiGraph), with node and edge attributes.
    """
    df['subject_id_full'] = df['subject_id_prefix'] + '::' + df['subject_id'].astype(str)
    df['object_id_full'] = df['object_id_prefix'] + '::' + df['object_id'].astype(str)
    # reorder?
    edges_df = df[['subject_id_full', 'object_id_full', 'predicate', 'Primary_Knowledge_Source', 'Knowledge_Source', 'publications']]
    graph = ig.Graph.DataFrame(edges_df, directed=directed, use_vids=False)
    node_attributes = {}
    for i, row in df.iterrows():
        if row['subject_id'] not in node_attributes:
            node_attributes[row['subject_id_full']] = {
                    'id': row['subject_id'],
                    'id_prefix': row['subject_id_prefix'],
                    'feature_name': row['subject_name'],
                    'category': row['subject_category']}
        if row['object_id'] not in node_attributes:
            node_attributes[row['object_id_full']] = {
                    'id': row['object_id'],
                    'id_prefix': row['object_id_prefix'],
                    'feature_name': row['object_name'],
                    'category': row['object_category']}
    for v in graph.vs:
        attributes = node_attributes[v['name']]
        v.update_attributes(attributes)
    return graph


def get_nodes_table(graph):
    """
    Returns a Pandas DataFrame of the nodes.
    """
    rows = []
    for v in graph.vs:
        row = {'id': v['name']}
        row.update(v.attributes())
        rows.append(row)
    return pd.DataFrame(rows)
