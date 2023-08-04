# TODO: get info on available graphs
import os
import random

import igraph as ig
import pandas as pd

def load_graph(filename):
    # if filename not in files, try opening it as a path
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

def get_names_to_ids(graph):
    """Returns a dict mapping node names to IDs (ignoring prefixes and categories so on)"""
    names_to_ids = {}
    for v in graph.vs:
        names_to_ids[v['feature_name']] = v['name']
    return names_to_ids

def get_names_to_ids_networkx(graph):
    """Returns a dict mapping node names to IDs (ignoring prefixes and categories so on)"""
    names_to_ids = {}
    for n, attrs in graph.nodes.items():
        names_to_ids[attrs['name']] = n
    return names_to_ids


def get_category_ids_to_nodes(graph, category):
    """
    Returns a dict that maps from identifiers in the specified category to graph node indices, for graphs that are not SPOKE.
    """
    identifiers_to_ids = {}
    for v in graph.vs:
        if v['category'] == category:
            identifiers_to_ids[v['id']] = v.index
    return identifiers_to_ids

def largest_component(graph):
    "Returns a subgraph containing the largest connected component of the given graph."
    components = graph.connected_components()
    sizes = components.sizes()
    largest_component = 0
    largest_component_size = 0
    for i, c in enumerate(sizes):
        if c > largest_component_size:
            largest_component_size = c
            largest_component = i
    subgraph = components.subgraph(largest_component)
    return subgraph


def nodes_in_category(graph, category, attr_name='category'):
    "Returns all nodes that are within a given category, as a list of igraph.Vertex objects."
    nodes_in_category = []
    for v in graph.vs:
        attrs = v.attributes()
        if attr_name in attrs and attrs[attr_name] == category:
            nodes_in_category.append(v)
    return nodes_in_category


def nodes_in_category_networkx(graph, category, attr_name='category'):
    nodes_in_category = []
    for n, attrs in graph.nodes.items():
        if attr_name in attrs and attrs[attr_name] == category:
            nodes_in_category.append(n)
    return nodes_in_category

def random_nodes_in_category(graph, category, n_nodes):
    """
    Returns a list of random node ids in the given category.
    """
    nodes_in_category = []
    for v in graph.vs:
        attrs = v.attributes()
        if 'category' in attrs and attrs['category'] == category:
            nodes_in_category.append(attrs['name'])
    return random.sample(nodes_in_category, n_nodes)

def random_nodes(graph, n_nodes):
    """
    Returns a list of random node ids.
    """
    return random.sample([v['name'] for v in graph.vs], n_nodes)


def random_nodes_in_category_networkx(graph, category, n_nodes):
    """
    Returns a list of random spoke ids in the given category.
    """
    nodes_in_category = []
    for n, attrs in graph.nodes.items():
        if 'category' in attrs and attrs['category'] == category:
            nodes_in_category.append((n, attrs['identifier']))
    return random.sample(nodes_in_category, n_nodes)


def degree_sample(graph, node_list, n_samples, dist):
    """
    Degree-based node sampling, to sample nodes such that they approximately match the given degree distribution.

    Args:
        graph - an igraph.Graph
        node_list - a list of vertices to be sampled from
        n_samples - the number of points to sample
        dist - a distribution that has a pdf function (could use kernel density estimation?)
    """
    import numpy as np
    prob_vals = []
    for node in node_list:
        degree = graph.degree(node)
        prob_vals.append(dist.pdf(degree)[0])
    prob_vals = np.array(prob_vals)
    prob_vals = prob_vals/prob_vals.sum()
    sampled_nodes = set([])
    while len(sampled_nodes) < n_samples:
        node = random.choices(node_list, prob_vals)[0]
        if node not in sampled_nodes:
            sampled_nodes.add(node)
    return sampled_nodes
