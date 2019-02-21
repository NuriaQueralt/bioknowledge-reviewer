# @name: graph.py
# @description: Module for graph building and management
# @version: 1.0
# @date: 16-02-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# All files before concatenation should have the same format:
# Edges format:
# 'subject_id',
# 'property_id',
# 'object_id',
# 'reference_uri',
# 'reference_supporting_text',
# 'reference_date',
# 'property_label',
# 'property_description',
# 'property_uri'

# Nodes format:
# 'id',
# 'semantic_groups',
# 'preflabel',
# 'synonyms',
# 'description'

# TODO: modularize from: http://localhost:8888/notebooks/workspace/ngly1-graph/curation/kylo/neo4j/networks/concatenate_network_files.ipynb
# TODO: each network type should be a class, then class methods: .get_nodes(), .get_edges(),...
# TODO: improve provenance, e.g.: adding to ref_supporting_text for monarch year, msigdb version, mondo version...

"""Module for functions to build the graph"""

import pandas as pd
import sys, os
import datetime

# VARIABLES
today = datetime.date.today()

# path to write data
path = os.getcwd() + "/graph"
if not os.path.isdir(path): os.makedirs(path)


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA

# BUILD NETWORK

# NETWORK MANAGEMENT FUNCTIONS


################ UTILS SECTION


def check_format(df,file_type='statements'):
    """This function checks if dataframe contains the expected columns before concatanation."""

    if file_type == 'concepts':
        try:
            df = df[['id', 'semantic_groups', 'preflabel', 'synonyms', 'description']]
        except:
            print('Concepts dataframe does not contain the expected columns. Raised error: ', sys.exc_info()[0])
            raise
        else:
            return df
    else:
        try:
            df = df[['subject_id', 'property_id', 'object_id', 'reference_uri',
            'reference_supporting_text', 'reference_date', 'property_label',
            'property_description', 'property_uri']]
        except:
            print('Statements dataframe does not contain the expected columns. Raised error: ', sys.exc_info()[0])
            raise
        else:
            return df

def print_graph(graph, filename):
    """This function save the graph into a CSV file."""

    # print output file
    path = os.getcwd() + '/graph'
    if not os.path.isdir(path): os.makedirs(path)
    pd.DataFrame(graph).to_csv('{}/{}_v{}.csv'.format(path,filename,today), index=False)
    # order columns
    #if file_type == 'concepts':
    #   df = pd.DataFrame(graph)[['id', 'semantic_groups', 'preflabel', 'synonyms', 'description']]
    # if file_type == 'statements':
    #   df = pd.DataFrame(graph)[['subject_id', 'property_id', 'object_id', 'reference_uri','reference_supporting_text', 'reference_date', 'property_label','property_description', 'property_uri']]
    #df.to_csv('{}/{}_v{}.csv'.format(path,filename,today), index=False)

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path,filename,today))


def graph_nodes():
    """This function generates graph nodes."""

    ## Edges
    #TODO: check with new files

    # load networks
    print('\nPreparing networks...')
    # curated_df = pd.read_csv('{}/curated_graph_edges_v2019-02-14.csv'.format(path))
    # monarch_path = os.getcwd() + "/monarch"
    # monarch_df = pd.read_table('{}/monarch_edges_v2019-02-15.tsv'.format(monarch_path))
    # rna = pd.read_csv('{}/rna_edges_v2019-01-25.csv'.format(path))
    # tf = pd.read_csv('{}/regulation_edges_v2019-01-29.csv'.format(path))
    print('Curated:')
    path = '/home/nuria/workspace/ngly1-graph/regulation/graph'
    curated_df = pd.read_csv('{}/curated_graph_edges_v2019-01-18.csv'.format(path))
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    monarch_path = '/home/nuria/workspace/ngly1-graph/monarch/1shell-animal/add-connections-to-net'
    monarch_df = pd.read_table('{}/monarch_edges_v2019-01-16.tsv'.format(monarch_path))
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    rna = pd.read_csv('{}/rna_edges_v2019-01-17.csv'.format(path))
    print(rna.shape)
    print(rna.columns)
    print('Regulatory:')
    tf = pd.read_csv('{}/regulation_edges_v2019-01-17.csv'.format(path), low_memory=False)
    print(tf.shape)
    print(tf.columns)

    # concat 1) curated 2) monarch 3) RNA-seq edges
    print('\nConcatenating into a graph...')
    statements = pd.concat([curated_df, monarch_df, rna], ignore_index=True, join="inner")
    print(statements.shape)

    # drop row duplicates
    print('\nDrop duplicated rows...')
    statements.drop_duplicates(keep='first', inplace=True)
    print(statements.shape)

    ## merge graph & tf
    # merge: 4 merges
    print('\nMerging tf-gene network to the graph...')
    # merge1: L_sub  &  tf_sub
    merge1 = pd.merge(statements, tf, how='inner', left_on='subject_id', right_on='subject_id',
                      suffixes=('_graph', '_tf'))

    # merge2: L_obj  &  tf_sub
    merge2 = pd.merge(statements, tf, how='inner', left_on='object_id', right_on='subject_id',
                      suffixes=('_graph', '_tf'))

    # merge3: L_sub  &  tf_obj
    merge3 = pd.merge(statements, tf, how='inner', left_on='subject_id', right_on='object_id',
                      suffixes=('_graph', '_tf'))

    # merge4: L_obj  &  tf_obj
    merge4 = pd.merge(statements, tf, how='inner', left_on='object_id', right_on='object_id',
                      suffixes=('_graph', '_tf'))

    # prepare merged edges: slice tf edges from merge
    # merge1
    merge1_clean = (merge1
    [['subject_id', 'property_id_tf', 'object_id_tf', 'reference_uri_tf',
      'reference_supporting_text_tf', 'reference_date_tf', 'property_label_tf',
      'property_description_tf', 'property_uri_tf']]
        .rename(columns={
        'property_id_tf': 'property_id',
        'object_id_tf': 'object_id',
        'reference_uri_tf': 'reference_uri',
        'reference_supporting_text_tf': 'reference_supporting_text',
        'reference_date_tf': 'reference_date',
        'property_label_tf': 'property_label',
        'property_description_tf': 'property_description',
        'property_uri_tf': 'property_uri'
    })
    )

    # merge2
    merge2_clean = (merge2
    [['subject_id_tf', 'property_id_tf', 'object_id_tf', 'reference_uri_tf',
      'reference_supporting_text_tf', 'reference_date_tf', 'property_label_tf',
      'property_description_tf', 'property_uri_tf']]
        .rename(columns={
        'subject_id_tf': 'subject_id',
        'property_id_tf': 'property_id',
        'object_id_tf': 'object_id',
        'reference_uri_tf': 'reference_uri',
        'reference_supporting_text_tf': 'reference_supporting_text',
        'reference_date_tf': 'reference_date',
        'property_label_tf': 'property_label',
        'property_description_tf': 'property_description',
        'property_uri_tf': 'property_uri'
    })
    )

    # merge3
    merge3_clean = (merge3
    [['subject_id_tf', 'property_id_tf', 'object_id_tf', 'reference_uri_tf',
      'reference_supporting_text_tf', 'reference_date_tf', 'property_label_tf',
      'property_description_tf', 'property_uri_tf']]
        .rename(columns={
        'subject_id_tf': 'subject_id',
        'property_id_tf': 'property_id',
        'object_id_tf': 'object_id',
        'reference_uri_tf': 'reference_uri',
        'reference_supporting_text_tf': 'reference_supporting_text',
        'reference_date_tf': 'reference_date',
        'property_label_tf': 'property_label',
        'property_description_tf': 'property_description',
        'property_uri_tf': 'property_uri'
    })
    )

    # merge4
    merge4_clean = (merge4
    [['subject_id_tf', 'property_id_tf', 'object_id', 'reference_uri_tf',
      'reference_supporting_text_tf', 'reference_date_tf', 'property_label_tf',
      'property_description_tf', 'property_uri_tf']]
        .rename(columns={
        'subject_id_tf': 'subject_id',
        'property_id_tf': 'property_id',
        'reference_uri_tf': 'reference_uri',
        'reference_supporting_text_tf': 'reference_supporting_text',
        'reference_date_tf': 'reference_date',
        'property_label_tf': 'property_label',
        'property_description_tf': 'property_description',
        'property_uri_tf': 'property_uri'
    })
    )

    ## concat merged edges to statements (<= curated+monarch+rna)
    # concat all 4 merges to merged edges
    merged = pd.concat([merge1_clean, merge2_clean, merge3_clean, merge4_clean], ignore_index=True, join="inner")

    # drop duplicates
    merged.drop_duplicates(inplace=True)
    print(merged.shape)

    # save graph
    print('\nSaving tf merged edges...')
    path = os.getcwd() + "/graph"
    merged.fillna('NA').to_csv('{}/regulation_graph_edges_v{}.csv'.format(path, today), index=False)

    # concat merged to statements
    statements = pd.concat([statements, merged], ignore_index=True, join="inner")
    print(statements.shape)

    # drop duplicates
    print('\nDrop duplicated rows...')
    statements.drop_duplicates(keep='first', inplace=True)
    print(statements.shape)

    ## Nodes
    # extracting nodes in the graph
    print('\nGenerating graph nodes...')
    st_nodes_l = pd.concat([statements.subject_id, statements.object_id], ignore_index=True)
    st_nodes_l.drop_duplicates(inplace=True)
    st_nodes_df = pd.DataFrame({'id': st_nodes_l})
    print(st_nodes_df.shape)

    return st_nodes_df

# BUILD GRAPH
def build_edges():
    """This function builds edges"""

    ## Edges
    #TODO: check with new files

    # load networks
    print('\nPreparing networks...')
    # curated_df = pd.read_csv('{}/curated_graph_edges_v2019-02-14.csv'.format(path))
    # monarch_path = os.getcwd() + "/monarch"
    # monarch_df = pd.read_table('{}/monarch_edges_v2019-02-15.tsv'.format(monarch_path))
    # rna = pd.read_csv('{}/rna_edges_v2019-01-25.csv'.format(path))
    # tf = pd.read_csv('{}/regulation_edges_v2019-01-29.csv'.format(path))
    print('Curated:')
    path = '/home/nuria/workspace/ngly1-graph/regulation/graph'
    curated_df = pd.read_csv('{}/curated_graph_edges_v2019-01-18.csv'.format(path))
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    monarch_df = pd.read_table('{}/monarch_edges_v2019-01-18.tsv'.format(path))
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    rna = pd.read_csv('{}/rna_edges_v2019-01-17.csv'.format(path))
    print(rna.shape)
    print(rna.columns)
    print('Regulatory:')
    tf_merged = pd.read_csv('{}/regulation_graph_edges_v2019-01-17.csv'.format(path))
    print(tf_merged.shape)
    print(tf_merged.columns)

    # concat 1) curated 2) monarch 3) RNA-seq edges
    print('\nConcatenating into a graph...')
    statements = pd.concat([curated_df, monarch_df, rna, tf_merged], ignore_index=True, join="inner")
    print(statements.shape)

    # drop row duplicates
    print('\nDrop duplicated rows...')
    statements.drop_duplicates(keep='first', inplace=True)
    print(statements.shape)

    # add property_uri for those without but with a curie property_id annotated
    curie_dct = {
        'ro': 'http://purl.obolibrary.org/obo/',
        'bfo': 'http://purl.obolibrary.org/obo/',
        'geno': 'http://purl.obolibrary.org/obo/',
        'dc': 'http://purl.org/dc/elements/1.1/',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
        'skos': 'http://www.w3.org/2004/02/skos/core#',
        'pato': 'http://purl.obolibrary.org/obo/',
        'sio': 'http://semanticscience.org/resource/',
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',
        'encode': 'https://www.encodeproject.org/search/?searchTerm='
    }
    for i, row in statements.iterrows():
        if ':' in str(row['property_uri']):
            property_uri = row['property_uri']
        elif ':' in str(row['property_id']):
            try:
                property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].replace(':',
                                                                                                                '_')
            except KeyError:
                property_uri = None
                print('There is a reference curie with and unrecognized namespace:', row['property_id'])
        else:
            property_uri = None
        statements.at[i, 'property_uri'] = property_uri

    # save graph
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    statements = statements[['subject_id', 'property_id', 'object_id', 'reference_uri',
                             'reference_supporting_text', 'reference_date', 'property_label',
                             'property_description', 'property_uri']]
    print(statements.shape)
    print(statements.columns)
    statements.fillna('NA').to_csv('{}/graph_edges_v{}.csv'.format(path, today), index=False)

    return statements


def build_nodes(statements):
    """This function builds graph nodes."""

    # TODO: check with new files

    # load networks
    # Note that i have to load the new files for curated and monarch nodes. That is because in the nb i included d2m and
    # g2p nodes (curated) and name attribute (monarch) in the fly.
    print('\nPreparing networks...')
    path = os.getcwd() + "/graph"
    curated_df = pd.read_csv('{}/curated_graph_nodes_v2019-02-14.csv'.format(path))
    monarch_path = os.getcwd() + "/monarch"
    # TODO: monarch is saved as csv but named as tsv=> fix monarch module
    monarch_df = pd.read_table('{}/monarch_nodes_v2019-02-15.tsv'.format(monarch_path), sep=',')
    #rna_df = pd.read_csv('{}/rna_nodes_v2019-01-25.csv'.format(path))
    #tf_df = pd.read_csv('{}/regulation_nodes_v2019-01-29.csv'.format(path))
    print('Curated:')
    path = '/home/nuria/workspace/ngly1-graph/regulation/graph'
    # curated_df = pd.read_csv('{}/curated_graph_nodes_v2019-01-18.csv'.format(path))
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    # monarch_df = pd.read_table('{}/monarch_nodes_v2019-01-18.tsv'.format(path))
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    rna_df = pd.read_csv('{}/rna_nodes_v2019-01-17.csv'.format(path))
    print(rna_df.shape)
    print(rna_df.columns)
    print('Regulatory:')
    tf_df = pd.read_csv('{}/regulation_nodes_v2019-01-17.csv'.format(path))
    print(tf_df.shape)
    print(tf_df.columns)

    ## Annotating nodes in the graph
    print('\nAnnotating nodes in the graph...')
    # extracting nodes in the graph
    st_nodes_l = pd.concat([statements.subject_id, statements.object_id], ignore_index=True)
    st_nodes_l.drop_duplicates(inplace=True)
    st_nodes_df = pd.DataFrame({'id': st_nodes_l})
    print('graph from e', st_nodes_df.shape)

    # annotating nodes
    curated_nodes = pd.merge(curated_df, st_nodes_df, how='inner', on='id')
    monarch_nodes = pd.merge(monarch_df, st_nodes_df, how='inner', on='id')
    rna_nodes = pd.merge(rna_df, st_nodes_df, how='inner', on='id')
    regulation_nodes = pd.merge(tf_df, st_nodes_df, how='inner', on='id')
    print('annotation check')
    print('curated', curated_nodes.shape)
    print('monarch', monarch_nodes.shape)
    print('rna', rna_nodes.shape)
    print('regulation', regulation_nodes.shape)

    # concat all, (importantly, concatenate first curated concepts with extended definitions)
    print('\nConcatenating all nodes...')
    nodes = pd.concat([curated_nodes, monarch_nodes, rna_nodes, regulation_nodes], ignore_index=True,
                      join="inner")
    print('graph ann', nodes.shape)
    diff = set(st_nodes_df.id) - set(nodes.id)
    print('diff', diff)

    # drop duplicated rows
    print('\nDrop duplicated rows...')
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    nodes.drop_duplicates(keep='first', inplace=True)
    print(nodes.shape)

    # drop duplicated nodes (keep first row (the curated), remove others (monarch))
    print('\nDrop duplicated nodes...')
    nodes.drop_duplicates(subset=['id'], keep='first', inplace=True)
    print(nodes.shape)

    # check
    if len(set(st_nodes_df.id)) != len(set(nodes.id)):
        print(
            '\nThere is a problem in the annotation of nodes.\nThe number of annotated nodes is different than the number of nodes in the graph.')
        print('Curated nodes not in the graph: {}'.format(set(curated_df.id) - set(curated_nodes.id)))
        print('Monarch nodes not in the graph: {}'.format(set(monarch_df.id) - set(monarch_nodes.id)))
        print('RNA-seq nodes not in the graph: {}'.format(set(rna_df.id) - set(rna_nodes.id)))
        print('Regulation nodes not in the graph: {}'.format(len(set(tf_df.id) - set(regulation_nodes.id))))
    else:
        print('\nAll graph nodes are annotated.')
        print('Regulation nodes not in the graph: {}'.format(len(set(tf_df.id) - set(regulation_nodes.id))))

    ## biothings
    # add attributes

    # all genes/proteins => add entrez|uniprot

    # save graph nodes
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    nodes = nodes[['id', 'semantic_groups', 'preflabel', 'synonyms', 'name', 'description']]
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    print(nodes.shape)
    print(nodes.columns)
    nodes.fillna('NA').to_csv('{}/graph_nodes_v{}.csv'.format(path, today), index=False)

    return nodes


# USER FUNCTIONS

def build(network_list):
    """This function formats, integrates and concats a list of networks passed."""

    # print of starting the process
    print('\nThe function "build()" is running, please keep calm and have some coffee...')

    # build the graph
    graph = dict()
    # for df in network_list:
    #     graph['edges'] = build_edges(df)
    #     graph['nodes'] = build_nodes(df)
    graph['edges'] = build_edges(df)
    graph['nodes'] = build_nodes(df)

    return graph

def edges(graph):
    """This function retrieves graph edges."""

    return graph.get('edges')


def nodes(graph):
    """This function retrieves graph edges."""

    return graph.get('nodes')


if __name__ == '__main__':
    #df = pd.read_table('./get-monarch-connections/monarch_connections.tsv')
    #monarch_edges(df)
    #nodes = monarch_nodes(df)
    #graph = build([df])
    #edges = edges(graph)
    #nodes = nodes(graph)
    #print_graph(nodes,'monarch_nodes')
    #print_graph(edges,'monarch_edges')

    # load networks and calculate graph nodes
    #graph_nodes_df = graph_nodes()
    edges = build_edges()
    nodes = build_nodes(edges)
    # check
    # print(graph_nodes_df.columns)
