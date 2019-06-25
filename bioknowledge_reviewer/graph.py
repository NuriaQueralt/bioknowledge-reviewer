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
import os
import datetime
from utils import *

# VARIABLES
today = datetime.date.today()

# path to write data
path = os.getcwd() + "/graph"
if not os.path.isdir(path): os.makedirs(path)


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA


# NETWORK MANAGEMENT FUNCTIONS

def print_graph(graph, filename):
    """
    This function saves the graph into a CSV file.
    :param graph: graph dictionary
    :param filename: file_name string
    :return: None object
    """

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


def graph_nodes(curation,monarch,transcriptomics,regulation,input_from_file=False):
    """
    This function generates graph nodes. The user can choose to input individual networks from file or \
    from the workflow.
    :param curation: curation graph edges
    :param monarch: monarch edges
    :param transcriptomics: rna graph edges
    :param regulation: regulation edges
    :param input_from_file: False (default value) or True
    :return: graph nodes list, regulation graph (merged) dataframe edges
    """

    print('\nThe function "graph_nodes()" is running...')

    ## Edges
    # load networks
    print('\nPreparing networks...')
    if input_from_file:
        if isinstance(curation,str) and isinstance(monarch,str) and isinstance(transcriptomics,str) and \
                isinstance(regulation,str):
            curated_df = get_dataframe_from_file(curation)
            monarch_df = get_dataframe_from_file(monarch)
            rna = get_dataframe_from_file(transcriptomics)
            tf = get_dataframe_from_file(regulation)
        else:
            print("Please, if you are providing the input from file then introduce the file path to the CSV file \
            , e.g. curation=str(/home/../file_name.csv). Otherwise, provide the objects and set the 'input_from_file' \
            argument to 'False'. Thanks!")
            raise
    else:
         curated_df = get_dataframe(curation)
         monarch_df = get_dataframe(monarch)
         rna = get_dataframe(transcriptomics)
         tf = get_dataframe(regulation)

    print('Curated:')
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    print(rna.shape)
    print(rna.columns)
    print('Regulatory:')
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
    print('\nThe regulation graph merged edges are saved at: {}/regulation_graph_edges_v{}.csv\n'.format(path, today))

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
    print('\nFinished graph_nodes().\n')

    return st_nodes_l, merged


# BUILD GRAPH

def build_edges(curation,monarch,transcriptomics,regulation,input_from_file=False):
    """
    This function builds the edges graph. The user can choose to input individual networks from file or \
    from the workflow.
    :param curation: curation graph edges object list
    :param monarch: monarch graph edges object list
    :param transcriptomics: rna graph edges object list
    :param regulation: regulation graph edges object list
    :param input_from_file: False (default value) or True
    :return: edges dataframe
    """

    print('\nThe function "build_edges()" is running...')
    ## Edges

    # load networks
    print('\nPreparing networks...')
    if input_from_file:
        if isinstance(curation,str) and isinstance(monarch,str) and isinstance(transcriptomics,str) and \
                isinstance(regulation,str):
            curated_df = get_dataframe_from_file(curation)
            monarch_df = get_dataframe_from_file(monarch)
            rna = get_dataframe_from_file(transcriptomics)
            tf_merged = get_dataframe_from_file(regulation)
        else:
            print("Please, if you are providing the input from file then introduce the file path to the CSV file \
            , e.g. curation=str(/home/../file_name.csv). Otherwise, provide the objects and set the 'input_from_file' \
            argument to 'False'. Thanks!")
            raise
    else:
         curated_df = get_dataframe(curation)
         monarch_df = get_dataframe(monarch)
         rna = get_dataframe(transcriptomics)
         tf_merged= get_dataframe(regulation)

    print('Curated:')
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    print(rna.shape)
    print(rna.columns)
    print('Regulatory:')
    print(tf_merged.shape)
    print(tf_merged.columns)

    # concat 1) curated 2) monarch 3) RNA-seq edges
    #TODO: check format
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
        elif ':' in str(row['property_id']) and str(row['property_id']).split(':')[0].lower() == 'skos':
            property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].split(':')[1]
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

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(statements.shape))
    print('* These are the edges attributes: {}'.format(statements.columns))
    print('* This is the first record:\n{}'.format(statements.head(1)))
    print('\nThe NGLY1 Deficiency knowledge graph edges are built and saved at:'
          ' {}/graph_edges_v{}.csv\n'.format(path, today))
    print('\nFinished build_edges().\n')

    return statements


def build_nodes(statements,curation,monarch,transcriptomics,regulation,input_from_file=False):
    """
    This function builds the nodes graph. The user can choose to input individual networks from file or \
    from the workflow.
    :param statements: graph edges dataframe
    :param curation: curation graph nodes object list
    :param monarch: monarch graph nodes object list
    :param transcriptomics: rna graph nodes object list
    :param regulation: regulation graph nodes object list
    :param input_from_file: False (default value) or True
    :return: nodes dataframe
    """

    print('\nThe function "build_nodes()" is running...')
    # load networks
    print('\nPreparing networks...')
    if input_from_file:
        if isinstance(curation,str) and isinstance(monarch,str) and isinstance(transcriptomics,str) and \
                isinstance(regulation,str):
            curated_df = get_dataframe_from_file(curation)
            monarch_df = get_dataframe_from_file(monarch)
            rna_df = get_dataframe_from_file(transcriptomics)
            tf_df = get_dataframe_from_file(regulation)
        else:
            print("Please, if you are providing the input from file then introduce the file path to the CSV file \
            , e.g. curation=str(/home/../file_name.csv). Otherwise, provide the objects and set the 'input_from_file' \
            argument to 'False'. Thanks!")
            raise
    else:
         curated_df = get_dataframe(curation)
         monarch_df = get_dataframe(monarch)
         rna_df = get_dataframe(transcriptomics)
         tf_df = get_dataframe(regulation)

    print('Curated:')
    print(curated_df.shape)
    print(curated_df.columns)
    print('Monarch:')
    print(monarch_df.shape)
    print(monarch_df.columns)
    print('Transcriptomics:')
    print(rna_df.shape)
    print(rna_df.columns)
    print('Regulatory:')
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
    #TODO: check format
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
            '\nThere is a problem in the annotation of nodes.\nThe number of annotated nodes '
            'is different than the number of nodes in the graph.')
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

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(nodes.shape))
    print('* These are the edges attributes: {}'.format(nodes.columns))
    print('* This is the first record:\n{}'.format(nodes.head(1)))
    print('\nThe NGLY1 Deficiency knowledge graph nodes are built and saved at: '
          '{}/graph_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished build_nodes().\n')

    return nodes


# USER FUNCTIONS

def _build(network_list):
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


def _edges(graph):
    """This function retrieves graph edges."""

    return graph.get('edges')


def _nodes(graph):
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
    # graph_nodes_df = graph_nodes()
    # print('graph nodes df:', graph_nodes_df.shape)
    curation_file = './graph/curated_graph_edges_v2019-03-04.csv'
    monarch_file = './monarch/monarch_edges_v2019-03-04.csv'
    rna_file = './graph/rna_edges_v2019-03-04.csv'
    tf_file = './graph/regulation_edges_v2019-03-04.csv'
    graph_nodes_list = graph_nodes(
        curation=curation_file,
        monarch=monarch_file,
        transcriptomics=rna_file,
        regulation=tf_file
    )
    print('graph nodes list len:', len(graph_nodes_list))
    print('graph nodes set len:', len(set(graph_nodes_list)))

    # build network
    # edges = build_edges()
    # nodes = build_nodes(edges)

    # check
    # print(graph_nodes_df.columns)
