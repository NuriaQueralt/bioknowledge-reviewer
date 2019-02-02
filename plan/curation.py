# @name: curation.py
# @description: Module for curation network preparation and management
# @version: 1.0
# @date: 01-03-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TO DO:

"""Module for the curation data"""

import os, glob, sys
import pandas as pd
from gsheets import Sheets
from biothings_client import get_client
sys.path.insert(0,'/home/nuria/soft/utils3/lib/')
import abravo_lib as utils
import datetime

# mondo annotation
sys.path.insert(0,'/home/nuria/workspace/utils3/ontologies')
import mondo_class as mondo


# VARIABLES
today = datetime.date.today()

# path to write data
path = os.getcwd() + "/curation"
if not os.path.isdir(path): os.makedirs(path)


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA

# import curated statements as edges and nodes
def read_data():
    """This function imports curated statements as edges and nodes."""

    # read edges
    csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/curated_v20180118'
    edges_df = pd.read_csv('{}/curated_edges_v2019-01-18.csv'.format(csv_path))
    print('\n* Number of curated edges:', len(edges_df))

    # read nodes
    nodes_df = pd.read_csv('{}/curated_nodes_v2019-01-18.csv'.format(csv_path))
    print('\n* Number of curated nodes:', len(nodes_df))

    return edges_df, nodes_df


# normalize ID CURIEs and types
def prepare_data_edges(edges_df):
    """This function normalizes data edges."""

    ## GENES: normalize to HGNC
    # biothings api + dictionaries
    # concepts
    concept_dct = dict()
    for i, row in edges_df.iterrows():
        # node for subject
        concept_dct[row['subject_id']] = 1
        # node for object
        concept_dct[row['object_id']] = 1
    len(concept_dct.keys())

    # api input
    entrez = list()
    diseases = set()
    for idx, row in concept_dct.items():
        if ':' in idx:
            if 'ncbigene' in idx.split(':')[0].lower():
                entrez.append(idx.split(':')[1])
            elif 'doid' in idx.split(':')[0].lower() or 'omim' in idx.split(':')[0].lower() or 'orphanet' in \
                    idx.split(':')[0].lower():
                diseases.add(idx)
    print(len(entrez))
    entrez = list(set(entrez))
    print(len(entrez))
    print(len(diseases))
    print(diseases)

    # api call
    mg = get_client('gene')
    df = mg.querymany(entrez, scopes='entrezgene', fields='HGNC', size=1, as_dataframe=True)
    #df.head(2)

    # build dictionary
    ids = df.reset_index().rename(columns={'query': 'entrez'}).copy()
    # deal with no mappings
    # ids['HGNC'] = ids.HGNC.apply(lambda x: x if type(x) == str else 'NA')
    # dictionary
    entrez2hgnc_dct = dict(zip(ids.entrez, ids.HGNC))
    print(entrez2hgnc_dct['173028'])
    print(entrez2hgnc_dct['358'])

    # map to hgnc
    lines = []
    for idx, row in edges_df.iterrows():
        # subject
        if ':' in row['subject_id']:
            if 'NCBIGene' in row['subject_id'].split(':')[0]:
                if str(entrez2hgnc_dct[row['subject_id'].split(':')[1]]) != 'nan':
                    row['subject_id'] = "HGNC:" + entrez2hgnc_dct[row['subject_id'].split(':')[1]]

        # object
        if ':' in row['object_id']:
            if 'NCBIGene' in row['object_id'].split(':')[0]:
                if str(entrez2hgnc_dct[row['object_id'].split(':')[1]]) != 'nan':
                    row['object_id'] = "HGNC:" + entrez2hgnc_dct[row['object_id'].split(':')[1]]

        lines.append((row))
    edges = pd.DataFrame.from_records(lines)
    #edges.head(2)

    #return edges


# add d2m and g2p edges
def prepare_curated_edges(edges, nodes_df):
    """This function prepares curated edges to build the graph."""

    ## DISEASES:
    # add d2m edges
    # manually: dict diseases to mondo
    d2m = {
        'OMIM:223900': 'MONDO:0009131',
        'DOID:2476': 'MONDO:0019064',
        'Orphanet:869': 'MONDO:0009279',
        'DOID:11589': 'MONDO:0009131',
        'OMIM:614653': 'MONDO:0013839',
        'OMIM:615510': 'MONDO:0014219',
        'Orphanet:314381': 'MONDO:0013839',
        'DOID:10595': 'MONDO:0015626',
        'OMIM:608984': 'MONDO:0012166',
        'DOID:5212': 'MONDO:0015286',
        'OMIM:615273': 'MONDO:0014109',
        'DOID:0060308': 'MONDO:0019502',
        'DOID:0060728': 'MONDO:0014109',
        'OMIM:231550': 'MONDO:0009279',
        'DOID:0050602': 'MONDO:0009279'
    }

    # add equivalentTo MONDO edges
    edges_l = list()
    for disease, mondo in d2m.items():
        edge = dict()
        edge['subject_id'] = disease
        edge['object_id'] = mondo
        edge['property_id'] = 'skos:exactMatch'
        edge['property_label'] = 'exact match'
        edge['property_description'] = 'NA'
        edge['property_uri'] = 'NA'
        edge['reference_uri'] = 'NA'
        edge['reference_supporting_text'] = 'NA'
        edge['reference_date'] = 'NA'
        edges_l.append(edge)

    d2m_edges_df = pd.DataFrame(edges_l)
    #d2m_edges_df

    edges = pd.concat([edges, d2m_edges_df], ignore_index=True, join="inner")
    #edges.tail(2)

    # add mondo nodes
    # build dictionary with nodes'description
    # test mondo module
    owl_f = '/home/nuria/workspace/ngly1-graph/ontologies/mondo.owl'
    tm = mondo.term(owl_f)
    print(tm.get_metadata_per_id(id='MONDO:0015286'))

    # extract metadata from the mondo ontology
    nodes_l = list()
    for disease, mondo in d2m.items():
        mondo_term = tm.get_metadata_per_id(id=mondo)
        node = dict()
        node['id'] = mondo_term['id']
        node['semantic_groups'] = 'DISO'
        node['preflabel'] = mondo_term['label']
        # node['name'] = mondo_term['label']
        node['synonyms'] = mondo_term['synonyms']
        node['description'] = mondo_term['definition']
        nodes_l.append(node)

    # add to nodes_df
    d2m_nodes_df = pd.DataFrame(nodes_l)
    d2m_nodes_df.drop_duplicates(inplace=True)
    #d2m_nodes_df

    nodes = pd.concat([nodes_df, d2m_nodes_df], ignore_index=True, join="inner")
    #nodes.tail(2)

    #return edges, nodes



# BUILD NETWORK

def build_edges(edges):
    """This function builds the edges network file."""

    #return


def build_nodes(edges):
    """This function builds the nodes network file."""

    #return

# NETWORK MANAGEMENT FUNCTIONS

def download_networks():
    """This function downloads curated network files (edges and nodes) from spreadsheets in Google Drive as csv files."""

    # create dir: curation/data
    data_path = os.getcwd() + '/curation/data'
    if not os.path.isdir(data_path): os.makedirs(data_path)

    # download to dir edges and nodes
    ###
    # url-sheet1=https://docs.google.com/spreadsheets/d/1ZF0cLyAN2_LPXbWNVss2yg2de7fvmHi21Zs7r1vFP54/edit#gid=0
    # url-sheet2=https://docs.google.com/spreadsheets/d/1ZF0cLyAN2_LPXbWNVss2yg2de7fvmHi21Zs7r1vFP54/edit#gid=841865829
    # INSTALL: https://pypi.python.org/pypi/gsheets/0.3
    # https://pypi.python.org/pypi/pygsheets (more mature but under dev)
    # https://developers.google.com/sheets/api/quickstart/python
    ###
    # access to my drive files
    print('\nConnecting Google Drive to start the download process of the curated networks..\n')
    driveFilesObject = Sheets.from_files(
        '/home/nuria/client_secrets.json',
        '/home/nuria/storage.json'
    )

    # create a list for curated network files spreadsheets (with Google Drive API?)
    spreadsheetsIds_dct = {
        '1pS3rT1hFShsCalu9bLGpAjQWp4y3_mgDxCtcHMAk2I8': 'ngly1-deficiency',
        '1z4PrO8AuNqyAOY3UMYxaMQ1JUcL5z-IbWsw54yQ8sYA': 'ngly1_human',
        '1thcKRGY1TnXepI8BJ6MYDhOgCayt3gEdPkOkFO4Nd4M': 'aqp1_human',
        '17kpta304URAxgd0NN4Uvyyu_qC1n2KObOrNAkFCQzRU': 'aqp1_mouse',
        '1ZF0cLyAN2_LPXbWNVss2yg2de7fvmHi21Zs7r1vFP54': 'glcnac_human',
        '17jXa5f_B74JaT8yuhExRNIozDRRNnPdcfM7yiV4BRtk': 'enns_etal',
        '1ZCfOdYtXn2mda2ybov8ibk0aXOtIBYDNw0RIoqkHldk': 'lam_etal'
    }
    spreadsheetIds_list = spreadsheetsIds_dct.keys()

    # download curated network files (edges and nodes)
    csv_name = lambda dct: '%s/%s.csv' % (data_path, dct.get('sheet'))
    #csv_name = lambda dct: '%s/%s-%s.csv' % (data_path, dct.get('title'), dct.get('sheet'))
    for sp_id in spreadsheetIds_list:
        driveFilesObject[sp_id].to_csv(make_filename=csv_name)

    return print('\nDownload process finished.\nFiles located at: {}.'.format(data_path))


def read_network():
    """This function concatenates and returns the curated network as a dataframe."""

    # concat all statements in the network
    df_l = []
    for file in glob.glob('{}/curation/data/*_edges.csv'.format(os.getcwd())):
        with open(file, 'r') as f:
            df_l.append(pd.read_table(f, sep=',', ))

    network_df = pd.concat(df_l, ignore_index=True, join="inner")
    #print(network_df.shape)
    #print(network_df.head(2))

    return network_df


def get_nodes_df(network_df):
    """This function extracts the network nodes."""

    # optimize network dataframe
    net_df = (network_df
    [['subject_id', 'property_id', 'object_id', 'subject_qid', 'property_pid',
      'object_qid', 'subject_type', 'object_type', 'subject_term',
      'object_term']]
    )
    #print(net_df.head(2))

    # drop row duplicates
    #net_raw_df = net_df.copy()
    net_df = net_df.drop_duplicates()
    #print('net_raw: {}, net_clean: {}'.format(net_raw_df.shape, net_df.shape))

    # transform to three columns [sub + obj id], [sub + obj qid], [sub + obj type]
    sub_df = (net_df[['subject_id', 'subject_qid', 'subject_type']]
                .rename(
                columns={'subject_id': 'node_id', 'subject_qid': 'node_qid', 'subject_type': 'node_type'})
            )
    obj_df = (net_df[['object_id', 'object_qid', 'object_type']]
                .rename(
                columns={'object_id': 'node_id', 'object_qid': 'node_qid', 'object_type': 'node_type'})
            )
    net_nodes_df = pd.concat([sub_df, obj_df], ignore_index=True)
    #print('total nodes in the network with duplicates: {}'.format(len(net_nodes_df)))
    nodes_df = net_nodes_df.dropna(subset=['node_id']).drop_duplicates().copy()
    #print('total nodes in the network: {}'.format(len(nodes_df)))
    #print(nodes_df.head())

    return nodes_df


def normalize_nodes(nodes_df):
    """This function normalizes IDs to Monarch scheme."""

    # get node types
    #print(nodes_df.node_type.value_counts())

    # normalize types
    nodes_df['node_type_normalized'] = (nodes_df.node_type
        .apply(lambda x:
               'protein' if str(x) == "protein (enzymes)" else
               'phenotype' if 'clinical feature' in str(x) else
               'phenotype' if 'laboratory findings' in str(x) else
               'anatomy' if 'anatomy' in str(x) else
               x.strip().replace(' ', '_') if ' ' in str(x) else x)
        )
    #print(nodes_df.node_type_normalized.value_counts())

    # get namespace id types
    nodes_df['id_type'] = nodes_df.node_id.apply(lambda x: x.split(':')[0] if str(x).find(':') != -1 else x)
    #print(nodes_df.id_type.value_counts())

    # normalize namespace id types according to Monarch id schemes used
    nodes_df['node_id_normalized'] = (nodes_df.node_id
        .apply(lambda x:
               'CL:0000738' if str(x) == "FMA_62852; CL_0000738" else
               'ClinVarVariant:50962' if 'HGVS' in str(x) else
               x.replace('Reactome', 'REACT') if 'Reactome' in str(x) else str(x).strip())
        )
    nodes_df['id_type_normalized'] = nodes_df.node_id_normalized.apply(
        lambda x: x.split(':')[0] if str(x).find(':') != -1 else x)
    #print(nodes_df.id_type_normalized.value_counts())

    # save df to file normalized list of nodes
    norm_nodes_df = (nodes_df
        [['node_id_normalized', 'node_qid', 'node_type_normalized', 'id_type_normalized']]
        .rename(columns={'node_id_normalized': 'node_id', 'node_type_normalized': 'node_type', 'id_type_normalized': 'id_type'})
        .copy()
        )
    norm_nodes_df.to_csv('{}/curation/network_nodes.tsv'.format(os.getcwd()), sep='\t', index=False, header=True)
    #print('\nNetwork nodes file created at: {}/curation/network_nodes.tsv\n'.format(os.getcwd()))

    # node counts:
    #print('Total nodes in the network: {}\n'.format(len(norm_nodes_df)))

    return norm_nodes_df


    ## First, get proteins as NCBIGene ID:
    # 1. get uniprot_id list in my network
def get_uniprot_list(df):
        """This function returns a UniProt ID list in the network."""

        uniprot_df = df[df.id_type == 'UniProt'].copy()
        uniprot_df['uniprot_id'] = (uniprot_df
                                        .loc[:, 'node_id']
                                        .apply(lambda x: x.split(':')[1])
                                        )
        uniprot_id_l = uniprot_df.uniprot_id
        # print(uniprot_df.head(2))
        return uniprot_id_l


    # 2. get {'uniprot_id':'gene_id'} dictionary from biothings
def get_uniprot2geneid_dict(uniprot_list):
    """This function returns the UniProt to NCBI gene ID dictionary from BioThings."""

    mg = get_client('gene')
    r_df = mg.querymany(uniprot_list, scopes='uniprot', fields='entrezgene', as_dataframe=True)
    #print(r_df.head(2))

    # get the dictionary from the dataframe
    p2g = r_df[['_id']].dropna()
    p2g_dict = {}
    for idx in p2g.index:
        uniprot = 'UniProt:' + str(idx)
        value = p2g.at[idx, '_id']
        if not isinstance(value, str):
            for val in value.tolist():
                ncbigene = 'NCBIGene:' + str(val)
                p2g_dict = utils.add_elem_dictionary2(p2g_dict, uniprot, ncbigene)
        else:
            ncbigene = 'NCBIGene:' + str(value)
            p2g_dict = utils.add_elem_dictionary2(p2g_dict, uniprot, ncbigene)

    return p2g_dict


    # 3. map network nodes from uniprot (node_id) to ncbigene (monarch_id) and
    # add df column 'monarch_id', which will be the nodes_list to input monarch apis
def map_uniprot2geneid(df, p2g_dct):
    """This function maps proteins (as UniProt ID) to genes (as NCBI gene IDs)."""

    nodes_df = df.reset_index(drop=True)
    nodes_df['monarch_id'] = (nodes_df
                                .node_id
                                .apply(lambda x: p2g_dct.get(x, None) if 'UniProt' in str(x) else x)
                            )
    #print(nodes_df.loc[nodes_df.node_type == 'protein'])

    # check node counts
    #print('\nTotal nodes in the network: {}'.format(len(nodes_df)))

    # save the protein normalized to gene id
    nodes_df.to_csv('{}/curation/network_nodes_monarch.tsv'.format(os.getcwd()), sep='\t', index=False, header=True)
    #print('\nNetwork nodes mapped to Monarch data model saved at: {}/curation/network_nodes_monarch.tsv\n'.format(os.getcwd()))
    return nodes_df


    ## Second, create the list of network nodes to query monarch from 'monarch_id' column:
    # 1. values in 'monarch_id': None, str(), list() -> normalize to str():
def get_nodes_as_monarch(df):
    """This function retuns nodes identified as Monarch can recognize: no protein IDs and all identified using vocabularies used in Monarch."""

    # get rid of None values:
    nodes_monarch_df = df[['monarch_id']].copy()
    nodes_monarch_df = nodes_monarch_df.dropna()
    #print('\nTotal nodes in the Monarch-converted network: {}'.format(len(nodes_monarch_df)))

    return nodes_monarch_df


def get_monarch_list(df):
    """This function returns a list of network nodes as Monarch expects from a dataframe to a list data structure."""

    # normalize all types to str():
    nodes_l = []
    for value in df.monarch_id:
        if isinstance(value, list):
            for node in value:
                nodes_l.append(node)
        elif ':' not in value:
            continue
        else:
            nodes_l.append(value)
    #print(len(nodes_monarch_df.monarch_id), len(nodes_l), len(set(nodes_l)))

    # save network nodes list
    with open('{}/curation/network_nodes_monarch_list.tsv'.format(os.getcwd()), 'w') as f:
        f.write("monarch_id\n")
        # before debugging:     f.write("\n".join(nodes_l))
        # it is just a stetic change because the next script in the pipeline also set() this list..
        f.write("\n".join(set(nodes_l)))
    #print('\nNetwork nodes list to query Monarch API saved at: {}/curation/network_nodes_monarch_list.tsv'.format(os.getcwd()))
    #print('\nTotal nodes to query Monarch: {}'.format(len(set(nodes_l))))
    nodes_s = set(nodes_l)

    return list(nodes_s)


def get_normalized_nodes(net_df):
    """This function returns a dataframe of all the nodes in the curated network normalized to the ID schemes used in Monarch."""

    # get the network dataframe
    # net_df = read_network()

    # get all network nodes as a dataframe
    nodes_df = get_nodes_df(net_df)

    # normalize node ids to Monarch ID schemes
    norm_df = normalize_nodes(nodes_df)

    return norm_df


def get_proteins_as_ncbigenes(df):
    """This function returns all network proteins identified by their genes using ncbi gene ids as in Monarch."""

    # get uniprot list
    uniprot_id_list = get_uniprot_list(df)

    # get uniprot2ncbigene dict from biothings
    uniprot2ncbigene_dict = get_uniprot2geneid_dict(uniprot_id_list)

    # map protein ids to gene ids
    nodes_df = map_uniprot2geneid(df, uniprot2ncbigene_dict)

    return nodes_df


def get_list_of_monarch_id_nodes(df):
    """This function returns a list of nodes under Monarch ID schemes and all proteins/genes identified by NCBI gene ID."""

    # get monarch id values
    monarch_df = get_nodes_as_monarch(df)

    # normalize to string type
    nodes_list = get_monarch_list(monarch_df)

    return nodes_list


def get_nodes(network_df):
    """This function get nodes list to query Monarch API."""

    # get all network nodes normalized to ID schemes used in Monarch
    nodes_df = get_normalized_nodes(network_df)

    # map all network proteins to their genes
    monarch_nodes_df = get_proteins_as_ncbigenes(nodes_df)

    # get list of all network nodes in Monarch IDs
    node_list = get_list_of_monarch_id_nodes(monarch_nodes_df)

    return node_list


if __name__ == '__main__':

    # prepare curated edges for monarch network interoperability
    # download curated network (web to local csv: edges and node descriptions)
    #download_networks()

    # read curated network(edges.csv)
    #curatedNetwork_df = read_network()

    # get nodes
    ##curatedNetworkNodes_df = get_nodes_df(curatedNetwork_df)
    ##curatedNetworkNodes_list = get_list_of_monarch_id_nodes(curatedNetworkNodes_df)
    #curatedNodes_list = get_nodes(curatedNetwork_df)

    # build network for the graph
    #edges_df, nodes_df = read_data()
    #data_edges = prepare_data_edges(edges_df)
    #prepare_curated_edges(data_edges,nodes_df)