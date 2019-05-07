# @name: curation.py
# @description: Module for curation network preparation and management
# @version: 1.0
# @date: 01-03-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TODO: function for version for wikibase (data will come from wikibase dump neo4j CSV)
"""Module for the curation data"""

import os,glob
import pandas as pd
from gsheets import Sheets
from biothings_client import get_client
#sys.path.insert(0,'/home/nuria/soft/utils3/lib/')
#import abravo_lib as utils
import utils
import datetime

# mondo annotation
#sys.path.insert(0,'/home/nuria/soft/utils3/ontologies')
import mondo_class as mondo


# VARIABLES
#timestamp
today = datetime.date.today()

# path to write data
path = os.getcwd() + "/curation"
if not os.path.isdir(path): os.makedirs(path)

# database version path
version='v20180118'

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


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA


def prepare_data_edges(curated_df):
    """
    This function pre-processes data edges from curation files retrieved from the biocurator: \
                  1) normalizes identifier CURIEs, and 2) uniformizes records attributes.

    :param curated_df: edges dataframe from curation edges files
    :return: pre-processed edges dataframe
    """

    print('\nThe function "prepare_data_edges()" is running...')
    print('\nPreparing curated network...')
    # concat curation edges
    #read_network()

    # ID curie normalization
    # subject_id
    curated_df['subject_id'] = (curated_df.subject_id
                                .apply(lambda x:
                                       'ClinVarVariant:50962' if 'HGVS' in str(x) else
                                       x.replace('Reactome', 'REACT') if 'Reactome' in str(x) else str(x).strip()
                                       )
                                )
    # object_id
    curated_df['object_id'] = (curated_df.object_id
                               .apply(lambda x:
                                      'ClinVarVariant:50962' if 'HGVS' in str(x) else
                                      x.replace('Reactome', 'REACT') if 'Reactome' in str(x) else str(x).strip()
                                      )
                               )

    # uniform fields
    curated_df = curated_df[['subject_id', 'property_id', 'object_id', 'reference_uri',
                             'reference_supporting_text', 'reference_date', 'property_label',
                             'property_description', 'property_uri']]

    # drop duplicates
    curated_df = curated_df.drop_duplicates()

    # save curated edges file at curation/
    #TODO: abstract this function
    print('\nSaving curated network at curation/...')
    path = os.getcwd() + "/curation"
    if not os.path.isdir(path): os.makedirs(path)
    curated_df.fillna('NA').to_csv('{}/curated_edges_v{}.csv'.format(path, today), index=False)
    print('\n*Curated edges data structure shape:', curated_df.shape)
    print('*Curated edges data structure fields:', curated_df.columns)
    print('\nThe curated edges are saved at: {}/curated_edges_v{}.csv\n'.format(path, today))
    print('\nFinished prepare_data_edges().\n')

    return curated_df


def prepare_data_nodes(curated_df):
    """
    This function pre-processes data nodes from curation files retrieved from the biocurator: \
                  1) normalizes identifier CURIEs, and 2) uniformizes records attributes.

    :param curated_df: nodes dataframe from curation nodes files
    :return: pre-processed nodes dataframe
    """

    print('\nThe function "prepare_data_nodes()" is running...')
    print('\nPreparing curated nodes...')
    # concat curation nodes
    #read_network()

    # ID curie normalization
    curated_df['id'] = (curated_df.id
                        .apply(lambda x:
                               'ClinVarVariant:50962' if 'HGVS' in str(x) else
                               x.replace('Reactome', 'REACT') if 'Reactome' in str(x) else str(x).strip()
                               )
                        )
    # uniform fields
    curated_df = curated_df[['id', 'semantic_groups', 'preflabel', 'synonyms', 'description']]

    # drop duplicates
    curated_df = curated_df.drop_duplicates()

    # save curated nodes file at curation/
    #TODO: abstract this function
    print('\nSaving curated nodes at curation/...')
    path = os.getcwd() + "/curation"
    if not os.path.isdir(path): os.makedirs(path)
    curated_df.fillna('NA').to_csv('{}/curated_nodes_v{}.csv'.format(path, today), index=False)
    print('\n*Curated nodes data structure shape:', curated_df.shape)
    print('*Curated nodes data structure fields:', curated_df.columns)
    print('\nThe curated nodes are saved at: {}/curated_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished prepare_data_nodes().\n')

    return curated_df


def read_data(csv_path, version):
    """
    This function imports curated network from edges and nodes csv files into dataframes.
    :param csv_path: string with the path to the network files, e.g '/home/nuria/workspace/ngly1-graph/curation'
    :param version: string with the files version, e.g. '2019-01-18'
    :return: edges and nodes dataframes read from files
    """

    print('\n Read data from curation/data/v2018 directory...')
    # read edges
    # csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/curated_v20180118'
    edges_df = pd.read_csv('{}/curated_edges_v{}.csv'.format(csv_path, version))
    print('\n* Number of curated edges:', len(edges_df))

    # read nodes
    nodes_df = pd.read_csv('{}/curated_nodes_v{}.csv'.format(csv_path, version))
    print('\n* Number of curated nodes:', len(nodes_df))

    return edges_df, nodes_df


def normalize_genes_to_graph(edges_df):
    """
    This function normalizes gene ID scheme. It performs Gene ID conversion \
    from curated to graph scheme. It replaces human entrez by HGNC ID, mouse entrez by MGI ID, worm entrez \
    by WormBase ID.
    :param edges_df: pre-processed curation edges dataframe
    :return: normalized edges dataframe, concept dictionary (where keys are all curated nodes)
    """

    print('\nMapping genes to HGNC ID...')
    ## GENES: normalize to HGNC
    # biothings api + dictionaries
    # concepts
    concept_dct = dict()
    for i, row in edges_df.iterrows():
        # node for subject
        concept_dct[row['subject_id']] = 1
        # node for object
        concept_dct[row['object_id']] = 1

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
    entrez = list(set(entrez))

    # api call
    mg = get_client('gene')
    df = mg.querymany(entrez, scopes='entrezgene', fields='HGNC', size=1, as_dataframe=True)

    # build dictionary
    ids = df.reset_index().rename(columns={'query': 'entrez'}).copy()
    entrez2hgnc_dct = dict(zip(ids.entrez, ids.HGNC))

    # map to hgnc
    lines = []
    for idx, row in edges_df.iterrows():
        # subject
        if ':' in row['subject_id']:
            if 'NCBIGene' in row['subject_id'].split(':')[0]:
                # human ncbi gene ids with HGNC ID
                if str(entrez2hgnc_dct[row['subject_id'].split(':')[1]]) != 'nan':
                    row['subject_id'] = "HGNC:" + entrez2hgnc_dct[row['subject_id'].split(':')[1]]
                # specific non human ncbi gene ids in the curated set
                elif row['subject_id'] == 'NCBIGene:173028':
                    row['subject_id'] = 'WormBase:WBGene00010160'
                elif row['subject_id'] == 'NCBIGene:11826':
                    row['subject_id'] = 'MGI:103201'

        # object
        if ':' in row['object_id']:
            if 'NCBIGene' in row['object_id'].split(':')[0]:
                # human ncbi gene ids with HGNC ID
                if str(entrez2hgnc_dct[row['object_id'].split(':')[1]]) != 'nan':
                    row['object_id'] = "HGNC:" + entrez2hgnc_dct[row['object_id'].split(':')[1]]
                # specific non human ncbi gene ids in the curated set
                elif row['object_id'] == 'NCBIGene:173028':
                    row['object_id'] = 'WormBase:WBGene00010160'
                elif row['object_id'] == 'NCBIGene:11826':
                    row['object_id'] = 'MGI:103201'

        lines.append((row))
    edges = pd.DataFrame.from_records(lines)

    return edges, concept_dct


def normalize_diseases_to_graph(edges_df):
    """
    This function normalizes disease ID scheme. It performs Disease ID mapping \
    from curated to graph scheme. It adds disease ID to MONDO ID edges.
    :param edges_df: pre-processed edges dataframe
    :return: normalized edges dataframe
    """


    ## DISEASES:
    print('\nAdding diseases to MONDO ID network...')
    # add d2m edges
    # add equivalentTo MONDO edges
    edges_l = list()
    for disease_id, mondo_id in d2m.items():
        edge = dict()
        edge['subject_id'] = disease_id
        edge['object_id'] = mondo_id
        edge['property_id'] = 'skos:exactMatch'
        edge['property_label'] = 'exact match'
        edge['property_description'] = 'NA'
        edge['property_uri'] = 'NA'
        edge['reference_uri'] = 'https://monarchinitiative.org/disease/' + mondo_id
        edge['reference_supporting_text'] = 'Manual extraction from Monarch Knowledge Graph.'
        edge['reference_date'] = '2018-04'
        edges_l.append(edge)

    d2m_edges_df = pd.DataFrame(edges_l)
    edges = pd.concat([edges_df, d2m_edges_df], ignore_index=True, join="inner")

    return edges


def normalize_genes_to_proteins_to_graph(curated_df,concept_dct):
    """
    This function normalizes genes to protein IDs. It performs Gene to Protein ID mapping \
    from curated to graph scheme. It adds gene ID to protein ID edges.
    :param curated_df: pre-processed edges dataframe
    :param concept_dct: concept dictionary
    :return: normalized edges dataframe
    """

    ## PROTEINS2GENES
    print('\nAdding gene to protein network...')
    # translate proteins to genes
    # biothings api + dictionaries
    # api input
    uniprot = list()
    for idx, row in concept_dct.items():
        if ':' in idx:
            if 'uniprot' in idx.split(':')[0].lower():
                uniprot.append(idx.split(':')[1])
    uniprot = list(set(uniprot))

    # api call
    mg = get_client('gene')
    df = mg.querymany(uniprot, scopes='uniprot', fields='HGNC', size=1, as_dataframe=True)

    # build dictionary
    ids = df.reset_index().rename(columns={'query': 'uniprot'}).copy()
    uniprot2hgnc_dct = dict(zip(ids.uniprot, ids.HGNC))

    # add equivalentTo edges
    edges_l = list()
    for uniprot, hgnc in uniprot2hgnc_dct.items():
        if str(uniprot2hgnc_dct[uniprot]) == 'nan':
            continue
        edge = dict()
        edge['subject_id'] = 'HGNC:' + hgnc
        edge['object_id'] = 'UniProt:' + uniprot
        edge['property_id'] = 'RO:0002205'
        edge['property_label'] = 'has gene product'
        edge['property_description'] = 'NA'
        edge['property_uri'] = 'NA'
        edge['reference_uri'] = 'http://mygene.info/clients/'
        edge['reference_supporting_text'] = 'Automatic extraction via the python client for mygene.info services.'
        edge['reference_date'] = today
        edges_l.append(edge)

    # add g2p network
    p2g_edges_df = pd.DataFrame(edges_l)
    curated_df = pd.concat([curated_df, p2g_edges_df], ignore_index=True, join="inner")

    # drop g2p duplicates
    print('\nDrop duplicated gene-protein relations...')
    # mark `has gene product` to be deleted if duplicated
    curated_df['g2p_mark'] = (
        [curated_df.at[idx, 'property_label'] if 'has gene product' in curated_df.at[idx, 'property_label'] else
         idx for idx in curated_df.index]
    )
    # keep first: keep the g2p manually added
    curated_df.drop_duplicates(subset=['subject_id', 'property_id', 'object_id', 'g2p_mark'], keep='first',
                               inplace=True)

    return curated_df


def prepare_curated_edges(edges_df):
    """
    This function prepares curated edges to the graph schema: 1) normalizes gene identifiers, 2) normalizes disease \
    identifiers, and 3) normalizes proteins to genes.
    :param edges_df: edges dataframe from the prepare_data_edges() function
    :return: edges dataframe
    """

    print('\nThe function "prepare_curated_edges()" is running...')
    print('\nPreparing curated edges to graph schema...')
    print('\nMapping genes to HGNC ID...')
    ## Normalize GENES: NCBI to HGNC
    edges, concept_dct = normalize_genes_to_graph(edges_df)

    ## Normalize DISEASES: add d2m network
    edges = normalize_diseases_to_graph(edges)

    ## Normalize PROTEINS2GENES: add g2p network
    edges = normalize_genes_to_proteins_to_graph(edges, concept_dct)
    print('\nFinished prepare_curated_edges().\n')

    return edges


def prepare_curated_nodes(curated_df):
    """
    This function prepares curated nodes to the graph schema: 1) normalizes gene identifiers, 2) normalizes disease \
    identifiers, and 3) normalizes proteins to genes.
    :param curated_df: nodes dataframe from the prepare_data_nodes() function
    :return: nodes dataframe
    """

    print('\nThe function "prepare_curated_nodes()" is running...')
    print('\nPreparing curated nodes to graph schema...')
    ## GENES: id normalization
    print('\nMapping genes to HGNC ID...')
    # biothings api + dictionaries
    # api input
    entrez = list()
    for i, row in curated_df.iterrows():
        if ':' in row['id']:
            if 'ncbigene' in row['id'].split(':')[0].lower():
                entrez.append(row['id'].split(':')[1])
    entrez = list(set(entrez))

    # api call
    print('\n* Querying BioThings to map Entrez gene IDs to HGNC IDs...')
    mg = get_client('gene')
    df = mg.querymany(entrez, scopes='entrezgene', fields='HGNC', size=1, as_dataframe=True)

    # build dictionary
    ids = df.reset_index().rename(columns={'query': 'entrez'}).copy()
    entrez2hgnc_dct = dict(zip(ids.entrez, ids.HGNC))

    # map to hgnc
    lines = []
    for idx, row in curated_df.iterrows():
        # subject
        if ':' in row['id']:
            if 'ncbigene' in row['id'].split(':')[0].lower():
                # human ncbi gene ids with HGNC ID
                if str(entrez2hgnc_dct[row['id'].split(':')[1]]) != 'nan':
                    row['id'] = "HGNC:" + entrez2hgnc_dct[row['id'].split(':')[1]]
                # specific non human ncbi gene ids in the curated set
                elif row['id'] == 'NCBIGene:173028':
                    row['id'] = 'WormBase:WBGene00010160'
                elif row['id'] == 'NCBIGene:11826':
                    row['id'] = 'MGI:103201'

        lines.append((row))
    curated_df = pd.DataFrame.from_records(lines)

    ## DISEASES: add mondo nodes
    # build dictionary with mondo nodes'description
    print('\nAdding diseases described by the MONDO ontology...')
    # import mondo owl terms
    # owl_f = '/home/nuria/workspace/ngly1-graph/ontologies/mondo.owl'
    owl_f = './ontologies/mondo.owl'
    tm = mondo.term(owl_f)

    # extract metadata from the mondo ontology
    nodes_l = list()
    for disease_id, mondo_id in d2m.items():
        mondo_term = tm.get_metadata_per_id(id=mondo_id)
        node = dict()
        node['id'] = mondo_term['id']
        node['semantic_groups'] = 'DISO'
        node['preflabel'] = mondo_term['label']
        node['synonyms'] = mondo_term['synonyms']
        node['description'] = mondo_term['definition']
        nodes_l.append(node)

    # add mondo nodes to curated_df
    d2m_nodes_df = pd.DataFrame(nodes_l)
    d2m_nodes_df.drop_duplicates(inplace=True)
    curated_df = pd.concat([curated_df, d2m_nodes_df], ignore_index=True, join="inner")

    ## ADD NAME ATTRIBUTE: gene name from BT, otherwise preflabel
    # biothings: annotate name to genes
    print('\nAdding Name attribute: gene names from BioThings...')
    # input: (preflabel) symbol,alias
    symbols = list(curated_df.preflabel)

    # query biothings
    print('\n* Querying BioThings to map gene symbols to name...')
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name', size=1, as_dataframe=True)

    # dictionary: {symbol:name}
    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    curated_s2n = dict(zip(ids.symbol, ids.name))

    # add name attribute
    curated_df['name'] = curated_df.preflabel.apply(lambda x: curated_s2n[x] if x in curated_s2n.keys() else x)

    ## PROTEINS: add and annotate encoding genes nodes of curated proteins
    # Annotating curated gene nodes from P-G edges
    print('\nPreparing encoding genes from ngly1 curated network...')
    # biothings api + dictionaries
    print('\nAdding BioThings annotation: gene symbol, name, synonyms, description...')
    # api input
    uniprot = set()
    for i, row in curated_df.iterrows():
        if ':' in row['id']:
            if 'uniprot' in row['id'].split(':')[0].lower():
                uniprot.add(row['id'].split(':')[1])

    # api call
    print('\n* Querying BioThings to map UniProt IDs to HGNC IDs, gene symbol, name, aliases, and description...')
    mg = get_client('gene')
    df = mg.querymany(uniprot, scopes='uniprot', fields='HGNC,symbol,name,alias,summary', size=1, as_dataframe=True)

    # build a list of nodes as list of dict, i.e a df, where a dict is a node
    nodes_l = list()
    for i, concept in df.iterrows():
        if str(concept['HGNC']) != 'nan':
            node = dict()
            node['id'] = 'HGNC:' + concept['HGNC']
            node['semantic_groups'] = 'GENE'
            node['preflabel'] = concept['symbol']
            node['name'] = concept['name']
            node['synonyms'] = '|'.join(list(concept['alias'])) if isinstance(concept['alias'], list) else concept[
                'alias']
            node['description'] = concept['summary']
            nodes_l.append(node)

    # structure as dataframe
    p2g_nodes_df = pd.DataFrame(nodes_l)
    p2g_nodes_df = p2g_nodes_df.fillna('NA')
    p2g_nodes_df.drop_duplicates(inplace=True)

    # add g2p network
    curated_df = pd.concat([curated_df, p2g_nodes_df], ignore_index=True, join="inner")

    # drop duplicates *new*
    curated_df = curated_df.drop_duplicates(subset=['id'], keep='first')
    print('\nFinished prepare_curated_nodes().\n')

    return curated_df


# BUILD NETWORK
def build_edges(edges_df):
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe from the prepare_curated_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    # give graph format
    edges_l = list()
    for row in edges_df.itertuples():
        edge = dict()
        edge['subject_id'] = row.subject_id
        edge['object_id'] = row.object_id
        edge['property_id'] = row.property_id
        edge['property_label'] = row.property_label
        edge['property_description'] = row.property_description
        edge['property_uri'] = row.property_uri
        edge['reference_uri'] = row.reference_uri
        edge['reference_supporting_text'] = row.reference_supporting_text
        edge['reference_date'] = row.reference_date
        edge['g2p_mark'] = row.g2p_mark
        edges_l.append(edge)

    # save curated graph edges to file at graph/
    #TODO: abstract this function
    print('\nSave curated graph edges file at graph/...')
    path = os.getcwd() + "/graph"
    if not os.path.isdir(path): os.makedirs(path)
    df = pd.DataFrame(edges_l)
    df = df[['subject_id','property_id','object_id','reference_uri','reference_supporting_text','reference_date', \
             'property_label','property_description','property_uri','g2p_mark']]
    df.fillna('NA').to_csv('{}/curated_graph_edges_v{}.csv'.format(path, today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe curation network edges are built and saved at: {}/curated_graph_edges_v{}.csv\n'.format(path, today))
    print('\nFinished build_edges().\n')

    return edges_l


def build_nodes(nodes_df):
    """
    This function builds the nodes network with the graph schema.
    :param nodes_df: nodes dataframe from the prepare_curated_nodes() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodges()" is running...')
    # give graph format
    nodes_l = list()
    for row in nodes_df.itertuples():
        node = dict()
        node['id'] = row.id
        node['semantic_groups'] = row.semantic_groups
        node['preflabel'] = row.preflabel
        node['name'] = row.name
        node['synonyms'] = row.synonyms
        node['description'] = row.description
        nodes_l.append(node)

    # save curated graph nodes to file at graph/
    #TODO: abstract this function
    print('\nSave curated graph nodes file at graph/...')
    path = os.getcwd() + "/graph"
    if not os.path.isdir(path): os.makedirs(path)
    df = pd.DataFrame(nodes_l)
    df = df[['id','semantic_groups','preflabel','synonyms','description','name']]
    df.fillna('NA').to_csv('{}/curated_graph_nodes_v{}.csv'.format(path, today), index=False)

    # print nodes info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe curation network nodes are built and saved at: {}/curated_graph_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished build_nodes().\n')

    return nodes_l


# NETWORK MANAGEMENT FUNCTIONS

def download_networks():
    """
    This function downloads curated network files (edges and nodes) from spreadsheets in Google Drive as csv files.
    """

    # create dir: curation/data
    data_path = os.getcwd() + '/curation/data/' + version
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
        '1pS3rT1hFShsCalu9bLGpAjQWp4y3_mgDxCtcHMAk2I8': 'ngly1_deficiency',
        '1z4PrO8AuNqyAOY3UMYxaMQ1JUcL5z-IbWsw54yQ8sYA': 'ngly1_human',
        '1thcKRGY1TnXepI8BJ6MYDhOgCayt3gEdPkOkFO4Nd4M': 'aqp1_human',
        '17kpta304URAxgd0NN4Uvyyu_qC1n2KObOrNAkFCQzRU': 'aqp1_mouse',
        '1ZF0cLyAN2_LPXbWNVss2yg2de7fvmHi21Zs7r1vFP54': 'glcnac_human',
        '17jXa5f_B74JaT8yuhExRNIozDRRNnPdcfM7yiV4BRtk': 'enns_2014',
        '1ZCfOdYtXn2mda2ybov8ibk0aXOtIBYDNw0RIoqkHldk': 'lam_2016'
    }
    spreadsheetIds_list = spreadsheetsIds_dct.keys()

    # download curated network files (edges and nodes)
    csv_name = lambda dct: '%s/%s.csv' % (data_path, dct.get('sheet'))
    #csv_name = lambda dct: '%s/%s-%s.csv' % (data_path, dct.get('title'), dct.get('sheet'))
    for sp_id in spreadsheetIds_list:
        driveFilesObject[sp_id].to_csv(make_filename=csv_name)

    return print('\nDownload process finished.\nFiles located at: {}.'.format(data_path))


def read_network(version=version):
    """
    This function concatenates and returns the curated network tsv files from drive as edges and nodes dataframes.
    :param version: string with data version
    :return: curated edges dataframe, curated nodes dataframe
    """

    print('\nThe function "read_network()" is running...')
    # concat all statements in the network
    print('\nReading and concatenating all curated statements in the network...')
    df_l = []
    #for file in glob.glob('{}/curation/data/*_edges.csv'.format(os.getcwd())):
    for file in glob.glob('{}/curation/data/{}/*_edges.tsv'.format(os.getcwd(), version)):
        with open(file, 'r') as f:
            #df_l.append(pd.read_table(f, sep=','))
            df_l.append(pd.read_table(f))

    network_df = pd.concat(df_l, ignore_index=True, join="inner")
    print('\n* Curation edge files concatenated shape:', network_df.shape)

    # concat all concepts in the network
    print('\nReading and concatenating all curated nodes in the network...')
    df_l = []
    for file in glob.glob('{}/data/{}/*_nodes.tsv'.format(path, version)):
        with open(file, 'r') as f:
            df_l.append(pd.read_table(f))

    nodes_df = pd.concat(df_l, ignore_index=True, join="inner")
    print('\n* Curation node files concatenated shape:', nodes_df.shape)
    print('\nFinished read_network().\n')

    return network_df, nodes_df


def _get_nodes_df(network_df):
    """
    This function extracts the network nodes from curated edges.
    :param network_df: curated edges dataframe
    :return: nodes dataframe
    """

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


def _normalize_nodes(nodes_df):
    """
    This function normalizes node IDs to Monarch scheme.
    :param nodes_df: nodes dataframe
    :return: normalized nodes dataframe
    """

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
    norm_nodes_df.to_csv('{}/curation/network_nodes.csv'.format(os.getcwd()), sep='\t', index=False, header=True)
    #print('\nNetwork nodes file created at: {}/curation/network_nodes.tsv\n'.format(os.getcwd()))

    # node counts:
    #print('Total nodes in the network: {}\n'.format(len(norm_nodes_df)))

    return norm_nodes_df


## First, get proteins as NCBIGene ID:
# 1. get uniprot_id list in my network
def _get_uniprot_list(df):
        """
        This function returns the UniProt ID list in the network.
        :param df: normalized nodes dataframe
        :return: uniprot list
        """

        uniprot_df = df[df.id_type == 'UniProt'].copy()
        uniprot_df['uniprot_id'] = (uniprot_df
                                        .loc[:, 'node_id']
                                        .apply(lambda x: x.split(':')[1])
                                        )
        uniprot_id_l = uniprot_df.uniprot_id
        # print(uniprot_df.head(2))
        return uniprot_id_l


# 2. get {'uniprot_id':'gene_id'} dictionary from biothings
def _get_uniprot2geneid_dict(uniprot_list):
    """
    This function returns the UniProt to NCBI gene ID dictionary from BioThings.
    :param uniprot_list: uniprot list
    :return: UniProt to Entrez ID dictionary
    """

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
def _map_uniprot2geneid(df, p2g_dct):
    """
    This function maps proteins (as UniProt ID) to genes (as NCBI Gene IDs).
    :param df: normalized nodes dataframe
    :param p2g_dct: UniProt to Entrez dictionary
    :return: nodes dataframe
    """

    nodes_df = df.reset_index(drop=True)
    nodes_df['monarch_id'] = (nodes_df
                                .node_id
                                .apply(lambda x: p2g_dct.get(x, None) if 'UniProt' in str(x) else x)
                            )
    #print(nodes_df.loc[nodes_df.node_type == 'protein'])

    # check node counts
    #print('\nTotal nodes in the network: {}'.format(len(nodes_df)))

    # save the protein normalized to gene id
    nodes_df.to_csv('{}/curation/network_nodes_monarch.csv'.format(os.getcwd()), sep='\t', index=False, header=True)
    #print('\nNetwork nodes mapped to Monarch data model saved at: {}/curation/network_nodes_monarch.tsv\n'.format(os.getcwd()))
    return nodes_df


## Second, create the list of network nodes to query monarch from 'monarch_id' column:
# 1. values in 'monarch_id': None, str(), list() -> normalize to str():
def _get_nodes_as_monarch(df):
    """
    This function returns nodes identified as Monarch can recognize: \
    no protein IDs and all identified using vocabularies used in Monarch.
    :param df: nodes dataframe
    :return: monarch nodes dataframe
    """

    # get rid of None values:
    nodes_monarch_df = df[['monarch_id']].copy()
    nodes_monarch_df = nodes_monarch_df.dropna()
    #print('\nTotal nodes in the Monarch-converted network: {}'.format(len(nodes_monarch_df)))

    return nodes_monarch_df


def _get_monarch_list(df):
    """
    This function returns a list of network nodes as Monarch expects from a dataframe to a list data structure.
    :param df: monarch nodes dataframe
    :return: monarch nodes list
    """

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


def _get_normalized_nodes(net_df):
    """
    This function returns a dataframe of all the nodes in the curated network \
    normalized to the ID schemes used in Monarch.
    :param net_df: curated edges dataframe
    :return: normalized nodes dataframe
    """

    # get the network dataframe
    # net_df = read_network()

    # get all network nodes as a dataframe
    nodes_df = _get_nodes_df(net_df)

    # normalize node ids to Monarch ID schemes
    norm_df = _normalize_nodes(nodes_df)

    return norm_df


def _get_proteins_as_ncbigenes(df):
    """
    This function returns all network proteins identified by their genes using NCBI Gene IDs as in Monarch.
    :param df: normalized nodes dataframe
    :return: nodes dataframe
    """

    # get uniprot list
    uniprot_id_list = _get_uniprot_list(df)

    # get uniprot2ncbigene dict from biothings
    uniprot2ncbigene_dict = _get_uniprot2geneid_dict(uniprot_id_list)

    # map protein ids to gene ids
    nodes_df = _map_uniprot2geneid(df, uniprot2ncbigene_dict)

    return nodes_df


def _get_list_of_monarch_id_nodes(df):
    """
    This function returns a list of nodes under Monarch ID scheme and all proteins/genes identified by NCBI Gene ID.
    :param df: nodes dataframe
    :return: nodes list
    """

    # get monarch id values
    monarch_df = _get_nodes_as_monarch(df)

    # normalize to string type
    nodes_list = _get_monarch_list(monarch_df)

    return nodes_list


def get_nodes(network_df):
    """
    This function gets nodes list to query Monarch API. It prepares all nodes from curation and normalizes \
    to Monarch ID scheme and translates all proteins to genes.
    :param network_df: curated edges dataframe
    :return: curated Monarch scheme nodes list
    """

    # get all network nodes normalized to ID schemes used in Monarch
    nodes_df = _get_normalized_nodes(network_df)

    # map all network proteins to their genes
    monarch_nodes_df = _get_proteins_as_ncbigenes(nodes_df)

    # get list of all network nodes in Monarch IDs
    node_list = _get_list_of_monarch_id_nodes(monarch_nodes_df)

    return node_list


if __name__ == '__main__':

    # prepare curated edges for monarch network interoperability
    # download curated network (web to local csv: edges and node descriptions)
    #download_networks()

    # read curated network(edges.csv)
    #curatedNetwork_df, nodes_df = read_network()

    # get nodes
    ##curatedNetworkNodes_df = get_nodes_df(curatedNetwork_df)
    ##curatedNetworkNodes_list = get_list_of_monarch_id_nodes(curatedNetworkNodes_df)
    #curatedNodes_list = get_nodes(curatedNetwork_df)

    ## build network for the graph
    ## read network from file (already preprocessed):
    ## prepare curated edges and nodes
    #csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/curated_v20180118'
    #edges_df, nodes_df = read_data(csv_path)
    #curated_graph_df = prepare_curated_edges(edges_df)
    #curated_graph_nodes_df = prepare_curated_nodes(nodes_df)

    ## build edges and nodes files
    #edges_l = build_edges(curated_graph_df)
    #nodes_l = build_nodes(curated_graph_nodes_df)

    ## checks
    #print('print edges', len(edges_l), len(edges_l[0].keys()))
    #print('print nodes', len(nodes_l), len(nodes_l[0].keys()))


    # build network for the graph
    # read network from drive sp version and concat all curated statements
    # curation_edges, curation_nodes = read_network(version='v20180118')

    # OR
    # download curated network (web to local csv: edges and node descriptions)
    download_networks()
    # read network from current drive version and concat all curated statements
    curation_edges, curation_nodes = read_network()

    # prepare data edges and nodes
    data_edges = prepare_data_edges(curation_edges)
    data_nodes = prepare_data_nodes(curation_nodes)

    # prepare curated edges and nodes
    curated_graph_df = prepare_curated_edges(data_edges)
    curated_graph_nodes_df = prepare_curated_nodes(data_nodes)


    # build edges and nodes files
    edges_l = build_edges(curated_graph_df)
    nodes_l = build_nodes(curated_graph_nodes_df)

    # checks
    print('print edges', len(edges_l), len(edges_l[0].keys()))
    print('print nodes', len(nodes_l), len(nodes_l[0].keys()))
