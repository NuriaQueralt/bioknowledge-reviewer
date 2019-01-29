# @name: regulation.py
# @description: Module for regulatory network preparation and management
# @version: 1.0
# @date: 21-01-2019
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TODO: run R inside python to prepare tftargets data https://stackoverflow.com/questions/49922118/running-r-script-in-python

"""Module for the regulation data"""


import datetime
import json
import os #, glob, sys
import gseapy as gs
#import pandas as pd
#from gsheets import Sheets
#from biothings_client import get_client
#sys.path.insert(0,'/home/nuria/soft/utils3/lib/')
#import abravo_lib as utils

# VARIABLES
today = datetime.date.today()

# path to write data


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA
# TODO: check functions

def unique_list(dictionary, key, value):
    """This function adds non redundant values into a list value in a passed key dictionary."""
    
    genes = dictionary.get(key,set())
    for gene in value:
        genes.add(gene)
    dictionary[key] = genes
    
    return dictionary

def add_gene(dictionary, key, element):
    """This function prepares tf-gene_list dictionary: {symbol: [entrez]}."""

    aux = dictionary.get(key)
    aux = set(aux)
    aux.add(element)
    dictionary[key] = list(aux)

    return dictionary

# check functions
def format_exp(dictionary, key, element):
    """This function creates a dictionary where values are lists of items that can be redundant."""

    aux = dictionary.get(key, [])
    aux.append(element)
    dictionary[key] = aux

    return dictionary


def check_msigdb_geneset_name_format(data):
    """This function checks the format of gene set names.
    Nomenclature for TF binding sites sequences: TF symbol + :
        * _01: binding site id (consecutive number)
        * _Q2: binding site id (quality)
        * _Q6_01: binding site id
        * _DECAMER_Q2: binding site id
        * _B: binding site id (consensus)
    """

    maxlen = 0
    format_name = {}

    msigdb = {}
    for geneset_name in data:
        gs_name_v = geneset_name.split('_')
        format_name = format_exp(format_name, len(gs_name_v), geneset_name)
        if maxlen < len(gs_name_v):
            maxlen = len(gs_name_v)

    print('\n* Gene set name max length: {}'.format(maxlen))
    print('* Gene set name lengths: {}'.format(format_name.keys()))
    # TODO: the following print does not work because cannot be transformed into a df
    #print('This is the first record: {}'.format(pd.DataFrame(format_name).head(1)))

    #return

## Prepare tftargets data: from regulation/tftargets/save_data.R
# TODO: function to prepare tftargets data with R

#def prepare_tftargets_data():
#    """This functions reads Tf-gene raw networks of tftarget sources and writes them into JSON files."""


    #return

## Prepare msigdb data: from regulation/msigdb/exploration.ipynb
def prepare_msigdb_data():
    """This function prepares MSigDB data."""

    # path to write data
    path = os.getcwd() + '/regulation/msigdb/out'
    if not os.path.exists(path): os.makedirs(path)

    # load C3:TFT data
    gmt_path = '/home/nuria/workspace/ngly1-graph/regulation/msigdb/data/c3.tft.v6.1.entrez.gmt'
    data = gs.parser.gsea_gmt_parser("{}".format(gmt_path), min_size=1, max_size=10000)
    print('\n* Number of Transcription Factor Targets (TFT) gene sets: {}'.format(len(data.keys())))

    ## save raw data
    # TODO: save raw data
    #path = os.getcwd() + '/regulation/msigdb/data'
    #if not os.path.isdir(path): os.makedirs(path)
    #pd.DataFrame(data).to_csv('{}/c3.tft.v6.1.entrez.gmt'.format(path), index=False)

    ## prepare tf-gene_list dictionary: {symbol: [entrez]}: compile all gene set names into TF symbol
    # select tf:
    unknown = list()
    unrecognized = list()
    redundant_tf = list()
    tf_tfbs = {}
    msigdb = {}
    for geneset_name in data:
        gs_name_v = geneset_name.split('_')
        # remove 'unknown' (motif without a link to a tf)
        # gene symbols: Ideally, symbols should be no longer than six characters in length.
        # GCTNWTTGK_UNKNOWN
        if 'unknown' in gs_name_v[-1].lower():
            unknown.append(geneset_name)
            continue
        # LEN = 4 GATTGGY_NFY_Q6_01
        elif len(gs_name_v) == 4:
            tf = gs_name_v[1]
        # LEN = 3 two tfids: motif+TFACID or TFACID
        # GCCATNTTG_YY1_Q6
        elif len(gs_name_v) == 3 and len(gs_name_v[0]) >= 6:
            tf = gs_name_v[1]
        # AP4_Q6_01
        elif len(gs_name_v) == 3 and len(gs_name_v[0]) < 6:
            tf = gs_name_v[0]
        # LEN = 2 GFI1_01
        elif len(gs_name_v) == 2:
            tf = gs_name_v[0]
        else:
            unrecognized.append(geneset_name)
            continue

        # save tf-genelist
        if not msigdb.get(tf):
            msigdb[tf] = data[geneset_name]
        else:
            for gene in data[geneset_name]:
                msigdb = add_gene(msigdb, tf, gene)
            redundant_tf.append(tf)

        # save tf-tfbs
        tf_tfbs = format_exp(tf_tfbs, tf, geneset_name)

    # save msigdb raw network data
    with open('{}/tfid_genelist_entrez_msigdb.json'.format(path), 'w') as f:
        json.dump(data, f, sort_keys=True, indent=2)

    with open('{}/tf_genelist_entrez_msigdb.json'.format(path), 'w') as f:
        json.dump(msigdb, f, sort_keys=True, indent=2)

    with open('{}/tf_tfid_entrez.json'.format(path), 'w') as f:
        json.dump(tf_tfbs, f, sort_keys=True, indent=2)

    with open('{}/unknown_entrez_tf.tsv'.format(path), 'w') as f:
        f.write('\n'.join(unknown))

    with open('{}/unrecognized_entrez_tfid.tsv'.format(path), 'w') as f:
        f.write('\n'.join(unrecognized))

    with open('{}/tf_with_multiple_tfid_entrez.tsv'.format(path), 'w') as f:
        f.write('\n'.join(set(redundant_tf)))

    return data


## Prepare regulation data: regulation/regulation.ipynb
def load_tf_gene_edges():
    """This function loads raw networks from JSON files into dict variables."""

    # load json files data
    json_tftargets_path = '/home/nuria/workspace/ngly1-graph/regulation/tftargets/data'
    json_msigdb_path = os.getcwd() + '/regulation/msigdb/out'
    tred = json.load(open('{}/tred.json'.format(json_tftargets_path)))
    encode = json.load(open('{}/encode.json'.format(json_tftargets_path)))
    neph2012 = json.load(open('{}/neph2012.json'.format(json_tftargets_path)))
    trrust = json.load(open('{}/trrust.json'.format(json_tftargets_path)))
    msigdb = json.load(open('{}/tf_genelist_entrez_msigdb.json'.format(json_msigdb_path)))

    # normalize neph2012 data structure
    neph_all = dict()
    for cell in neph2012:
        for tf, genes in neph2012[cell].items():
            neph_all = unique_list(neph_all, tf, genes)
    neph = {key: list(neph_all[key]) for key in neph_all}

    # prepare data in a tuple
    data = (tred, encode, neph, trrust, msigdb)

    return data


# normalize to entrez, hgnc ids
def normalize_to_hgnc_id():
    """This function normalizes network gene IDs from symbol/entrez to HGNC."""

    # dict from symbol

    # dict from entrez

    # normalize to hgnc, else: entrez/symbol

    #return
    

# save edges
def prepare_data_edges():
    """This function prepares the expression dataset as edges."""

    #return data_edges


# prepare regulation edges to build the graph
def prepare_regulation_edges():
    """This function prepares and compiles all individual data edges into regulation edges to build the graph."""

    #return edges


# BUILD NETWORK

## build edges and nodes files
# build edges
def build_edges(edges):
    """This function builds the edges network file."""

    #return


# build nodes
def build_nodes(edges):
    """This function builds the nodes network file."""

    #return


# NETWORK MANAGEMENT FUNCTIONS


def print_nodes(nodes, filename):
    """This function save nodes into a CSV file."""

    # print output file

    #return

if __name__ == '__main__':
    # check TF name label format
    #data = prepare_msigdb_data()
    #check_msigdb_geneset_name_format(data)
    # prepare msigdb data
    prepare_msigdb_data()
    load_tf_gene_edges()
