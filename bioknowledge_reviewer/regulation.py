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
import os
import gseapy as gs
from biothings_client import get_client
import gzip
import pandas as pd

# VARIABLES
today = datetime.date.today()

# path to write data
path = os.getcwd() + "/regulation"
if not os.path.isdir(path): os.makedirs(path)

# path to write graph
graph = os.getcwd() + '/graph'
if not os.path.isdir(graph): os.makedirs(graph)


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA
# TODO: check functions

def unique_list(dictionary, key, value):
    """
    This function adds non redundant values from an input list into the values list of a passed key and dictionary.
    :param dictionary: dictionary to update
    :param key: key string
    :param value: input list
    :return: updated dictionary
    """
    
    genes = dictionary.get(key,set())
    for gene in value:
        genes.add(gene)
    dictionary[key] = genes
    
    return dictionary


def add_gene(dictionary, key, element):
    """
    This function prepares tf-gene_list dictionary: {symbol: [entrez]}.
    :param dictionary: dictionary to update
    :param key: key string
    :param element: value string
    :return: updated dictionary
    """

    aux = dictionary.get(key)
    aux = set(aux)
    aux.add(element)
    dictionary[key] = list(aux)

    return dictionary


def add_elem_dictionary(dictionary, key, elem, repet = False):
    """
    This function adds elements to a dictionary key. The value list can be set \
    whether to include redundant values (repet=True) or not (repet=False).
    :param dictionary: dictionary to update
    :param key: key string
    :param elem: value string
    :param repet: False (default value) or True
    :return: updated dictionary
    """

    if key in dictionary:
        aux = dictionary.get(key)
        if repet:
            aux.append(elem)
            dictionary[key] = aux
        else:
            if not elem in aux:
                aux.append(elem)
                dictionary[key] = aux
    else:
        dictionary[key] = [elem]

    return dictionary


# check functions
def format_exp(dictionary, key, element):
    """
    This function creates a dictionary where values are lists of items that can be redundant.
    :param dictionary: dictionary to update
    :param key: key string
    :param element: value string
    :return: updated dictionary
    """

    aux = dictionary.get(key, [])
    aux.append(element)
    dictionary[key] = aux

    return dictionary


def check_msigdb_geneset_name_format(data):
    """
    This function checks the format of gene set names.
    Nomenclature for TF binding sites sequences: TF symbol + :
        * _01: binding site id (consecutive number)
        * _Q2: binding site id (quality)
        * _Q6_01: binding site id
        * _DECAMER_Q2: binding site id
        * _B: binding site id (consensus)
    :param data: MSigDB raw gene_set_name-geneset edges data
    :return: None object
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


## Prepare individual raw TF-gene networks
# Prepare tftargets data: from regulation/tftargets/save_data.R
# TODO: function to prepare tftargets data with R

#def prepare_tftargets_data():
#    """This functions reads Tf-gene raw networks of tftarget sources and writes them into JSON files."""


    #return


def prepare_msigdb_data(gmt_path):
    """
    This function prepares MSigDB raw data for TF-gene integration.
    It saves the MSigDB raw TF-gene network as JSON file.
    :param gmt_path: path to C3:TFT entrez GMT data file string
    :return: None object
    """

    print('\nThe function "prepare_msigdb_data()" is running...')
    # path to write data
    path = os.getcwd() + '/regulation/msigdb/out'
    if not os.path.exists(path): os.makedirs(path)

    # load C3:TFT data
    # gmt_path = '/home/nuria/workspace/ngly1-graph/regulation/msigdb/data/c3.tft.v6.1.entrez.gmt'
    data = gs.parser.gsea_gmt_parser("{}".format(gmt_path), min_size=1, max_size=10000)
    print('\n* Number of Transcription Factor Targets (TFT) gene sets: {}'.format(len(data.keys())))

    ## save raw data
    # TODO: save raw data as gmt (data is of type dict())
    # data_path = os.getcwd() + '/regulation/msigdb/data'
    # if not os.path.isdir(data_path): os.makedirs(data_path)
    # pd.DataFrame(data).to_csv('{}/c3.tft.v6.1.entrez.gmt'.format(data_path), index=False)

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
    print('\nThe MSigDB raw network is saved at: {}/tf_genelist_entrez_msigdb.json. '
          'Other reporting files are also saved at the same directory.\n'.format(path))
    print('\nFinished prepare_msigdb_data().\n')
    #
    # return data


## Prepare regulation data
def load_tf_gene_edges():
    """
    This function loads individually each database raw network from JSON files into a dict variable.
    Each dictionary contains TFs as keys and target genes list as values.
    :return: tred, encode, neph, trrust and msigdb individual raw networks dictionaries (in this order) as tuple
    """

    print('\nThe function "load_tf_gene_edges()" is running...')
    # load json files data
    #json_tftargets_path = '/home/nuria/workspace/ngly1-graph/regulation/tftargets/data'
    json_tftargets_path = './regulation/tftargets/data'
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
    print('\nFinished load_tf_gene_edges().\n')

    return data


# normalize to entrez, hgnc ids
def get_gene_id_normalization_dictionaries(data):
    """
    This function gets gene ID dictionaries to normalize network gene IDs to gene symbol, entrez and HGNC IDs. \
    In the raw networks there is a mismatch of ID schemes:
        * **symbols**: TF, trrust.genes
        * **entrez**: tred.genes, encode.genes, neph.genes, msigdb.genes
    :param data: (tred, encode, neph, trrust, msigdb) individual raw networks dictionaries (in this order) as tuple \
    from the load_tf_gene_edges() function
    :return: (symbol2entrez_dict, symbol2hgnc_dict, entrez2hgnc_dict, entrez2symbol_dict) dictionaries tuple
    """

    print('\nThe function "get_gene_id_normalization_dictionaries()" is running...')
    # individual raw networks
    tred = data[0]
    encode = data[1]
    neph = data[2]
    trrust = data[3]
    msigdb = data[4]

    ## dict from symbol
    # symbols input list
    trrust_genes = {gene for tf in trrust for gene in trrust.get(tf)}
    symbols = list(tred.keys() | encode.keys() | neph.keys() | trrust.keys() | msigdb.keys() | trrust_genes)
    #print('symbols:', len(symbols))

    ## ID dictionaries: symbols2entrez and symbol2hgnc
    # query biothings for symbol2entrez and symbol2hgnc
    print('\n* Querying BioThings to map gene symbols to HGNC and Entrez IDs...')
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='entrezgene,HGNC', size=1, as_dataframe=True)
    #print('symbols to entrez/hgnc: ',df.shape)

    # not found
    missing = df[['notfound']].copy()
    missing = missing.reset_index().rename(columns={'query': 'symbol'})
    missing = missing[missing['notfound'] == True][['symbol']]
    missing_symbol_l = [symbol for symbol in missing['symbol']]

    # save not found
    #print('\n * Missing symbols:', len(missing_symbol_l))
    with open('{}/not_found_symbols.list'.format(path), 'w') as f:
        f.write('\n'.join(missing_symbol_l))
    print('\nSaving not found gene symbols at: {}/not_found_symbols.list\n'.format(path))

    # prepare ids dataframe for dictionary construction
    ids = (df.reset_index()
           .rename(columns={'query': 'symbol', 'HGNC': 'hgnc', 'entrezgene': 'entrez'})
           [['symbol', 'hgnc', 'entrez']]
           .copy()
           )
    #print('ids:',ids.shape)

    # build dictionaries
    symbol2entrez_dict = dict(zip(ids.symbol, ids.entrez))
    symbol2hgnc_dict = dict(zip(ids.symbol, ids.hgnc))

    # not null ID allowed: associate symbol for those genes without entrez: symbol => entrez > symbol
    # add namespaces
    for symbol in symbol2entrez_dict:
        if isinstance(symbol2entrez_dict[symbol], float):
            symbol2entrez_dict[symbol] = symbol
        else:
            entrez = symbol2entrez_dict[symbol]
            symbol2entrez_dict[symbol] = "NCBIGene:" + entrez

    # not null ID allowed: associate entrez else symbol for those genes without hgnc: symbol =>  hgnc > entrez > symbol
    # add namespaces
    for symbol in symbol2hgnc_dict:
        if isinstance(symbol2hgnc_dict.get(symbol), float):
            symbol2hgnc_dict[symbol] = symbol2entrez_dict[symbol]
        else:
            hgnc = symbol2hgnc_dict[symbol]
            symbol2hgnc_dict[symbol] = "HGNC:" + hgnc


    ## dict from entrez
    # entrez input list
    tred_genes = {gene for tf in tred for gene in tred.get(tf)}
    encode_genes = {gene for tf in encode for gene in encode.get(tf)}
    neph_genes = {gene for tf in neph for gene in neph.get(tf)}
    msigdb_genes = {gene for tf in msigdb for gene in msigdb.get(tf)}
    entrez = list(tred_genes | encode_genes | neph_genes | msigdb_genes)
    #print('entrez:', len(entrez))

    ## ID dictionaries: entrez2hgnc, entrez2symbol
    # query biothings for entrez2hgnc, entrez2symbol
    print('\n* Querying BioThings to map Entrez to HGNC IDs and gene symbols...')
    mg = get_client('gene')
    df = mg.querymany(entrez, scopes='entrezgene', fields='HGNC,symbol', size=1, as_dataframe=True)
    #print('entrez to hgnc/symbol: ', df.shape)

    # not found
    missing = df[['notfound']].copy()
    missing = missing.reset_index().rename(columns={'query': 'entrez'})
    missing = missing[missing['notfound'] == True][['entrez']]
    missing_entrez_l = [entrez for entrez in missing['entrez']]

    # save not found
    #print('\n * Missing entrez:', len(missing_entrez_l))
    with open('{}/not_found_entrez.list'.format(path), 'w') as f:
        f.write('\n'.join(missing_entrez_l))
    print('\nSaving not found Entrez gene IDs at: {}/not_found_entrez.list\n'.format(path))

    # prepare ids dataframe for dictionary construction
    ids = (df.reset_index()
           .rename(columns={'query': 'entrez', 'HGNC': 'hgnc'})
           [['entrez', 'hgnc', 'symbol']]
           .copy()
           )
    # add namespaces
    ids['entrez_id'] = ids.entrez.apply(lambda x: "NCBIGene:" + x if type(x) == str else x)
    ids['hgnc_id'] = ids.hgnc.apply(lambda x: "HGNC:" + x if type(x) == str else x)
    #print('ids:', ids.shape)

    # build dictionaries
    entrez2hgnc_dict = dict(zip(ids.entrez_id, ids.hgnc_id))
    entrez2symbol_dict = dict(zip(ids.entrez_id, ids.symbol))

    # not null ID allowed: associate entrez for those genes without hgnc: entrez=> hgnc > entrez
    for entrez in entrez2hgnc_dict:
        if isinstance(entrez2hgnc_dict.get(entrez), float):
            entrez2hgnc_dict[entrez] = entrez

    # not null ID allowed: associate symbol, else 'NA' because it's a dict to annotate: entrez=> symbol or 'NA'
    for entrez in entrez2symbol_dict:
        if isinstance(entrez2symbol_dict.get(entrez), float):
            entrez2symbol_dict[entrez] = 'NA'

    dicts = (symbol2entrez_dict, symbol2hgnc_dict, entrez2hgnc_dict, entrez2symbol_dict)
    print('\nFinished get_gene_id_normalization_dictionaries().\n')

    return dicts


# save edges
def prepare_data_edges(data,dicts):
    """
    This function prepares each individual regulatory dataset as edges. It normalizes and stores them \
     separately as tftargets and msigdb edges.
    Edges (network) data structure:

    | source | dataset | tf_source_id | gene_source_id | source_uri | s_entrez_id | s_hgnc_id | s_symbol | p_id | \
    p_label | o_entrez_id | o_hgnc_id | o_symbol | reference_id | reference_date |
    |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

    1. `source`: str tftargets
    2. `dataset`: str tred
    3. `tf_source_id`: str transcriptor factor tftargets id
    4. `gene_source_id`: str gene tftargets id
    5. `source_uri`: str 'https://github.com/slowkow/tftargets'
    6. `s_entrez_id`: str subject entrez
    7. `s_hgnc_id`: str subject hgnc id
    8. `s_symbol`: str subject label
    9. `p_id`:  str property id 'RO:0002434'
    10. `p_label`: str property label 'interacts with'
    11. `o_entrez_id`: str object entrez
    12. `o_hgnc_id`: str object hgnc id
    13. `o_symbol`: str object label
    14. `reference_id`: str pmid
    15. `reference_date`: str yyyy-mm-dd
    :param data: (tred, encode, neph, trrust, msigdb) individual raw networks dictionaries (in this order) as tuple \
    from the load_tf_gene_edges() function
    :param dicts: (symbol2entrez_dict, symbol2hgnc_dict, entrez2hgnc_dict, entrez2symbol_dict) dictionaries tuple \
    from the get_gene_id_normalization_dictionaries() function
    :return: (tftargets, msigdb) edges dataframes tuple
    """

    print('\nThe function "prepare_data_edges()" is running...')
    # raw networks
    tred = data[0]
    encode = data[1]
    neph = data[2]
    trrust = data[3]
    msigdb = data[4]

    # gene_id mapping dictionaries
    symbol2entrez_dict = dicts[0]
    symbol2hgnc_dict = dicts[1]
    entrez2hgnc_dict = dicts[2]
    entrez2symbol_dict = dicts[3]

    ## tftargets network
    # REFERENCES: add ref_uri to trrust statements
    # build trrust dict: {'tf:gene': 'PMID:'}
    st2ref = {}
    #references_path = '/home/nuria/workspace/ngly1-graph/regulation/tftargets/data-raw/TRRUST'
    references_path = './regulation/tftargets/data-raw/TRRUST'
    for line in gzip.open('{}/trrust_rawdata.txt.gz'.format(references_path), 'rt').readlines():
        tf, gene, mor, ref_list = line.strip().split('\t')
        st = tf + ':' + gene
        ref_list = ref_list.split(';') if ';' in ref_list else [ref_list]
        for ref in ref_list:
            add_elem_dictionary(st2ref, st, ref)

    for st in st2ref:
        # st2ref[st] = 'https://www.ncbi.nlm.nih.gov/pubmed/'+','.join(st2ref[st])
        st2ref[st] = 'PMID:' + ';'.join(st2ref[st])

    # build tftargets networks: tred, encode, neph, trrust
    tftargets_path = path + '/tftargets'
    if not os.path.isdir(tftargets_path): os.makedirs(tftargets_path)
    #with open('./tftargets/out/tftargets_edges.csv', 'w') as f:
    with open('{}/tftargets_edges.csv'.format(tftargets_path), 'w') as f:
        f.write(
            "source,dataset,tf_source_id,gene_source_id,source_uri,s_entrez_id,s_hgnc_id,s_symbol,p_id,p_label,\
            o_entrez_id,o_hgnc_id,o_symbol,reference_id,reference_date\n")
        source = "tftargets"
        source_uri = "https://github.com/slowkow/tftargets"
        p_id = "RO:0002434"
        p_label = "interacts with"
        # tred
        dataset = "tred"
        reference_id = "PMID:17202159"
        reference_date = "2007-01-01"
        for tf in tred:
            for gene in tred[tf]:
                f.write(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
                        .format(
                        source,
                        dataset,
                        tf,
                        gene,
                        source_uri,
                        symbol2entrez_dict[tf],
                        symbol2hgnc_dict[tf],
                        tf,
                        p_id,
                        p_label,
                        "NCBIGene:" + str(gene),
                        entrez2hgnc_dict["NCBIGene:" + str(gene)],
                        entrez2symbol_dict["NCBIGene:" + str(gene)],
                        reference_id,
                        reference_date
                    )
                )
        # encode
        dataset = "encode_ENCFF001UUQ"
        reference_id = "ENCODE:ENCFF001UUQ"
        reference_date = "2012-08-28"
        for tf in encode:
            for gene in encode[tf]:
                f.write(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
                        .format(
                        source,
                        dataset,
                        tf,
                        gene,
                        source_uri,
                        symbol2entrez_dict[tf],
                        symbol2hgnc_dict[tf],
                        tf,
                        p_id,
                        p_label,
                        "NCBIGene:" + str(gene),
                        entrez2hgnc_dict["NCBIGene:" + str(gene)],
                        entrez2symbol_dict["NCBIGene:" + str(gene)],
                        reference_id,
                        reference_date
                    )
                )
        # neph
        dataset = "neph2012"
        reference_id = "PMID:22959076"
        reference_date = "2012-09-14"
        for tf in neph:
            for gene in neph[tf]:
                f.write(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
                        .format(
                        source,
                        dataset,
                        tf,
                        gene,
                        source_uri,
                        symbol2entrez_dict[tf],
                        symbol2hgnc_dict[tf],
                        tf,
                        p_id,
                        p_label,
                        "NCBIGene:" + str(gene),
                        entrez2hgnc_dict["NCBIGene:" + str(gene)],
                        entrez2symbol_dict["NCBIGene:" + str(gene)],
                        reference_id,
                        reference_date
                    )
                )
        # trrust
        dataset = "trrust"
        # reference_id = "PMID:26066708"
        reference_date = "NA"  # "2015-06-12"
        for tf in trrust:
            for gene in trrust[tf]:
                reference_id = st2ref[tf + ':' + gene]
                f.write(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
                        .format(
                        source,
                        dataset,
                        tf,
                        gene,
                        source_uri,
                        symbol2entrez_dict[tf],
                        symbol2hgnc_dict[tf],
                        tf,
                        p_id,
                        p_label,
                        symbol2entrez_dict[gene],
                        symbol2hgnc_dict[gene],
                        gene,
                        reference_id,
                        reference_date
                    )
                )
    print('\nThe tftargets edges are saved at: {}/tftargets_edges.csv\n'.format(tftargets_path))

    ## msigdb network
    # REFERENCES: add ref_uri to msigdb statements
    # build msigdb dict: {'tf:gene': 'reference'}
    unknown = list()
    unrecognized = list()
    st2ref = {}
    #references_path = '/home/nuria/workspace/ngly1-graph/regulation/msigdb/data'
    references_path = './regulation/msigdb/data'
    for line in open("{}/c3.tft.v6.1.entrez.gmt".format(references_path)).readlines():
        geneset_name = line.strip().split('\t')[0]
        gs_name_v = geneset_name.split('_')
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
        ref_uri = line.strip().split('\t')[1]
        genelist = line.strip().split('\t')[2:]
        for gene in set(genelist):
            st = tf + ':' + gene
            add_elem_dictionary(st2ref, st, ref_uri)

    # convert reference_uri list to str
    for st in st2ref:
        if len(st2ref[st]) > 1:
            st2ref[st] = '|'.join(st2ref[st])
            # print(st, st2ref[st])
        else:
            # st2ref[st] = st2ref[st]
            for ref in st2ref[st]:
                st2ref[st] = ref

    # build msigdb network
    msigdb_path = path + '/msigdb'
    if not os.path.isdir(msigdb_path): os.makedirs(msigdb_path)
    #with open('./msigdb/out/msigdb_edges.csv', 'w') as f:
    with open('{}/msigdb_edges.csv'.format(msigdb_path), 'w') as f:
        f.write(
            "source,dataset,tf_source_id,gene_source_id,source_uri,s_entrez_id,s_hgnc_id,s_symbol,p_id,p_label,\
            o_entrez_id,o_hgnc_id,o_symbol,reference_id,reference_date\n")
        source = "msigdb"
        source_uri = "http://software.broadinstitute.org/gsea/msigdb"
        p_id = "RO:0002434"
        p_label = "interacts with"
        # c3:tft
        dataset = "c3:tft"
        # reference_id = "NA"
        reference_date = "NA"
        for tf in msigdb:
            for gene in msigdb[tf]:
                reference_id = st2ref[tf + ':' + gene]
                f.write(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
                        .format(
                        source,
                        dataset,
                        tf,
                        gene,
                        source_uri,
                        symbol2entrez_dict[tf],
                        symbol2hgnc_dict[tf],
                        tf,
                        p_id,
                        p_label,
                        "NCBIGene:" + str(gene),
                        entrez2hgnc_dict["NCBIGene:" + str(gene)],
                        entrez2symbol_dict["NCBIGene:" + str(gene)],
                        reference_id,
                        reference_date
                    )
                )
    print('\nThe MSigDB edges are saved at: {}/msigdb_edges.csv\n'.format(msigdb_path))

    tftargets = pd.read_csv('{}/tftargets/tftargets_edges.csv'.format(path), low_memory=False)
    msigdb = pd.read_csv('{}/msigdb/msigdb_edges.csv'.format(path))
    data_edges = (tftargets, msigdb)
    print('\nFinished prepare_data_edges().\n')

    return data_edges


# prepare regulation edges to build the graph
def prepare_regulation_edges(data_edges):
    """
    This function prepares and compiles all individual data edges into regulation edges to build the graph.
    :param data_edges: (tftargets, msigdb) edges dataframes tuple
    :return: network dataframe
    """

    print('\nThe function "prepare_regulation_edges()" is running...')
    # load individual regulation networks
    tftargets = data_edges[0]
    msigdb = data_edges[1]
    # tftargets = pd.read_csv('{}/tftargets/tftargets_edges.csv'.format(path), low_memory=False)
    # msigdb = pd.read_csv('{}/msigdb/msigdb_edges.csv'.format(path))
    #print('tftargets:',tftargets.shape)
    #print('msigdb:',msigdb.shape)

    # concat tftargets and msigdb tg-gene edges
    edges = pd.concat([tftargets, msigdb], ignore_index=True)
    #print('edges:',edges.shape)

    # drop duplicates
    edges.drop_duplicates(inplace=True)
    #print(len(edges))

    # select gene ID to HGNC and rename columns to graph schema
    edges = (edges
             .rename(columns={'s_hgnc_id': 'subject_id',
                              'p_id': 'property_id',
                              'p_label': 'property_label',
                              'o_hgnc_id': 'object_id'
                              })
             )
    print('\nFinished prepare_regulation_edges().\n')

    return edges


# BUILD NETWORK

def build_edges(edges):
    """
    This function builds the edges network with the graph schema.
    :param edges: network dataframe from the prepare_regulation_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    # give graph format
    curie_dct = {
        'ro': 'http://purl.obolibrary.org/obo/',
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',
        'encode': 'https://www.encodeproject.org/search/?searchTerm='
    }

    edges_l = list()
    for i, row in edges.iterrows():
        # property uri: http://purl.obolibrary.org/obo/RO_0002434
        property_uri = 'NA'
        if ':' in row['property_id']:
            property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].replace(':', '_')

        # reference_uri: https://www.ncbi.nlm.nih.gov/pubmed/25416956
        # capture nan or None values, i.e. all possible nulls
        # if (isinstance(row['reference_id'], float) and str(row['reference_id']).lower() == 'nan') or \
        # row['reference_id'] is None:
        #    row['reference_id'] = 'NA'
        if ':' not in str(row['reference_id']):
            reference_uri = row['source_uri']  # row['reference_id']
        elif 'http://www.broadinstitute.org/gsea/msigdb/cards/' in str(row['reference_id']):
            reference_uri = row['reference_id']
        else:
            try:
                reference_uri = curie_dct[row['reference_id'].split(':')[0].lower()] + row['reference_id'].split(':')[
                    1].replace(';', ',')
                # if 'pmid' in row['reference_id'].lower(): print(reference_uri)
            except KeyError:
                reference_uri = row['reference_id']
                print('There is a reference curie with and unrecognized namespace:', row['reference_id'])
        # build list of edges as list of dict, i.e a df, where a dict is an edge
        edge = dict()
        edge['subject_id'] = row['subject_id']
        edge['object_id'] = row['object_id']
        edge['property_id'] = row['property_id']
        edge['property_label'] = row['property_label']
        edge['property_description'] = 'NA'
        edge['property_uri'] = property_uri
        edge['reference_uri'] = reference_uri
        edge['reference_supporting_text'] = 'This edge comes from the {} dataset in "{}" source.'.format(
            row['dataset'].upper(), row['source'])  # 'NA'
        edge['reference_date'] = row['reference_date']
        edges_l.append(edge)

    # save edges file
    pd.DataFrame(edges_l).fillna('NA').to_csv('{}/regulation_edges_v{}.csv'.format(graph,today), index=False)

    # print edges info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe regulation network edges are built and saved at: {}/regulation_edges_v{}.csv\n'.format(graph,today))
    print('\nFinished build_edges().\n')

    return edges_l


def build_nodes(edges):
    """
    This function builds the nodes network with the graph schema.
    :param edges: network dataframe from the prepare_regulation_edges() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodes()" is running...')
    # retrieve node attributes from biothings and build dictionary
    # from biothings we retrieve: name (new attribute for short description), alias (synonyms), summary (description)
    # symbols in this case come from the original source. otherwise are gonna be retrieved from biothings as well.
    # build concept dict: {id:symbol}
    concept_dct = dict()
    for i, row in edges.iterrows():
        # node for subject
        concept_dct[row['subject_id']] = {
            'preflabel': row['s_symbol'],
            'name': None,
            'synonyms': None,
            'description': None
        }
        # node for object
        concept_dct[row['object_id']] = {
            'preflabel': row['o_symbol'],
            'name': None,
            'synonyms': None,
            'description': None
        }
    print('\n* Total number of nodes: {}'.format(len(concept_dct.keys())))

    ### retrieve node attributes from gene symbol using BioThings API
    ## But first, trap genes without symbol, i.e. discontinued entrez
    # some entrez don't exist anymore >> no attributes
    print('\n* Trap genes without gene symbol, i.e. genes with discontinued entrez ID...')
    concepts_no_symbol_set = set()
    for concept in concept_dct:
        if str(concept_dct[concept]['preflabel']) == 'nan':
            concepts_no_symbol_set.add(concept)
    print('* Number of concepts without gene symbol:',len(concepts_no_symbol_set))

    # checked that all are entrez
    print('* Check that all genes without gene symbol are identified by entrez ID...')
    concepts_l = list(concepts_no_symbol_set)
    print('* Number of concepts without gene symbol by namespace: ',pd.DataFrame({'id': concepts_l}).id.apply(lambda x: x.split(':')[0]).value_counts())

    # retrieve symbol for every ncbi gene
    entrez = list()
    for concept in concept_dct:
        if str(concept_dct[concept]['preflabel']) == 'nan':
            entrez.append(concept.split(':')[1])

    # query biothings for retired entrez to symbol
    print('\n* Querying BioThings to map retired Entrez to gene symbols...')
    mg = get_client('gene')
    df = mg.querymany(entrez, scopes='entrezgene,retired', fields='symbol', size=1, as_dataframe=True)

    # build ncbi2symbol dictionary
    e2s_df = df.reset_index().rename(columns={'query': 'entrez'}).copy()
    e2s = dict(zip(e2s_df.entrez, e2s_df.symbol))

    # add symbol for retired entrez
    for concept in concept_dct:
        if str(concept_dct[concept]['preflabel']) == 'nan' and concept.split(':')[0] == 'NCBIGene':
            concept_dct[concept]['preflabel'] = e2s[concept.split(':')[1]]

    ## retrieve node attributes from biothings: alias, name and description
    # input list for api: since by id we have hgnc/entrez or symbol, i am gonna use symbol
    symbols = list()
    symbols2 = list()
    for idx, symbol in concept_dct.items():
        # id = key.split(':')[1] if ':' in key else key
        if str(symbol['preflabel']) != 'nan':
            symbols.append(symbol['preflabel'])
            # withdrawn entrez
        else:
            symbols2.append(symbol['preflabel'])

    # api call
    print('\n* Querying BioThings to retrieve node attributes...')
    symbols = list(set(symbols))
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='alias,name,summary', size=1, as_dataframe=True)
    #print(df.shape)
    #print(len(concept_dct.keys()))

    # dictionaries {id: {name:, alias:, summary:}}
    #print(len(concept_dct))
    for symbol, row in df.iterrows():
        # associate concept to symbol
        for concept in concept_dct:
            # print(concept, symbol, concept_dct[concept]['preflabel'], row)
            # 88 concepts without symbol (16992-88=16904 with symbol)
            if concept_dct[concept]['preflabel'] == symbol:
                # add attributes
                # print(concept, symbol, row['name'], row)
                concept_dct[concept]['name'] = row['name']
                concept_dct[concept]['synonyms'] = row['alias']
                concept_dct[concept]['description'] = row['summary']

    # build a list of nodes as list of dict, i.e a df, where a dict is a node
    nodes_l = list()
    for concept in concept_dct:
        # node for subject
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = 'GENE'
        node['preflabel'] = concept_dct[concept]['preflabel']
        node['name'] = concept_dct[concept]['name']
        node['synonyms'] = '|'.join(list(concept_dct[concept]['synonyms'])) if isinstance(
            concept_dct[concept]['synonyms'], list) else concept_dct[concept]['synonyms']
        node['description'] = concept_dct[concept]['description']
        nodes_l.append(node)

    # save nodes file
    pd.DataFrame(nodes_l).fillna('NA').to_csv('{}/regulation_nodes_v{}.csv'.format(graph,today), index=False)
    #print(len(nodes_l))

    # print nodes info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe regulation network nodes are built and saved at: {}/regulation_nodes_v{}.csv\n'.format(graph,today))
    print('\nFinished build_nodes().\n')

    return nodes_l


# NETWORK MANAGEMENT FUNCTIONS
#TODO: prepare graph functions to get list of nodes, edges..

def _print_nodes(nodes, filename):
    """This function save nodes into a CSV file."""

    # print output file

    #return


if __name__ == '__main__':
    ## check TF name label format
    #data = prepare_msigdb_data()
    #check_msigdb_geneset_name_format(data)
    # prepare msigdb data
    gmt_path = '/home/nuria/workspace/ngly1-graph/regulation/msigdb/data/c3.tft.v6.1.entrez.gmt'
    prepare_msigdb_data(gmt_path)
    # prepare individual networks
    data = load_tf_gene_edges()
    dicts = get_gene_id_normalization_dictionaries(data)
    data_edges = prepare_data_edges(data, dicts)
    # prepare regulation network
    network = prepare_regulation_edges(data_edges)
    # build regulation network
    regulation_edges = build_edges(network)
    regulation_nodes = build_nodes(network)
