# @name: graph.py
# @version: 1.0
# @author: Núria Queralt Rosinach
# @date: 16-02-2018
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

# TO DO: modularize from:
# http://localhost:8888/notebooks/workspace/ngly1-graph/curation/kylo/neo4j/networks/concatenate_network_files.ipynb

"""Module for functions to build the graph"""

import pandas as pd
import sys, os
import datetime

# VARIABLES
today = datetime.date.today()


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


################ CURATION SECTION

def create_edges():
    """This function creates edge data."""

    return edges_df

def create_nodes():
    """This function creates edge data."""

    return nodes_df

################ MONARCH SECTION
# create monarch edges and nodes data structures
# Functionality from http://localhost:8889/notebooks/workspace/ngly1-graph/ngly1-graph-v20180130/build-graph/add_connections_to_net.ipynb

# edges_df = pd.read_table('./get-monarch-connections/monarch_connections.tsv')

def monarch_edges(edges_df):
    """This function creates edge data."""

    # add attribute/columns: 'rel_term_id', 'rel_term_label', rel_term_iri to file
    # authoritative url
    uriPrefixes_dct = {
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  # 'http://identifiers.org/pubmed/',
        'react': 'http://reactome.org/content/detail/',  # 'http://identifiers.org/reactome/',
        'zfin': 'http://zfin.org/',
        'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  # 'http://identifiers.org/go.ref/GO_REF:',
        'mgi': 'http://www.informatics.jax.org/accession/MGI:',  # 'http://identifiers.org/mgi/MGI:'
        'flybase': 'http://flybase.org/reports/',
        'wormbase': 'http://www.wormbase.org/resources/paper/',
        'isbn-13': 'ISBN-13:',
        'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
        'isbn-10': 'ISBN-10:'
    }
    # source/database
    dbPrefixes_dct = {
        'na': 'NA',
        'mgi': 'http://www.informatics.jax.org/',
        'fb': 'http://flybase.org/',
        'rgd': 'http://rgd.mcw.edu/',
        'zfin': 'http://zfin.org/',
        'sgd': 'https://www.yeastgenome.org/',
        'hgnc': 'https://www.genenames.org/'
    }
    # prepare dataframe = [{} ...{}], where every row = {} = concept
    ref_text = 'NA'
    ref_date = 'NA'
    ### temporary
    path = os.getcwd() + '/graph'
    if not os.path.isdir(path): os.makedirs(path)
    ###
    # prepare dataframe = [{} ...{}], where every row = {} = edge
    edges_l = list()
    #for edge in edge_df:
    for edge_tuple in edges_df:
        edge = list(edge_tuple)
        ref_uri_l = list()
        # expand to uri or NA
        pmid_l = list()
        for ref in edge[-1].strip().split('|'):
            # NA or database
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("In monarch_edges() method, update 'dbPrefixes_dct' variable with '{}'".format(ref))
                    print(edge)
                ref_uri_l.append(ref_uri)
            # publications
            else:
                pref, uriId = ref.split(':')
                # separate pmid from non pmid
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("In monarch_edges() method, update 'uriPrefixes_dct' variable with '{}'".format(pref))
                        print(edge)
                    ref_uri_l.append(ref_uri)
        # create multi-term pubmed url
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)
        # write the associations + list of references as uri or NA
        sub_id = 'NA' if edge[0] is None else edge[0]
        rel_id = 'NA' if edge[2] is None else edge[2]
        obj_id = 'NA' if edge[4] is None else edge[4]
        rel_label = 'NA' if edge[3] is None else edge[3]
        rel_def = 'NA' 
        if ':' in rel_id:
            rel_iri = 'http://purl.obolibrary.org/obo/' + rel_id.replace(':', '_')
        else:
            rel_iri = rel_id
        #print(sub_id,rel_id,obj_id, rel_iri)
        # define concept rows as dict for the data structure
        edge = dict()
        edge['subject_id'] = sub_id
        edge['property_id'] = rel_id
        edge['object_id'] = obj_id
        edge['reference_uri'] = ref_uri_list
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] = rel_iri
        edges_l.append(edge)


    return edges_l


def monarch_nodes(nodes_df):
    """This function creates node data."""

    # semantic groups dictionary
    # collide concepts
    concept_dct = dict()
    for node_tuple in nodes_df:
        fields = list(node_tuple)
        sid = fields[0]
        oid = fields[4]
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    #print(len(concept_dct.keys()))

    # list of concept prefixes
    conceptPrefix_dct = dict()
    for concept in concept_dct:
        conceptPrefix_dct[concept.split(':')[0]] = 1
    #print(conceptPrefix_dct.keys())

    # build conceptPrefix2semantic dict
    conceptPrefix2semantic_dct = dict()
    for prefix in conceptPrefix_dct:
        prefix = prefix.lower()
        if 'variant' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'VARI'
        elif 'phenotype' in prefix or 'mondo' in prefix or 'omim' in prefix or 'doid' in prefix or 'mesh' in prefix or 'hp' in prefix or 'mp' in prefix or 'fbcv' in prefix or 'fbbt' in prefix or 'apo' in prefix or 'aqtltrait' in prefix or 'zp' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'DISO'
        elif 'gene' in prefix or 'hgnc' in prefix or 'ensembl' in prefix or 'mgi' in prefix or 'flybase' in prefix or 'wormbase' in prefix or 'xenbase' in prefix or 'zfin' in prefix or 'sgd' in prefix or 'rgd' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENE'
        elif 'react' in prefix or 'kegg-path' in prefix or 'go' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'PHYS'
        elif 'uberon' in prefix or 'cl' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'ANAT'
        elif 'geno' in prefix or 'coriell' in prefix or 'monarch' in prefix or 'mmrrc' in prefix or 'person' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENO'
        else:
            conceptPrefix2semantic_dct[prefix] = 'CONC'

    # concept attribute dictionaries: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()
    for node_tuple in nodes_df:
        fields = list(node_tuple)
        # id: integration of sub and obj IDs in a unique data structure
        sid = fields[0]
        slab = fields[1]
        oid = fields[4]
        olab = fields[5]
        # build the concept data structure
        concept_dct[sid] = {'preflabel': slab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(sid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'definition': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(oid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'definition': 'NA'}

    # prepare dataframe = [{} ...{}], where every row = {} = concept
    nodes_l = list()
    for concept in concept_dct:
        # semantic_groups
        semantic = concept_dct.get(concept).get('semantic_groups')
        # preflabel
        preflabel = concept_dct.get(concept).get('preflabel')
        # synonyms
        synonyms = concept_dct.get(concept).get('synonyms')
        # definition
        definition = concept_dct.get(concept).get('definition')
        # define concept rows as dict for the data structure
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = semantic
        node['preflabel'] = preflabel
        node['synonyms'] = synonyms
        node['description'] = definition
        nodes_l.append(node)

    return nodes_l


################ INTEGRATION
# if there is integration with curation, i.e. with proteins, add g2p edges
# Functionality from http://localhost:8889/notebooks/workspace/ngly1-graph/ngly1-graph-v20180130/build-graph/add_connections_g2p.ipynb


################ CONCAT

################ USER FUNCTIONS

def build(network_list):
    """This function formats, integrates and concats a list of networks passed."""

    # print of starting the process
    print('\nThe function "build()" is running, please keep calm and have some coffee...')

    # build the graph
    graph = dict()
    for df in network_list:
        graph['edges'] = monarch_edges(df)
        graph['nodes'] = monarch_nodes(df)

    return graph

def edges(graph):
    """This function retrieves graph edges."""

    return graph.get('edges')


def nodes(graph):
    """This function retrieves graph edges."""

    return graph.get('nodes')


if __name__ == '__main__':
    df = pd.read_table('./get-monarch-connections/monarch_connections.tsv')
    #monarch_edges(df)
    #nodes = monarch_nodes(df)
    graph = build([df])
    edges = edges(graph)
    nodes = nodes(graph)
    #print_graph(nodes,'monarch_nodes')
    #print_graph(edges,'monarch_edges')
