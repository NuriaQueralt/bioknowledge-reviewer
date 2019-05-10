# @name: monarch.py
# @description: Module for Monarch network preparation and management
# @version: 2.0
# @date: 22-01-2019
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TODO: Add printings for starting execution of functions
# TODO: add private functions
"""Module for Monarch data"""

import requests
import sys,os
import json
import datetime
import pandas as pd
from biothings_client import get_client
from tqdm import tqdm


# VARIABLES
# timestamp
today = datetime.date.today()

# path to read modules
#sys.path.insert(0, './graph')

# path to write data
# TODO: write on one path or the other, monarch/ or graph/. name different monarch and monarch_graph networks
# TODO: use the function print to print networks
# TODO: save read data into a 'data/' dir
# case = monarch network
path = os.getcwd() + '/monarch'
if not os.path.isdir(path): os.makedirs(path)

# case = graph connectivity
#graph = os.getcwd() + '/graph'
#if not os.path.isdir(graph): os.makedirs(graph)


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA
# check network schema
# TODO: check functions and document

def read_connections(filename):
    """
    This function reads monarch_connections CSV file.
    :param filename: complete path to the monarch connections csv file string
    :return: monarch edges dataframe
    """

    # read monarch_connections or get-monarch-connections/monarch_connections.[tsv|csv]
    # monarch network
    # csv_path = ''
    # monarch graph network
    #csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/monarch_connections_regulation_graph.tsv'
    # csv_path = '/home/nuria/workspace/ngly1-graph/monarch/1shell-animal/get-monarch-connections/monarch_connections.tsv'
    path = os.getcwd() + '/monarch'
    csv_path = path + '/' + filename
    #TODO: do not parse correctly labels due to commas within str labels
    network_df = pd.read_csv('{}'.format(csv_path))
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))

    return network_df


# retrieve subnetwork from Monarch knowledge graph

def hit_monarch_api(node = 'HGNC:17646', rows = 2000):
    """
    This function performs api calls to Monarch to retrieve out and in edges from a query node.
    It retrieves entities plus associations via the BioLink API service.
    It hits two endpoints:
        * association/from - for out edges
        * association/to - for in edges
    It returns out and in edges.

    :param node: node id to query (string). Default: 'HGNC:17646'.
    :param rows: the maximum number of results to return (integer). Default: 2000.
    :return: two api response objects: 'out' and 'in' response objects, in this order.
    """
    
    # api address
    biolink = 'https://api.monarchinitiative.org/api/association'
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': rows}
    # out edges: from/
    r_out = requests.get('{}/from/{}'.format(biolink,node),params=parameters)

    # in edges: to/
    r_in = requests.get('{}/to/{}'.format(biolink,node),params=parameters)

    return r_out, r_in


def get_edges_objects(r_out, r_in):
    """
    This function prepares the api object responses from Monarch.
    It returns four lists, one for subjects, relations, objects, and references.
    Subjects, relations and objects are lists of dictionaries, where each dictionary is a node.
    References list lists strings, where each string is a chain of references for each edge.

    :param r_out: BioLink API 'out' response object
    :param r_in: BioLink API 'in' response object
    :return: subjects, relations, objects and references lists (in this order)
    """

    # variables
    sub_l = list()
    rel_l = list()
    obj_l = list()
    ref_l = list()

    # compose list of dictionaries
    for associations in [r_out.json()['associations'], r_in.json()['associations']]:
        for association in associations:
            pub_l = list()
            sub_l.append(association['subject'])
            rel_l.append(association['relation'])
            obj_l.append(association['object'])
            # add references to each association as a list of strings
            if association['publications']:
                for publication in association['publications']:
                    pub_l.append(publication['id'])
            else:
                pub_l.append('NA')
            ref_l.append('|'.join(pub_l))

    return sub_l, rel_l, obj_l, ref_l


def get_edges(sub_l, rel_l, obj_l, ref_l, attribute='id'):
    """
    This function builds edges using a user-specified attribute for each node.
    It returns a set of edges, where edges are tuples.

    :param sub_l: subjects (objects) list from the get_edges_objects() function
    :param rel_l: relations (objects) list from the get_edges_objects() function
    :param obj_l: objects (objects) list from the get_edges_objects() function
    :param ref_l: references (strings) list from the get_edges_objects() function
    :param attribute: object attribute, default 'id'
    :return: edges (as tuples) set
    """
    edges = set()
    # compose tuple
    for i in range(len(sub_l)):
        sub = sub_l[i][attribute]
        rel = rel_l[i][attribute]
        obj = obj_l[i][attribute]
        ref = ref_l[i]
        edges.add((sub, rel, obj, ref))

    return edges


def keep_edges(keep, new):
    """
    This function adds edges from a new set to a keep set.
    :param keep: edges set
    :param new: edges set
    :return: updated edges set
    """

    for edge in new:
        keep.add(edge)

    return keep


def keep_nodes(keep, edges, seed):
    """
    This function collects nodes from the edges that are not in the nodes query list to keep nodes set.
    It filters out: PMID nodes, and nodes related by provenance:
        * None
        * 'dc:source'
        * 'IAO:0000136' or 'is about'
        * 'IAO:0000142' or 'mentions'
    i.e., not biologically related
    It returns a set of nodes.

    :param keep: nodes set
    :param edges: edges set
    :param seed: query nodes list
    :return: updated nodes set
    """

    for (sub, rel, obj, ref) in edges:
        # if ':.' in (sub or obj):
        #    continue
        # if 'Coriell' in (sub or obj):
        #    continue
        # if 'MMRRC' in (sub or obj):
        #    continue
        # if 'MONARCH' in (sub or obj):
        #    continue
        if 'PMID' in sub or 'PMID' in obj:
            continue
        if rel == None:
            rel = 'None'
        if 'dc:source' in rel:
            continue
        if 'IAO:0000136' in rel:  # is about
            continue
        if 'IAO:0000142' in rel:  # mentions
            continue
        if sub not in seed:
            keep.add(sub)
        if obj not in seed:
            keep.add(obj)

    return keep


def get_neighbours(seed):
    """
    This function gets the first layer of neighbours and relations.
    :param seed: query nodes list
    :return: nodes set, edges set (in this order)
    """

    keepNodes = set()
    keepEdges = set()
    seedNodes = set(seed)
    for node in tqdm(seedNodes):
        try:
            r_out, r_in = hit_monarch_api(node)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            keepEdges = keep_edges(keepEdges, edges)
            keepNodes = keep_nodes(keepNodes, edges, seedNodes)

        except (ValueError, KeyError):
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keepNodes, keepEdges


def filter_edges(nodes, edges):
    """
    This function filters down edges with both nodes in a nodes list.

    :param nodes: nodes list
    :param edges: edges set
    :return: filtered edges set
    """
    nodes = set(nodes)
    keep = set()
    for (start, pred, stop, ref) in edges:
        if {start, stop} <= nodes:
            keep.add((start, pred, stop, ref))

    return keep


def add_attributes(sub_l, rel_l, obj_l, edges):
    """
    This function adds 'label' attribute to each entity in the edge.
    :param sub_l: subjects (object) list
    :param rel_l: relations (object) list
    :param obj_l: objects (object) list
    :param edges: edges set
    :return: metaedges set
    """

    metaedges = set()
    for (sub_id, rel_id, obj_id, refs) in edges:
        for i in range(len(sub_l)):
            if sub_l[i]['id'] == sub_id and rel_l[i]['id'] == rel_id and obj_l[i]['id'] == obj_id:
                metaedges.add((sub_l[i]['id'],
                               sub_l[i]['label'],
                               rel_l[i]['id'],
                               rel_l[i]['label'],
                               obj_l[i]['id'],
                               obj_l[i]['label'],
                               refs)
                              )
                break
    return metaedges


def keep_node_type(edges, seed, nodeType='ortho'):
    """
    This function keeps specific node types for objects in edges.

    :param edges: edges set
    :param seed: the query nodes list
    :param nodeType: Introduce node type to keep (string): 'ortho' for orthologs or 'pheno' \
    for phenotypes/diseases, default is 'ortho'
    :return: nodes set
    """

    propertyList = ['RO:HOM0000017', 'RO:HOM0000020']
    if nodeType == 'pheno':
        propertyList = ['RO:0002200', 'RO:0002607', 'RO:0002326', 'GENO:0000840']

    keep = set()
    for (sub, rel, obj, ref) in edges:
        if rel == None:
            continue
        if rel in propertyList:
            if sub not in seed:
                keep.add(sub)
            if obj not in seed:
                keep.add(obj)

    return keep


def get_connections(nodes):
    """
    This function returns associations retrieved from Monarch among a list of query nodes."
    :param nodes: the query nodes list
    :return: edges set
    """""

    keep = set()
    for node in tqdm(nodes):
        try:
            r_out, r_in = hit_monarch_api(node, 1000)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            filteredEdges = filter_edges(nodes, edges)
            metaFilteredEdges = add_attributes(sub_l, rel_l, obj_l, filteredEdges)
            keep = keep_edges(keep, metaFilteredEdges)

        except (ValueError, KeyError):
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keep


# NETWORK MANAGEMENT FUNCTIONS

def get_neighbours_list(seed_list):
    """
    This function returns the first explicit layer of neighbours from a list of query nodes.
    :param seed_list: biomedical entities list, where each entity is the identifier string like 'HGNC:17646'
    :return: neighbours list
    """

    # print executing function
    print('\nThe function "get_neighbours_list()" is running. Its runtime may take some minutes. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)
    print('\nFinished get_neighbours_list().\n')

    return list(neighbours)


def get_orthopheno_list(seed_list):
    """
    This function returns orthologs-phenotypes nodes in ortho-pheno relationships for a list of query genes.
    :param seed_list: gene list, where each gene is the identifier string like 'HGNC:17646'
    :return: orthopheno list
    """

    # print executing function
    print('\nThe function "get_orthopheno_list()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)

    # keep orthologs in the first layer
    orthologs = keep_node_type(relations, seed_list)

    # get second layer from orthologs
    neighbours, relations = get_neighbours(orthologs)

    # keep phenotypes in the second layer
    phenotypes = keep_node_type(relations, orthologs, 'pheno')

    nodes = set()
    nodes.update(orthologs, phenotypes)
    print('\nFinished get_orthopheno_list().\n')

    return list(nodes)


def extract_edges(gene_list):
    """
    This function returns the Monarch network from a list of query nodes. It retrieves connectivity from Monarch, i.e. \
    edges from Monarch between query nodes.
    :param gene_list: gene list
    :return: edges (as tuples) set
    """

    # print executing function
    print('\nThe function "extract_edges()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the edges retrieved '
          'and you should start over the execution of this function.')

    # set network nodes: gene list provided by the user
    nodes = set(gene_list)

    # get connections
    network = get_connections(nodes)
    print('\nFinished extract_edges(). To save the retrieved Monarch edges use the function "print_network()".\n')

    return network


def _print_network2(network, filename):
    """This function saves the Monarch expanded network into a CSV file. this function save connections file format into
    get-monarch-connections/"""
    # TODO: this function only writes down the network connections format from monarch,
    #  originally stored into get-monarch-connections/

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write(
            'subject_id,subject_label,relation_id,relation_label,object_id,object_label,reference_id_list\n')
        for edge in network:
            edge = ['None' if t is None else '"{}"'.format(str(t)) for t in edge]
            f.write('{}\n'.format(','.join(edge)))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def print_network(network, filename):
    """
    This function saves the Monarch network into a CSV file. It only prints connections file format only.
    :param network: monarch edges dataframe
    :param filename: file name without extension string, e.g. 'monarch_connections'
    :return: None object
    """

    # transform set of tuples to list of dictionaries
    edges = list()
    for tuple in network:
        row = dict()
        row['subject_id'] = tuple[0]
        row['subject_label'] = tuple[1]
        row['relation_id'] = tuple[2]
        row['relation_label'] = tuple[3]
        row['object_id'] = tuple[4]
        row['object_label'] = tuple[5]
        row['reference_id_list'] = tuple[6]
        edges.append(row)

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    pd.DataFrame(edges).fillna('None').to_csv('{}/{}_v{}.csv'.format(path, filename, today), index=False)

    return print("\nSaving Monarch edges at: '{}/{}_v{}.csv'...\n".format(path, filename, today))


def print_nodes(nodes, filename):
    """
    This function saves Monarch nodes into a CSV file.
    :param nodes: nodes list
    :param filename: file name without path and extension, e.g. 'monarch_nodes'
    :return: None object
    """

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write('{}\n'.format('\n'.join(list(nodes))))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def expand_edges(seed_list):
    """
    This function returns the Monarch network expanded with the first layer of neighbors from a list of query nodes.
    This function builds monarch 1shell network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    # get 1shell list of nodes or neighbors
    neighbours, relations = get_neighbours(seed_list)

    # network nodes:  seed + 1shell
    nodes = set(seed_list).union(neighbours)

    # get connections for network nodes
    network = get_connections(nodes)

    return network


def orthopheno_expand_edges(seed_list):
    """
    This function returns the Monarch network expanded with the orthologs and the ortholog associated phenotypes from
     a list of query nodes.
    This function builds monarch 1shell-animal network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    # get ortholog-phenotypes list
    orthophenoList = get_orthopheno_list(seed_list)

    # network nodes:  seed + neighbors + orthologs + phenotypes
    nodes = set(seed_list).union(set(orthophenoList))

    # get connections for network nodes
    network = get_connections(nodes)

    return network


# BUILD NETWORK

def build_edges(edges_df):
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    ## variables
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['relation_id'] = tuple[2]
            record['relation_label'] = tuple[3]
            record['object_id'] = tuple[4]
            record['object_label'] = tuple[5]
            record['reference_id_list'] = tuple[6]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)


    # generate static variable: uriPrefixes_dct (url references)
    uriPrefixes_dct = {
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  # 'http://identifiers.org/pubmed/',
        'react': 'http://reactome.org/content/detail/',  # 'http://identifiers.org/reactome/',
        'zfin': 'http://zfin.org/',
        'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  # 'http://identifiers.org/go.ref/GO_REF:',
        'mgi': 'http://www.informatics.jax.org/accession/MGI:',  # 'http://identifiers.org/mgi/MGI:'
        'flybase': 'http://flybase.org/reports/',
        'wormbase': 'http://www.wormbase.org/resources/paper/',
        'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
        'isbn-10': 'ISBN-10:',
        'isbn-13': 'ISBN-13:',
        #'isbn-10': 'https://www.wikidata.org/wiki/Special:BookSources/',
        #'isbn-13': 'https://www.wikidata.org/wiki/Special:BookSources/'
        'mondo': 'http://purl.obolibrary.org/obo/MONDO_',  # http://purl.obolibrary.org/obo/MONDO_0009026
        'rgd': 'https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=', \
        # https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=1600115
        'omim': 'http://omim.org/entry/',  # http://omim.org/entry/61527
        'sgd_ref': 'https://db.yeastgenome.org/reference/',  # https://db.yeastgenome.org/reference/S000124036
        'genereviews': 'https://www.ncbi.nlm.nih.gov/books/',  # https://www.ncbi.nlm.nih.gov/books/NBK1526/
        'omia': 'http://omia.angis.org.au/',  # http://omia.angis.org.au/000214/9913
        'hgnc': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:', \
        # https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:7132
        'orpha': 'Go to ORPHANET web site (https://www.orpha.net/) and in the search field introduce the Orpha number: '
                 'ORPHA:' # no entry in monarch for that edge
    }

    # generate static variable: dbPrefixes_dct (source/database references)
    dbPrefixes_dct = {
        'na': 'NA',
        'nan': 'NA',
        'mgi': 'http://www.informatics.jax.org/',
        'fb': 'http://flybase.org/',
        'rgd': 'http://rgd.mcw.edu/',
        'zfin': 'http://zfin.org/',
        'sgd': 'https://www.yeastgenome.org/',
        'hgnc': 'https://www.genenames.org/',
        'xenbase': 'http://www.xenbase.org/'
    }

    # provenance variables
    ref_text = 'This edge comes from the Monarch Knowledge Graph 2019.'  # 'NA'
    ref_date = 'NA'


    ## build graph schema network edges data structure and save edges file
    #csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/monarch_connections_regulation_graph.tsv'
    #csv_path = '/home/nuria/workspace/ngly1-graph/monarch/1shell-animal/get-monarch-connections/monarch_connections.tsv'
    # edges_df = pd.read_table('./get-monarch-connections/monarch_connections.tsv')
    # prepare dataframe = [ {} ... {} ], where every row = {} = concept
    edges_l = list()
    for edge in edges_df.itertuples():
        # edge or row is a tuple (named and ordered attributes)
        #print("edge_tuple:", edge)
        # edge.reference_id_list >> can be 1) np.nan (float type) or 2) str without "|" 3) str with "|"
        ref_s = str(edge.reference_id_list)

        # prepare reference_uri_list attribute
        ref_uri_l = list()
        # expand to uri or NA
        pmid_l = list()
        # reference_id list iteration
        for ref in ref_s.strip().split('|'):
            # NA or database
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("Warning:")
                    print('Detected a new reference database in Monarch not yet implemented in this module. '
                          'The new source should be added to the dictionary of databases.'
                          'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                    print("In the build_edges() method, update 'dbPrefixes_dct' dictionary with '{}'".format(ref))
                    print('The edge that includes this new reference database is {}'.format(edge))
                    print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                          'until the dictionary is updated.')
                ref_uri_l.append(ref_uri)
            # publication uri: pubmed_id or url
            else:
                pref, uriId = ref.split(':')
                # separate pmid from non pmid and detect:
                # pubmed_id
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                # url
                elif ref.lower().startswith('http'):
                    ref_uri_l.append(ref)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("Warning:")
                        print('Detected a new reference source in Monarch not yet implemented in this module. '
                              'The new source should be added to the dictionary of sources.'
                              'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                        print("In the build_edges() method, update 'uriPrefixes_dct' dictionary with '{}'".format(pref))
                        print('The edge that includes this new reference source is {}'.format(edge))
                        print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                              'until the dictionary is updated.')
                    ref_uri_l.append(ref_uri)
        # create multi-term pubmed url
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)

        # prepare edge attributes: sub_id, obj_id, rel_id, rel_label, rel_def, rel_iri
        sub_id = 'NA' if edge.subject_id is None or str(edge.subject_id) == 'nan' else edge.subject_id
        rel_id = 'NA' if edge.relation_id is None or str(edge.relation_id) == 'nan' else edge.relation_id
        obj_id = 'NA' if edge.object_id is None or str(edge.object_id) == 'nan' else edge.object_id
        rel_label = 'NA' if edge.relation_label is None or str(edge.relation_label) == 'nan' else edge.relation_label
        rel_def = 'NA'
        if ':' in rel_id:
            rel_iri = 'http://purl.obolibrary.org/obo/' + rel_id.replace(':', '_')
        else:
            rel_iri = rel_id

        # build the data structure = list of edges as list of dict, where a dict is an edge
        edge = dict()
        edge['subject_id'] = sub_id
        edge['object_id'] = obj_id
        edge['property_id'] = rel_id
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] = rel_iri
        edge['reference_uri'] = ref_uri_list
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edges_l.append(edge)

    # save edges file
    #TODO: abstract this function
    df = pd.DataFrame(edges_l)
    print('df',df.shape)
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', \
             'property_label', 'property_description', 'property_uri']]
    df.fillna('NA').to_csv('{}/monarch_edges_v{}.csv'.format(path,today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe Monarch network edges are built and saved at: {}/monarch_edges_v{}.csv\n'.format(path,today))
    print('\nFinished build_edges().\n')

    return edges_l


def build_nodes(edges_df):
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodes()" is running...')
    # build semantic groups dictionary
    # collide concepts in a concept dict
    concept_dct = dict()
    # read edges from input file
    #header = 1
    ##for row in open('../monarch/1shell-animal-hgnc/get-monarch-connections/monarch_connections.tsv').readlines():
    ##csv_path = '/home/nuria/workspace/ngly1-graph/regulation/graph/monarch_connections_regulation_graph.tsv'
    #csv_path = '/home/nuria/workspace/ngly1-graph/monarch/1shell-animal/get-monarch-connections/monarch_connections.tsv'
    #for row in open('{}'.format(csv_path)).readlines():
    #    if header:
    #        header = 0
    #        continue
    # read edges from variable
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['relation_id'] = tuple[2]
            record['relation_label'] = tuple[3]
            record['object_id'] = tuple[4]
            record['object_label'] = tuple[5]
            record['reference_id_list'] = tuple[6]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        #fields = row.strip('\n').split('\t')
        #fields = list(edge_tuple)
        sid = edge.subject_id #fields[0]
        oid = edge.object_id #fields[4]
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))

    # list of concept prefixes with dict
    conceptPrefix_dct = dict()
    for concept in concept_dct:
        conceptPrefix_dct[concept.split(':')[0]] = 1
    print('Number of nodes CURIEs: {}'.format(len(conceptPrefix_dct.keys())))
    print('List of nodes CURIEs: {}'.format(conceptPrefix_dct.keys()))

    # build conceptPrefix2semantic dict
    conceptPrefix2semantic_dct = dict()
    for prefix in conceptPrefix_dct:
        prefix = prefix.lower()
        if 'variant' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'VARI'
        elif 'phenotype' in prefix or 'mondo' in prefix or 'omim' in prefix or 'doid' in prefix or 'mesh' in prefix or 'hp' in prefix or 'mp' in prefix or 'fbcv' in prefix or 'fbbt' in prefix or 'zp' in prefix or 'apo' in prefix or 'trait' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'DISO'
        elif 'gene' in prefix or 'hgnc' in prefix or 'ensembl' in prefix or 'mgi' in prefix or 'flybase' in prefix or 'wormbase' in prefix or 'xenbase' in prefix or 'zfin' in prefix or 'rgd' in prefix or 'sgd' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENE'
        elif 'react' in prefix or 'kegg-path' in prefix or 'go' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'PHYS'
        elif 'uberon' in prefix or 'cl' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'ANAT'
        elif 'geno' in prefix or 'coriell' in prefix or 'monarch' in prefix or 'mmrrc' in prefix or '' in prefix \
                or 'bnode' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENO'
        else:
            conceptPrefix2semantic_dct[prefix] = 'CONC'

    # build concept attributes dict: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()
    # read edges from file
    #header = 1
    ##for row in open('../monarch/1shell-animal-hgnc/get-monarch-connections/monarch_connections.tsv').readlines():
    #for row in open('{}'.format(csv_path)).readlines():
    #    if header:
    #        header = 0
    #        continue
    #    fields = row.strip('\n').split('\t')
    # read edges from variable
    for edge in edges_df.itertuples():
        #fields = list(edge_tuple)
        # id: integration of sub and obj IDs in a unique data structure
        sid = edge.subject_id #fields[0]
        slab = edge.subject_label #fields[1]
        oid = edge.object_id #fields[4]
        olab = edge.object_label #fields[5]
        # build the concept data structure
        concept_dct[sid] = {'preflabel': slab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(sid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(oid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'description': 'NA'}

    # build graph schema network nodes data structure and save nodes file
    # biothings: annotate name,synonyms,description to genes
    print('\nAdding BioThings annotation: gene name, synonyms, description...')
    # input: (preflabel) symbol,alias
    symbols = list()
    for concept in concept_dct:
        if isinstance(concept_dct[concept]['semantic_groups'], list):
            for label in concept_dct[concept]['semantic_groups']:
                if 'GENE' in label:
                    symbols.append(concept_dct[concept]['preflabel'])
        else:
            if 'GENE' in concept_dct[concept]['semantic_groups']:
                symbols.append(concept_dct[concept]['preflabel'])
    print('symbols:', len(symbols))

    # query biothings
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name,alias,summary', size=1, as_dataframe=True)

    # dictionary: {symbol:name}
    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
    ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
    monarch_s2n = dict(zip(ids.symbol, ids.name))
    monarch_s2s = dict(zip(ids.symbol, ids.synonyms))
    monarch_s2d = dict(zip(ids.symbol, ids.description))

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        preflabel = concept_dct[concept]['preflabel']
        concept_dct[concept]['synonyms'] = monarch_s2s[preflabel] if preflabel in monarch_s2s.keys() else 'NA'
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['preflabel'] = preflabel
        node['name'] = monarch_s2n[preflabel] if preflabel in monarch_s2n.keys() else preflabel
        node['synonyms'] = '|'.join(list(concept_dct[concept]['synonyms'])) if isinstance(
            concept_dct[concept]['synonyms'], list) else concept_dct[concept]['synonyms']
        node['description'] = monarch_s2d[preflabel] if preflabel in monarch_s2d.keys() else 'NA'
        nodes_l.append(node)

    # save nodes file
    #TODO: abstract the print function
    df = pd.DataFrame(nodes_l)
    df = df[['id', 'semantic_groups', 'preflabel', 'synonyms', 'description', 'name']]
    #TODO: check why i am saving as csv but naming the file tsv
    df.fillna('NA').to_csv('{}/monarch_nodes_v{}.csv'.format(path,today), index=False)

    # print info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe Monarch network nodes are built and saved at: {}/monarch_nodes_v{}.csv\n'.format(path,today))
    print('\nFinished build_nodes().\n')

    return nodes_l


if __name__ == '__main__':
    #geneList = ['OMIM:615273']  # NGLY1 deficiency
    #network = monarch_expand(geneList)
    # build monarch network
    #seedList = ['HGNC:17646','HGNC:633']
    #neighbourList = get_neighbours_list(seedList)
    #print(len(neighbourList)) # 353
    #orthophenoList = get_orthopheno_list(seedList)
    #print(len(neighbourList), len(orthophenoList))
    #geneList = sum([seedList,neighbourList,orthophenoList], [])
    #print(len(geneList))
    #network = monarch_extract(geneList)
    #print(head(network)
    #print(len(network))
    #print_network(network, 'monarch_orthopeno_network'
    #read_connections()
    #build_edges()
    #build_nodes()
    #external function print
    # TODO: external prints (use print_network2() as the default)
    #nodes = build_nodes()
    #print_network(nodes,'nodes_print1_prova')
    # build monarch graph from monarch connections network
    filepath = '/home/nuria/workspace/graph-hypothesis-generation-lib/plan/monarch/monarch_v2019-03-01.csv'
    monarch_connections = read_connections(filepath) # OR monarch_network = read_connections()
    print('### len of monarch_connections input:',len(monarch_connections))
    monarch_edges = build_edges(monarch_connections)
    print('### len of monarch_edges output:',len(monarch_edges))
    monarch_nodes = build_nodes(monarch_connections)
    print('### len of monarch_nodes output:',len(monarch_nodes))

    # print network
    # file_name = 'monarch_prova2_dataframe_v2019-03-01.csv'
    # monarch_connections = read_connections(file_name)
    # build_edges(monarch_connections)
    # build_nodes(monarch_connections)
    # path = '/home/nuria/workspace/graph-hypothesis-generation-lib/plan/monarch'
    # pd.read_table('{}/monarch_v2019-03-01.csv'.format(path), sep=',')