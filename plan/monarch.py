# @name: monarch.py
# @description: Module for Monarch network preparation and management
# @version: 2.0
# @date: 22-01-2019
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TODO: Add printings for starting execution of functions
# TODO: Debugg: Why json.decoder.JSONDecodeError?,t
# TODO: add documentation
"""Module for Monarch data"""


import requests
import sys, os
import json
import datetime
import pandas as pd

# VARIABLES
today = datetime.date.today()


# CHECK NETWORK SCHEMA AND NORMALIZE TO GRAPH SCHEMA
# check network schema
# TODO: check functions

## from monarch/add-connections-to-net or regulation/monarch
# TODO: build the following functions according monarch.ipynb and add-connections-to-net.ipynb notebooks

# read monarch edges or connections or get-monarch-connections/monarch_connections.[tsv|csv]

# generate statement file/monarch network edges file with graph schema: add-connections-to-net/monarch_edges_v{}.[tsv|csv]

  # generate static variable: uriPrefixes_dct (url references)

  # generate static variable: dbPrefixes_dct (source/database references)

  # build graph schema network edges

# generate concept file/monarch network nodes file with graph schema: add-connections-to-net/monarch_nodes_v{}.[tsv|csv]

  # build semantic groups dictionary

    # collide concepts in a concept dict

    # list of concept prefixes with dict

    # build conceptPrefix2semantic dict

  # build concept attributes dict

  # build graph schema network nodes



# BUILD NETWORK
# preamble functions in notebook

def hit_monarch_api(node = 'HGNC:17646', rows = 2000):
    """
    This function performs api calls to Monarch to retrieve out and in edges from a query node.
    It retrieves entities plus associations via the BioLink API service.
    It hits two endpoints:
        * association/from - for out edges
        * association/to - for in edges
    It returns out and in edges.

    :param node: node id to query (string). Default: 'NCBIGene:55768'.
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

    :param r_out:
    :param r_in:
    :return:
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
    This function builds edges with an attribute for each node.
    It returns a set of edges, where edges are tuples.

    :param sub_l:
    :param rel_l:
    :param obj_l:
    :param ref_l:
    :param attribute:
    :return:
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
    """This function adds edges in a set."""

    for edge in new:
        keep.add(edge)

    return keep


def keep_nodes(keep, edges, seed):
    """
    This function collects nodes from the edges.
    It filters out: PMID nodes, and nodes related by provenance:
        * None
        * 'dc:source'
        * 'IAO:0000136' or 'is about'
        * 'IAO:0000142' or 'mentions'
    i.e., not biologically related
    It returns a set of nodes.

    :param keep:
    :param edges:
    :param seed:
    :return:
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
    """This function gets the first layer of neighbours and relations."""

    keepNodes = set()
    keepEdges = set()
    seedNodes = set(seed)
    for node in seedNodes:
        try:
            r_out, r_in = hit_monarch_api(node)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            keepEdges = keep_edges(keepEdges, edges)
            keepNodes = keep_nodes(keepNodes, edges, seedNodes)

        except json.decoder.JSONDecodeError:
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keepNodes, keepEdges


def filter_edges(nodes, edges):
    """
    This function filters down edges with both nodes in the nodes set.

    :param nodes:
    :param edges:
    :return:
    """
    nodes = set(nodes)
    keep = set()
    for (start, pred, stop, ref) in edges:
        if {start, stop} <= nodes:
            keep.add((start, pred, stop, ref))

    return keep


def add_attributes(sub_l, rel_l, obj_l, edges):
    """ This function adds attributes to each entity in the edge."""

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

    :param edges:
    :param seed:
    :param nodeType: Introduce node type to keep (string): 'ortho' or 'pheno' (orthologs or phenotypes/diseases, respectively).
    :return:
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
    """This function returns associations retrieved from Monarch among a list of query nodes."""

    keep = set()
    for node in nodes:
        try:
            r_out, r_in = hit_monarch_api(node, 1000)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            filteredEdges = filter_edges(nodes, edges)
            metaFilteredEdges = add_attributes(sub_l, rel_l, obj_l, filteredEdges)
            keep = keep_edges(keep, metaFilteredEdges)

        except json.decoder.JSONDecodeError:
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keep


# From here, lib new functions...
# NETWORK MANAGEMENT FUNCTIONS

def get_neighbours_list(seed_list):
    """This function returns the 1st layer of neighbours from a list of nodes."""

    # print executing function
    print('\nThe function "get_neighbours_list()" is running, please keep calm and have some coffee...')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)

    return list(neighbours)


def get_orthopheno_list(seed_list):
    """This function returns ortho-pheno nodes for a list of nodes."""

    # print executing function
    print('\nThe function "get_orthopheno_list()" is running, please keep calm and have some coffee...')

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

    return list(nodes)


def extract_edges(gene_list):
    """This function returns the Monarch network from a list of query nodes. It retrieves connectivity from Monarch.
    The network variable that returns is a set of tuples (edges).

    Note that i changed the name of this function in v2.0 of the module. This function is named as 'expand_edges' in
    version 1.0 of the module, and used as 'expand_edges' in the graph building notebook
    'graph_hypothesis_1shell_v11012018'.
    """

    # print executing function
    print('\nThe function "expand_edges()" is running, please keep calm and have some coffee...')

    # set network nodes: gene list provided by the user
    nodes = set(gene_list)

    # get connections
    network = get_connections(nodes)

    return network


def print_network(network, filename):
    """This function saves the Monarch expanded network into a CSV file."""

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write(
            'subject_id,subject_label,relation_id,relation_label,object_id,object_label,reference_id_list\n')
        for edge in network:
            edge = ['None' if t is None else t for t in edge]
            f.write('{}\n'.format(','.join(edge)))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))

def print_network2(network, filename):
    """This function saves the Monarch expanded network into a CSV file."""

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
    pd.DataFrame(edges).to_csv('{}/{}_dataframe_v{}.csv'.format(path, filename, today), index=False)

    return print("\nFile '{}/{}_dataframe_v{}.csv' saved.".format(path, filename, today))

def print_nodes(nodes, filename):
    """This function saves Monarch nodes into a CSV file."""

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write('{}\n'.format('\n'.join(list(nodes))))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))

# New functions to debug.

def expand_edges(seed_list):
    """This function returns the Monarch network expanded with the first layer of neighbors from a list of query nodes.
    This function builds monarch 1shell network.
    This function receives a list of nodes and returns a network or edges from Monarch."""

    # get 1shell list of nodes or neighbors
    neighbours, relations = get_neighbours(seed_list)

    # network nodes:  seed + 1shell
    nodes = set(seedList).union(neighbours)

    # get connections for network nodes
    network = get_connections(nodes)

    return network

def orthopheno_expand_edges(seed_list):
    """This function returns the Monarch network expanded with the orthologs and the ortholog associated phenotypes from
     a list of query nodes.
    This function builds monarch 1shell-animal network.
    This function receives a list of nodes and returns a network or edges from Monarch."""

    # get ortholog-phenotypes list
    orthophenoList = get_orthopheno_list(seed_list):

    # network nodes:  seed + neighbors + orthologs + phenotypes
    nodes = set(seed_list).union(set(orthophenoList))

    # get connections for network nodes
    network = get_connections(nodes)

    return network

if __name__ == '__main__':
    #geneList = ['OMIM:615273']  # NGLY1 deficiency
    #network = monarch_expand(geneList)
    # control 1: one gene --> makes sense? what are the results?
    # control 2: one disease --> makes sense? what are the results?
    # control 3: node list --> makes sense? what are the results?
    # what are the json decoding errors
    seedList = ['HGNC:17646','HGNC:633']
    neighbourList = get_neighbours_list(seedList)
    print(len(neighbourList)) # 353
    #orthophenoList = get_orthopheno_list(seedList)
    #print(len(neighbourList), len(orthophenoList))
    #geneList = sum([seedList,neighbourList,orthophenoList], [])
    #print(len(geneList))
    #network = monarch_extract(geneList)
    #print(head(network)
    #print(len(network))
    #print_network(network, 'monarch_orthopeno_network'
