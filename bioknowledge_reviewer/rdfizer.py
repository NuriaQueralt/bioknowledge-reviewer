# -*- coding: utf-8 -*-
# @name: rdfizer.py
# @description: Module to turn the Neo4j graph into RDF graph. Adapted from Carmen Reep's 'rdf2vec.py' module.
# @version: 1.0
# @date: 03-05-2023
# @author: Núria Queralt Rosinach
# @email: n.queralt_rosinach@lumc.nl
"""Module to turn the Neo4j graph into RDF graph"""

import os
import datetime
import csv
import pandas as pd
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS

# VARIABLES
today = datetime.date.today()
# based on the biolink datamodel [https://github.com/biolink/biolink-model][https://bioportal.bioontology.org/ontologies/BIOLINK/]
uriType_dict = {
    'anatomical entity': 'http://semanticscience.org/resource/SIO_001262',
    'disease': 'http://semanticscience.org/resource/SIO_010299',
    'molecular function': 'http://purl.obolibrary.org/obo/GO:0003674',
    'biological process': 'http://purl.obolibrary.org/obo/GO:0008150',
    'genotype': 'http://semanticscience.org/resource/SIO_001079',
    'gene': 'http://semanticscience.org/resource/SIO_010035',
    'phenotype': 'http://semanticscience.org/resource/SIO_010056',
    'variant': 'http://semanticscience.org/resource/SIO_001381',
    'pathway': 'http://semanticscience.org/resource/SIO_001107',
    'cellular component': 'http://purl.obolibrary.org/obo/GO:0005575',
    'drug': 'http://semanticscience.org/resource/SIO_010038',
}

# FUNCTIONS
def csv_to_rdf(edges_file_str, nodes_file_str, uri_type_dict):
    ''' This function turns a csv graph file into an RDF graph and saves it as a turtle file.
    Drug-gene interactions are removed to avoid bias in the drug-gene prediction.
    :param edges_file_str: the string of the edges csv of the graph to turn into rdf
    :param nodes_file_str: the string of the nodes csv of the graph to turn into rdf
    :param uri_type_dict: a dictionary of node types in the graph to uri links of these types
    :return: RDF graph in turtle format
    '''
    # open the graph files
    graph_nodes_file = pd.read_csv(nodes_file_str)

    # initialise the graph
    output_graph = Graph()

    input_file = csv.DictReader(open(edges_file_str))

    for row in input_file:
        # convert row from an OrderedDict to a regular dict
        row = dict(row)

        # structure edges in the dict
        subject_id = row['subject_id']
        property_uri = row['property_uri']
        object_id = row['object_id']
        property_label = row['property_label']

        # # if link is from drug to gene, remove link in rdf graph to avoid bias while predicting
        # subject_type = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'semantic_groups'].iloc[0]
        # object_type = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'semantic_groups'].iloc[0]
        # if (subject_type == 'drug') and (object_type == 'gene'):
        #     continue

        # get the uris for each node
        subject_uri = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'uri'].iloc[0]
        object_uri = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'uri'].iloc[0]

        # get subject and object type and turn into uri using uri_type_dict
        subject_type = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'semantic_groups'].iloc[0]
        try:
            subject_type_uri = uri_type_dict[subject_type]
        except:
            continue
        object_type = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'semantic_groups'].iloc[0]
        try:
            object_type_uri = uri_type_dict[object_type]
        except:
            continue

        # get subject and object label
        subject_label = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'preflabel'].iloc[0]
        object_label = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'preflabel'].iloc[0]

        # add types to the graph
        output_graph.add((URIRef(subject_uri), RDF.type, URIRef(subject_type_uri)))
        output_graph.add((URIRef(object_uri), RDF.type, URIRef(object_type_uri)))

        # add labels to the graph
        output_graph.add((URIRef(subject_uri), RDFS.label, Literal(subject_label)))
        output_graph.add((URIRef(object_uri), RDFS.label, Literal(object_label)))

        # add properties to the graph
        # output_graph.add((URIRef(property_uri), RDF.type, RDF.Property))
        # output_graph.add((URIRef(property_uri), RDFS.comment, Literal(property_label)))

        # add triples to the graph
        output_graph.add((URIRef(subject_uri), URIRef(property_uri), URIRef(object_uri)))

    # serialize and save the graph as turtle format
    # path to rdf/
    path = os.getcwd() + "/rdf"
    if not os.path.isdir(path): os.makedirs(path)
    output_graph.serialize(destination=f'./rdf/rdf_monarch_v{today}.ttl', format='ttl', encoding="utf-8")

    return output_graph.serialize(format="ttl")


if __name__ == '__main__':

    # create the rdf graph (where drug-gene interactions are not included)
    graph_edges_str = './carmen-graph/monarch/monarch_edges_disease_v2022-07-26.csv'
    graph_nodes_str = './carmen-graph/monarch/monarch_nodes_disease_v2022-07-26.csv'
    # create the RDF graph save as .ttl file
    csv_to_rdf(graph_edges_str, graph_nodes_str, uriType_dict)

