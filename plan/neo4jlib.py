# @name: neo4jlib.py
# @description: Module for Neo4j data structure and management
# @version: 1.0
# @date: 22-02-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# All files before concatenation should have the same format:
# Statements format:
# :START_ID,
# :TYPE,
# :END_ID,
# reference_uri,
# reference_supporting_text,
# reference_date,
# property_label,
# property_description,
# property_uri

# Concepts format:
# id:ID,
# :LABEL,
# preflabel,
# synonyms:IGNORE,
# description

# TO DO:
# code
# modularize
# debugging
# documentation
"""Module for Neo4j"""

import pandas as pd
import datetime
import os, sys
import subprocess

# VARIABLES
today = datetime.date.today()

# CHECK GRAPH SCHEMA AND NORMALIZE TO NEO4J FORMAT

# BUILD NETWORK

# NETWORK MANAGEMENT FUNCTIONS


#### UTILS

def get_dataframe(object):
    """This function converts a list_of_dictionaries object into a dataframe."""

    try:
        df = pd.DataFrame(object)
    except ValueError:
        raise
    else:
        return df

def get_dataframe_from_file(filename):
    """This function opens a file and returns a dataframe."""

    try:
        df = pd.read_csv('{}.csv'.format(filename))
    except OSError:
        print('cannot open: ', filename)
        raise
    else:
        return df


def check_format(df, file_type='statements'):
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


def save_neo4j_files(object, neo4j_path, file_type = 'statements'):#(graph_l, ):
    """This function save the neo4j graph files in CSV format into the neo4j import directory."""

    graph_version = 'v{}'.format(today)
    path_to_import = neo4j_path + '/import/ngly1'
    path_to_version = neo4j_path + '/import/ngly1/' + graph_version
    for dir in path_to_import, path_to_version:
        if not os.path.isdir(dir): os.makedirs(dir)
    # save filling null and with sep=','
    if file_type == 'statements':
        object.to_csv('{}/ngly1_statements.csv'.format(path_to_import),
                          index=False, na_rep='NA')
        object.to_csv('{}/ngly1_statements.csv'.format(path_to_version),
                          index=False, na_rep='NA')
        return print("\nFile '{}/ngly1_statements.csv' saved.".format(path_to_import))
    elif file_type == 'concepts':
        object.to_csv('{}/ngly1_concepts.csv'.format(path_to_import),
                        index=False, na_rep='NA')
        object.to_csv('{}/ngly1_concepts.csv'.format(path_to_version),
                    index=False, na_rep='NA')
        return print("\nFile '{}/ngly1_concepts.csv' saved.".format(path_to_import))
    else:
        return print('The user should provide the "file_type" argument with any of the [statements or concepts] value.')


### FORMAT NEO4J

def get_statements(edges):
    """This function returns the statements neo4j file format."""

    # check header
    statements = check_format(edges)
    # format header
    statements = (
        statements
            .rename(columns={
            'subject_id': ':START_ID',
            'property_id': ':TYPE',
            'object_id': ':END_ID',
            'property_description': 'property_description:IGNORE'
        }
        )
    )

    return statements


def get_concepts(nodes):
    """This function returns the concepts neo4j file format."""

    # check header
    concepts = check_format(nodes, file_type='concepts')
    # format header
    concepts = (
        concepts
            .rename(columns={
            'id': 'id:ID',
            'semantic_groups': ':LABEL',
            'synonyms': 'synonyms:IGNORE'
        }
        )
    )

    return concepts


def do_import(neo4j_path):
    """This function executes the import of the graph into the neo4j server instance."""

    try:
        #graph_version = 'v{}'.format(today)
        path_to_import = neo4j_path + '/import/ngly1'
        # stop neo4j
        cmd = '{}/bin/neo4j stop'.format(neo4j_path)
        subprocess.call(cmd, shell=True)
        # rm any database in the database dir
        cmd = 'rm -rf {}/data/databases/graph.db/*'.format(neo4j_path)
        subprocess.call(cmd, shell=True)
        # cd import dir files path
        cmd = 'cd {}'.format(path_to_import)
        subprocess.call(cmd, shell=True)
        #  neo4j-import
        cmd = '{}/bin/neo4j-import --into {}/data/databases/graph.db --id-type string --nodes {}/ngly1_concepts.csv --relationships {}/ngly1_statements.csv'.format(neo4j_path, neo4j_path, path_to_import, path_to_import)
        subprocess.call(cmd, shell=True)
        # start neo4j from database dir
        cmd = 'cd {}/data/databases/graph.db'.format(neo4j_path)
        subprocess.call(cmd, shell=True)
        cmd = '{}/bin/neo4j start'.format(neo4j_path)
        subprocess.call(cmd, shell=True)
    except:
        print('error: {}'.format(sys.exc_info()[0]))
        raise
    else:
        return print('\nThe graph is imported into the server. The server is running.\n')


if __name__ == '__main__':
    ## get edges and files for neo4j
    #edges = get_dataframe('./graph/monarch_edges_v2018-02-24')
    #nodes = get_dataframe('./graph/monarch_nodes_v2018-02-23')
    #statements = get_statements(edges)
    #concepts = get_concepts(nodes)

    ## import the graph into neo4j
    # save files into neo4j import dir
    neo4j_path = '/home/nuria/workspace/ngly1-graph/neo4j-community-3.0-animalModel-v3'
    #save_neo4j_files(statements, neo4j_path, file_type = 'statements')
    #save_neo4j_files(concepts, neo4j_path, file_type = 'concepts')

    # import graph into neo4j
    do_import(neo4j_path)
