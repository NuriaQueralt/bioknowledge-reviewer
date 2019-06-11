# @name: hypothesis.py
# @description: Module for hypothesis generation
# @version: 2.0
# @date: 24-02-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# This module is v2.0 because it is integrated with the rest of the graph library.
# This module reproduces orthopheno query on ngly1 graph created in the pipeline:
#      Database: neo4j-animal-v3 (localhost: http-7476, bolt-7689)
#      Seed: 'HGNC:17646','HGNC:633'  # HGNC ID scheme DOES NOT WORK
#      Query: clean orthopheno query
#      Output: json in hypothesis/ directory
#      PROBLEM: concept type
#      Sanity check:

#TODO:
#   * add check function to ensure that neo4j is up after import
#   * add check function to ensure that a query does not kill neo4j and that is up
#   * clean query (eliminate g1, ..)
#   * open query topology can be more open...
#   * add specific query topology as a function parameter
"""Module for hypothesis generation"""

from neo4j import GraphDatabase
import sys,os
import json
import yaml
import datetime
import neo4j.exceptions


# VARIABLES
today = datetime.date.today()


# QUERY MANAGEMENT FUNCTIONS


def parse_path( path ):
    """
    This function parses neo4j results.
    :param path: neo4j path object
    :return: parsed path dictionary
    """

    out = {}
    out['Nodes'] = []
    for node in path['path'].nodes:
        n = {}
        n['idx'] = node.id
        n['label'] = list(node.labels)[0]
        n['id'] = node.get('id')
        n['preflabel'] = node.get('preflabel')
        n['name'] = node.get('name')
        n['description'] = node.get('description')
        out['Nodes'].append(n)
    out['Edges'] = []
    for edge in path['path'].relationships:
        e = {}
        e['idx'] = edge.id
        e['start_node'] = edge.start_node.id
        e['end_node'] = edge.end_node.id
        e['type'] = edge.type
        e['preflabel'] = edge.get('property_label')
        e['references'] = edge.get('reference_uri')
        out['Edges'].append(e)
    return out


def query(genelist, queryname='', pwdegree='50', phdegree='20', format='json', port='7687'):
    """
    This function queries the graph database with the following specific ortho-phenotype query topology: \
    '(source:GENE)-[:`RO:HOM0000020`]-(:GENE)--(ds:DISO)--(:GENE)-[:`RO:HOM0000020`]-(:GENE)--(pw:PHYS)--(target:GENE)'\
    From an input gene list, it queries all possible source-target pairwise queries.
    :param genelist: gene list
    :param queryname: output file name string
    :param pwdegree: pathway node degree threshold string ('50' as default value)
    :param phdegree: phenotype node degree threshold string ('20' as default value)
    :param format: yaml or json output file format string  ('json' as default value)
    :param port: neo4j bolt port string ('7687' as default value)
    :return: None object
    """

    print('\nThe function "query()" is running...')
    # initialize neo4j
    try:
        driver = GraphDatabase.driver("bolt://localhost:{}".format(port), auth=("neo4j", "ngly1"))
    except neo4j.exceptions.ServiceUnavailable:
        raise

    # query topology
    query_topology = """
    (source:GENE)-[:`RO:HOM0000020`]-(:GENE)--(ds:DISO)--(:GENE)-[:`RO:HOM0000020`]-(:GENE)--(pw:PHYS)--(target:GENE)
    """

    # ask the driver object for a new session
    with driver.session() as session:
        # create outdir
        if not os.path.isdir('./hypothesis'): os.makedirs('./hypothesis')
        sys.path.insert(0, '.')

        # output filename
        filename = 'query_{}_pwdl{}_phdl{}_paths_v{}'.format(queryname, pwdegree, phdegree, today)

        # run query
        outputAll = list()
        for gene1 in genelist:
            for gene2 in genelist:
                if gene2 == gene1:
                    continue
                source = gene1
                target = gene2
                query = """
                MATCH path=""" + query_topology + """
                MATCH (source { id: '""" + source + """'}), (target { id: '""" + target + """'})
                // no loops or only one pass per node
                WHERE ALL(x IN nodes(path) WHERE single(y IN nodes(path) WHERE y = x))
                WITH ds, pw, path,
                     // count node degree
                     max(size( (pw)-[]-() )) AS pwDegree,
                     max(size( (ds)-[]-() )) AS dsDegree,
                     // mark general nodes to filter path that contain them out
                     [n IN nodes(path) WHERE n.preflabel IN ['cytoplasm','cytosol','nucleus','metabolism','membrane','protein binding','visible','viable','phenotype']] AS nodes_marked,
                     // mark promiscuous edges to filter path that contain them out,
                     [r IN relationships(path) WHERE r.property_label IN ['interacts with','in paralogy relationship with','in orthology relationship with','colocalizes with']] AS edges_marked
                // condition to filter paths that do contain marked nodes and edges out
                WHERE size(nodes_marked) = 0 AND size(edges_marked) = 0 AND pwDegree <= """ + pwdegree + """ AND dsDegree <= """ + phdegree + """
                RETURN path
                """
                result = session.run(query)
                pair = {}
                pair['source'] = source
                pair['target'] = target
                # parse query results
                output = list()
                counter = 0
                for record in result:
                    path_dct = parse_path(record)
                    output.append(path_dct)
                    counter += 1
                    if (counter % 100000 == 0):
                        sys.stderr.write("Processed " + str(counter) + "\n")
                pair['paths'] = output
                outputAll.append(pair)

                # print output
                path = os.getcwd() + '/hypothesis'
                if not os.path.isdir(path): os.makedirs(path)
                if ( format == "yaml" ):
                    with open('{}/{}.yaml'.format(path, filename), 'w') as f:
                        yaml.dump(outputAll, f, default_flow_style=False)
                    # print(yaml.dump(outputAll, default_flow_style=False))
                    print('\nThe query results are saved at: {}/{}.yaml\n'.format(path, filename))
                elif ( format == "json" ):
                    with open('{}/{}.json'.format(path, filename), 'w') as f:
                        json.dump(outputAll, f)
                    # print(json.dumps(outputAll))
                    print('\nThe query results are saved at: {}/{}.json\n'.format(path, filename))
                else:
                    sys.stderr.write("Error.\n")

    return sys.stderr.write("\nHypothesis generator has finished. {} QUERIES completed.\n".format(len(outputAll)))


def open_query(genelist, queryname='', format='json', port='7687'):
    """
    This function performs open queries to the graph database with the following query topology:\
    '(source:GENE)-[:`RO:HOM0000020`]-(:GENE)--(:DISO)--(:GENE)-[:`RO:HOM0000020`]-(:GENE)--(:PHYS)--(target:GENE)'
    From an input gene list, it queries all possible source-target pairwise queries.
    :param genelist: gene list
    :param queryname: output file name string
    :param format: yaml or json output file format string  ('json' as default value)
    :param port: neo4j bolt port string ('7687' as default value)
    :return: None object
    """

    print('\nThe function "open_query()" is running...')
    # initialize neo4j
    try:
        driver = GraphDatabase.driver("bolt://localhost:{}".format(port), auth=("neo4j", "ngly1"))
    except neo4j.exceptions.ServiceUnavailable:
        raise

    # query topology
    query_topology = """
    (source:GENE)-[:`RO:HOM0000020`]-(:GENE)--(:DISO)--(:GENE)-[:`RO:HOM0000020`]-(:GENE)--(:PHYS)--(target:GENE)
    """

    # ask the driver object for a new session
    with driver.session() as session:
        # create outdir
        if not os.path.isdir('./hypothesis'): os.makedirs('./hypothesis')
        sys.path.insert(0, '.')

        # output filename
        filename = 'query_{}_paths_v{}'.format(queryname, today)

        # run query
        outputAll = list()
        for gene1 in genelist:
            for gene2 in genelist:
                if gene2 == gene1:
                    continue
                source = gene1
                target = gene2
                query = """
                MATCH path=""" + query_topology + """
                MATCH (source { id: '""" + source + """'}), (target { id: '""" + target + """'})
                // no loops or only one pass per node
                //WHERE ALL(x IN nodes(path) WHERE single(y IN nodes(path) WHERE y = x))
                RETURN path
                """
                result = session.run(query)
                pair = {}
                pair['source'] = source
                pair['target'] = target
                # parse query results
                output = list()
                counter = 0
                for record in result:
                    path_dct = parse_path(record)
                    output.append(path_dct)
                    counter += 1
                    if (counter % 100000 == 0):
                        sys.stderr.write("Processed " + str(counter) + "\n")
                pair['paths'] = output
                outputAll.append(pair)

                # print output
                path = os.getcwd() + '/hypothesis'
                if not os.path.isdir(path): os.makedirs(path)
                if ( format == "yaml" ):
                    with open('{}/{}.yaml'.format(path, filename), 'w') as f:
                        yaml.dump(outputAll, f, default_flow_style=False)
                    # print(yaml.dump(outputAll, default_flow_style=False))
                    print('\nThe query results are saved at: {}/{}.yaml\n'.format(path, filename))
                elif ( format == "json" ):
                    with open('{}/{}.json'.format(path, filename), 'w') as f:
                        json.dump(outputAll, f)
                    # print(json.dumps(outputAll))
                    print('\nThe query results are saved at: {}/{}.json\n'.format(path, filename))
                else:
                    sys.stderr.write("Error.\n")

    return sys.stderr.write("\nHypothesis generator has finished. {} QUERIES completed.\n".format(len(outputAll)))


if __name__ == '__main__':
    seed = list([
        'HGNC:17646',  # NGLY1 human gene
        'HGNC:633'  # AQP1 human gene
    ])
    #seed=['HGNC:17646', 'HGNC:633']
    #seed=['NGLY1', 'AQP1']
    #query(seed,'ngly1_aqp1_port7689',port='7689')
    # seed=['NCBIGene:55768', 'NCBIGene:358']
    #query(seed, 'ngly1_aqp1_port7688', port='7688')
    #open_query(seed, queryname='ngly1_aqp1', port='7687')
