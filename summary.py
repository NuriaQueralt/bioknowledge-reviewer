# @name: summary.py
# @description: Module for hypothesis summarization
# @version: 2.0
# @date: 28-02-2018
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# This module is v2.0 because it is integrated with the rest of the graph library.
# It prints summaries in CSV format for consistency with the rest of the pipeline.

# TODO:
#   * change filename for outputs
#   * do a function than return table objects into dataframes, sort dataframe values by columns and print a table
#   * do a control for queries with 0 paths as result
#   * clean metapaths functions and distinguish work with dataframes
#   * do a execution block with the reduce function
#   * do private functions
"""Module for path summarization"""

import json
import pandas as pd
import os
import datetime


# VARIABLES
today = datetime.date.today()


# MANAGEMENT FUNCTIONS

def print_summaries(summary, filename):
    """
    This function saves the graph into a CSV file.
    :param summary: summary object list
    :param filename: filename string
    :return: None object
    """

    # print output file
    path = os.getcwd() + '/summaries'
    if not os.path.isdir(path): os.makedirs(path)
    pd.DataFrame(summary).to_csv('{}/{}_v{}.csv'.format(path,filename,today), index=False)

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path,filename,today))


def query_parser(query):
    """
    This function parses the query result object.
    :param query: result query object
    :return: parsed query object
    """

    print('\nThe function "query_parser()" is running...')
    source = query.get('source')
    target = query.get('target')
    #print('The query is {}-{}'.format(source, target))
    metapath_dict = dict()
    metapaths = list()
    entities = list()
    nodes = list()
    edges = list()
    path_idx = 0
    metapath_idx = 0
    for path in query.get('paths'):
        path_idx += 1
        nodes_list = path.get('Nodes')
        for node in nodes_list:
            node['object_order'] = nodes_list.index(node) * 2
            node['path_idx'] = path_idx
            nodes.append(node)
        edges_list = path.get('Edges')
        for edge in edges_list:
            edge['object_order'] = edges_list.index(edge) * 2 + 1
            edge['path_idx'] = path_idx
            edge['label'] = edge['preflabel']
            edge['id'] = edge['type']
            edges.append(edge)
        # metapath idx
        objects = list()
        [objects.append(object) for object_list in [nodes_list, edges_list] for object in object_list]
        #
        #metapath_list = list()
        #[metapath_list.append(object['label']) for i in range(len(objects)) for object in objects if
        # object['object_order'] == i]
        #metapath_str2 = '-'.join(metapath_list)
        #print('str2: {}'.format(metapath_str))
        #
        metapath_list = list()
        [metapath_list.append(object) for i in range(len(objects)) for object in objects if object['object_order'] == i]
        label_list = list()
        [label_list.append(object['label']) for object in metapath_list]
        metapath_str = '-'.join(label_list)
        #print('str: {}'.format(metapath_str2))
        if metapath_str not in metapath_dict:
            metapath_idx += 1
            metapath_dict[metapath_str] = metapath_idx
        path_entities_l = list()
        #for object in objects:
        for object in metapath_list:
            object['metapath_idx'] = metapath_dict[metapath_str]
        # create a uniform path object
            entity = dict()
            entity['idx'] = object['idx']
            entity['label'] = object['label']
            entity['id'] = object['id']
            entity['preflabel'] = object['preflabel']
            entity['path_idx'] = object['path_idx']
            entity['metapath_idx'] = object['metapath_idx']
            entity['object_order'] = object['object_order']
            path_entities_l.append(entity)
            entities.append(entity)
        # create metapath obj in path records
        metapath = {
            'metapath_idx': entity['metapath_idx'],
            'path_idx': entity['path_idx'],
            'label': metapath_str,
            'entities': path_entities_l,
            'length': len(path_entities_l)
        }
        metapaths.append(metapath)

    # add metapaths, entities, nodes and edges attributes
    query['metapaths'] = metapaths
    query['entities'] = entities
    query['nodes'] = nodes
    query['edges'] = edges
    print('\nFinished query_parser().\n')

    return query


def path_load(filename):
    """
    This function loads paths from a json file to a digital object.
    :param filename: path to JSON file string
    :return: loaded json object
    """
    return json.load(open('{}'.format(filename),'r'))


def path_count(query):
    """
    This function returns the total number of paths per query passed.
    :param query: query object
    :return: number of paths integer
    """
    return len(query.get('paths'))


def metapath_summarization(query):
    """
    This function returns metapaht summary objects.
    :param query: query object
    :return: metapath summary object, metapath label summary object
    """

    # functions from entities
    df = pd.DataFrame(query.get('entities'))
    metapath_df = df.groupby(['metapath_idx','path_idx'])['idx'].count().reset_index().groupby(['metapath_idx','idx'])['path_idx'].count().reset_index().rename(columns={'idx': 'metapath_length', 'path_idx': 'metapath_count'})
    metapath_label_df = df.groupby(['metapath_idx','path_idx'])['label'].apply(lambda a: '\t'.join(a)).reset_index().groupby(['metapath_idx','label'])['path_idx'].apply(lambda a: 'path').reset_index()[['metapath_idx','label']].rename(columns={'label': 'metapath_label'})
    return metapath_df.to_dict('records'), metapath_label_df.to_dict('records')


def metapath_count(metapath_idx, metapaths_l, attribute='count'):
    """
    This function returns two metapath counts depending on the attribute setting: \
    the total number of metapaths ( attribute = 'count' ) or the metapath length ( attribute = 'length' )
    :param metapath_idx: metapath index integer
    :param metapaths_l: metapaths list
    :param attribute: 'count' (as default value) or 'length'
    :return: count integer
    """

    for mp in metapaths_l:
        if mp.get('metapath_idx') == metapath_idx:
            if attribute == 'length':
                return mp.get('metapath_length')
            else:
                return mp.get('metapath_count')


def metapath_label(metapath_idx, metapaths_l):
    """
    This function returns the metapath_label
    :param metapath_idx: metapath index integer
    :param metapaths_l: metapaths list
    :return: metapath label string
    """

    for mp in metapaths_l:
        if mp.get('metapath_idx') == metapath_idx:
            return mp.get('metapath_label')


def metapath_count2(metapath_idx, metapaths_l):
        """
        This function returns the metapath_count.
        :param metapath_idx: metapath index integer
        :param metapaths_l: metapaths list
        :return: count integer
        """

    #for mp in metapaths_l:
        return len(list(filter(lambda x: x.get('metapath_idx') == metapath_idx, metapaths_l)))


def metapath_label2(metapath_idx, metapaths_l):
    """
    This function returns the metapath_label.
    :param metapath_idx: metapath index integer
    :param metapaths_l: metapaths list
    :return: metapath label string
    """

    for mp in metapaths_l:
        if mp.get('metapath_idx') == metapath_idx:
            return mp.get('label').replace('-','\t')


def metapaths(data):
    """
    This function prepares metapath summary table.
    :param data: parsed queries list
    :return: None object
    """

    print('\nThe function "metapaths()" is running...')
    for query in data:
        #print()
        #print(query.get('source'), query.get('target'))

        # table 1: user profile
        # format : path_count | metapath_idx | metapath_count | obj_1 LABEL | ... | obj_n LABEL (rows = metapaths)
        # rows = metapaths
        # functions:
        ## from 'entities'
        #metapaths_counts_l, metapaths_label_l = metapath_summarization(query)
        ##print('path_count\tmetapath_idx\tmetapath_count\tmetapath_label')
        #for metapath_idx in {entity.get('metapath_idx') for entity in query.get('entities')}:
        #    print()
        #    #print('{}\t{}\t{}\t{}\n'.format(path_count(query), metapath_idx, metapath_count(metapath_idx, metapaths_counts_l, attribute='count'), metapath_label(metapath_idx, metapaths_label_l)))

        # from 'metapaths'
        #print('path_count\tmetapath_idx\tmetapath_count\tmetapath_label')
        print('\nPrinting summaries...\n')
        metapath_dct = dict()
        metapath_l = list()
        for metapath_idx in { mp.get('metapath_idx') for mp in query.get('metapaths') }:
            metapath_dct['path_count'] = path_count(query)
            metapath_dct['metapath_idx'] = metapath_idx
            metapath_dct['metapath_count'] = metapath_count2(metapath_idx, query.get('metapaths'))
            metapath_dct['metapath_label'] = metapath_label2(metapath_idx, query.get('metapaths'))
            metapath_l.append(dict(metapath_dct))
            #print('{}\t{}\t{}\t{}\n'.format(path_count(query), metapath_idx,
            #                                metapath_count2(metapath_idx, query.get('metapaths')),
            #                                metapath_label2(metapath_idx, query.get('metapaths'))))
        print_summaries(metapath_l,
                        filename='query_source:{}_target:{}_summary_metapaths'.format(
                            query.get('source'), query.get('target')))

        # table 2: bioinformatic profile
        # format : metapath_idx | object_order | object_type | metapath_counts | metapath_length (rows = entities)
        # rows = query.get('entities')
        # functions: { metapath_idx: metapath_counts }, { metapath_idx: metapath_length }
        ## from 'entities'
        ##print('table 2 - entities')
        #for entity in query.get('entities'):
        #    print()
        #    #print(entity.get('metapath_idx'), entity.get('object_order'),
        #    #      entity.get('label'), metapath_count(entity.get('metapath_idx'), metapaths_counts_l, attribute='count'),
        #    #      metapath_count(entity.get('metapath_idx'), metapaths_counts_l, attribute='length'))

        # from 'metapaths'
        #print('table 2 - metapaths')
        entity_dct = dict()
        entity_l = list()
        for mp in query.get('metapaths'):
            metapath_idx = mp.get('metapath_idx')
            for entity in mp.get('entities'):
                entity_dct['metapath_idx'] = metapath_idx
                entity_dct['object_order'] = entity.get('object_order')
                entity_dct['object_type'] = entity.get('label')
                entity_dct['metapath_count'] = metapath_count2(metapath_idx, query.get('metapaths'))
                entity_dct['metapath_length'] = len(mp.get('entities'))
                entity_l.append(dict(entity_dct))
                #print(mp.get('metapath_idx'), entity.get('metapath_idx'), entity.get('object_order'), entity.get('label'), metapath_count2(mp.get('metapath_idx'), query.get('metapaths')), len(mp.get('entities')), mp.get('length'))

        print_summaries(entity_l,
                        filename='query_source:{}_target:{}_summary_entities_in_metapaths'.format(
                            query.get('source'), query.get('target')))

    print('\nFinished metapaths().\n')

def node(idx, nodes_l):
    """
    This function returns the node/edge dictionary data object.
    :param idx: node/edge index integer
    :param nodes_l: nodes/edges list
    :return: node/edge dictionary
    """

    for nodus in nodes_l:
        if nodus.get('idx') == idx:
            return nodus


def edges_count(paths):
    """
    This function returns edge patterns such as 'start-edge-end' in path results as a list of tuples.
    It also returns the map between the tuple {start, edge, end} and the pattern.
    :param paths: paths query object
    :return: edges count dictionary
    """

    edges_l = list()
    edge_pattern_dct = dict()
    for path in paths:
        for node_i in range(len(path.get('Nodes'))):
            start = path.get('Nodes')[node_i].get('preflabel') + '::' + path.get('Nodes')[node_i].get('id')
            edge = path.get('Edges')[node_i].get('preflabel') + '::' + path.get('Edges')[node_i].get('type')
            end = path.get('Nodes')[node_i + 1].get('preflabel') + '::' + path.get('Nodes')[node_i + 1].get('id')
            key = start + '__' + edge + '__' + end
            edge_pattern_dct[key] = {start, edge, end}
            edges_l.append({start, edge, end})
            if node_i == len(path.get('Nodes')) - 2:
                break

    # count
    count_dct = dict()
    for key in edge_pattern_dct:
        count_dct[key] = len(list(filter(lambda edge: edge == edge_pattern_dct[key], edges_l)))

    return count_dct


def nodes_count(idx, nodes_l, attribute = 'idx'):
        """
        This function returns node/edge counts. With the attribute parameter, it counts by index setting 'idx' \
        or by label setting 'label'.
        :param idx: node/edge index integer
        :param nodes_l: nodes/edges list
        :param attribute: 'idx' (as default value) or 'label'
        :return: count integer
        """

        if attribute == 'idx':
            return len(list(filter(lambda x: x.get('idx') == idx, nodes_l)))
        if attribute == 'label':
            return len(list(filter(lambda x: x.get('label') == idx, nodes_l)))


def nodes(data):
    """
    This function prepares nodes summary table.
    :param data: parsed queries list
    :return: None object
    """

    print('\nThe function "nodes()" is running...')
    print('\nPrinting summaries...\n')
    for query in data:
        #print()
        #print(query.get('source'), query.get('target'))

        # table 1: user profile
        # format : path_count | node_type | node_value | node_count (rows = nodes)
        #print('path_count\tnode_type\tnode_value\tnode_count')
        node_sum_dct = dict()
        nodes_sum_l = list()
        for node_idx in {item.get('idx') for item in query.get('nodes')}:
            node_dct = node(node_idx, query.get('nodes'))
            node_sum_dct['path_count'] = path_count(query)
            node_sum_dct['node_type'] = node_dct.get('label')
            node_sum_dct['node_value'] = node_dct.get('preflabel') + '::' + node_dct.get('id')
            node_sum_dct['node_count'] = nodes_count(node_idx, query.get('nodes'), attribute = 'idx')
            nodes_sum_l.append(dict(node_sum_dct))
            #print('{}\t{}\t{}\t{}\n'.format(path_count(query), node_dct.get('label'),
            #                                   node_dct.get('preflabel') + '::' + node_dct.get('id'),
            #                                   nodes_count(node_idx, query.get('nodes'), attribute = 'idx')))

        print_summaries(nodes_sum_l,
                        filename='query_source:{}_target:{}_summary_nodes'.format(
                            query.get('source'), query.get('target')))

    print('\nFinished nodes().\n')

    #return nodes_sum_l


def node_types(data):
    """
    This function prepares node types summary table.
    :param data: parsed queries list
    :return: None object
    """

    print('\nThe function "node_types()" is running...')
    print('\nPrinting summaries...\n')
    for query in data:
        #print()
        #print(query.get('source'), query.get('target'))

        # table 1: user profile
        # format : path_count | node_type | node_type_count (rows = node types)
        #print('path_count\tnode_type\tnode_type_count')
        node_sum_dct = dict()
        nodes_sum_l = list()
        for node_label in {item.get('label') for item in query.get('nodes')}:
            node_sum_dct['path_count'] = path_count(query)
            node_sum_dct['node_type'] = node_label
            node_sum_dct['node_type_count'] = nodes_count(node_label, query.get('nodes'), attribute = 'label')
            nodes_sum_l.append(dict(node_sum_dct))
            #print('{}\t{}\t{}\n'.format(path_count(query), node_label,
            #                            nodes_count(node_label, query.get('nodes'), attribute = 'label')))

        print_summaries(nodes_sum_l,
                        filename='query_source:{}_target:{}_summary_node_types'.format(
                            query.get('source'), query.get('target')))

    print('\nFinished node_types().\n')
    #return nodes_sum_l


def edges(data):
    """
    This function prepares edges summary table."
    :param data: parsed queries list
    :return: None object
    """

    print('\nThe function "edges()" is running...')
    print('\nPrinting summaries...\n')
    for query in data:
        #print()
        #print(query.get('source'), query.get('target'))

        # table 1: user profile
        # format : path_count | node_type | node_value | node_count (rows = nodes)
        #print('path_count\tsubject_value\tedge_type\tobject_value\tedge_count')
        edge_sum_dct = dict()
        edges_sum_l = list()
        edge_sum_dct['path_count'] = path_count(query)
        edge2count_dict = edges_count(query.get('paths'))
        for edge_pattern in edge2count_dict:
            start,edge,end = edge_pattern.split('__')
            edge_sum_dct['subject_value'] = start
            edge_sum_dct['edge_type'] = edge
            edge_sum_dct['object_value'] = end
            edge_sum_dct['edge_count'] = edge2count_dict[edge_pattern]
            edges_sum_l.append(dict(edge_sum_dct))
            #print('{}\t{}\t{}\t{}\t{}\n'.format(path_count(query), edge_sum_dct.get('subject_value'),
            #                                    edge_sum_dct.get('edge_type'), edge_sum_dct.get('object_value'),
            #                                    edge_sum_dct['edge_count']))

        print_summaries(edges_sum_l,
                        filename='query_source:{}_target:{}_summary_edges'.format(
                            query.get('source'), query.get('target')))

    print('\nFinished edges().\n')
    #return edges_sum_l


def edge_types(data):
    """
    This function prepares edge types summary table.
    :param data: parsed queries list
    :return: None object
    """

    print('\nThe function "edge_types()" is running...')
    print('\nPrinting summaries...\n')
    for query in data:
        #print()
        #print(query.get('source'), query.get('target'))

        # table 1: user profile
        # format : path_count | node_type | node_type_count (rows = node types)
        #print('path_count\tedge_type\tedge_type_count')
        edge_sum_dct = dict()
        edges_sum_l = list()
        for edge_label in {item.get('label') for item in query.get('edges')}:
            edge_sum_dct['path_count'] = path_count(query)
            edge_sum_dct['edge_type'] = edge_label
            edge_sum_dct['edge_type_count'] = nodes_count(edge_label, query.get('edges'), attribute = 'label')
            edges_sum_l.append(dict(edge_sum_dct))
            #print('{}\t{}\t{}\n'.format(path_count(query), edge_label,
            #                            nodes_count(edge_label, query.get('edges'), attribute = 'label')))
            #

        print_summaries(edges_sum_l, filename='query_source:{}_target:{}_summary_edge_types'.format(query.get('source'),query.get('target')))

    print('\nFinished edge_types().\n')
    #return edges_sum_l


if __name__ == "__main__":

    #data = path_load('out/query_1_pwdl20_phdl20_paths_v2018-02-25')
    ##data = path_load('out/q1_1_in0_pwdl50_phdl20_paths')

    #data_parsed = list()
    #funcs = [metapaths, nodes, node_types, edges, edge_types]
    #for query in data:
     #   query_parsed = query_parser(query)
        ##metapath(query_parsed)
        ##map(lambda x: x(query_parsed), funcs)
       # data_parsed.append(query_parsed)
    #metapaths(data_parsed)
    ##nodes(data_parsed)
    ##node_types(data_parsed)
    ##edges(data_parsed)
    ##edge_types(data_parsed)
    ##for query in data_parsed:
    ##    map(lambda x: x(query), funcs)
    ## prints

    # 2019
    # data = path_load('./hypothesis/query_ngly1_aqp1_paths_v2019-02-21')
    data = path_load('./hypothesis/query_ngly1_aqp1_paths_v2019-03-06')
    data_parsed = list()
    for query in data:
        query_parsed = query_parser(query)
        data_parsed.append(query_parsed)

    metapaths(data_parsed)

