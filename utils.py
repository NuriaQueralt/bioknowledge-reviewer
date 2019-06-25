# @name: utils.py
# @description: Module for utils
# @version: 1.0
# @date: 21-01-2019
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TODO: add print functions
"""Module for utils"""


import datetime
import pandas as pd

# VARIABLES
today = datetime.date.today()


# FUNCTIONS

def get_dataframe(object):
    """
    This function converts a list_of_dictionaries object into a dataframe.
    :param object: graph object list
    :return: dataframe
    """

    try:
        df = pd.DataFrame(object)
    except ValueError:
        raise
    else:
        return df


def get_dataframe_from_file(filename):
    """
    This function opens a file and returns a dataframe.
    :param filename: CSV path_to_file_name string
    :return: dataframe
    """

    try:
        df = pd.read_csv('{}'.format(filename), low_memory=False)
    except OSError:
        print('cannot open: ', filename)
        print('Please, provide the correct file path and file name and the file in CSV format.')
        raise
    else:
        return df


def check_format(df, file_type='statements'):
    """
    This function checks if dataframe contains the expected columns before concatenation.
    :param df: dataframe
    :param file_type: concepts or statements (default value) string
    :return:
    """

    if file_type == 'concepts':
        try:
            df = df[['id', 'semantic_groups', 'preflabel', 'synonyms', 'name', 'description']]
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


def add_elem_dictionary2(dictionary, key, elem, repet = False):
    """
    This functions adds an element to a passed key and dictionary. \
    With the repet paramenter the user can set if value elements can be redundant in the value list ( repet = True )\
    or not ( repet = False ).
    :param dictionary: dictionary to update
    :param key: key string
    :param elem: new element string
    :param repet: False (as default value) or True
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


def add_one_dictionary2(dictionary, key):
    """
    This functions increases by one the counter value in the passed key and dictionary.
    :param dictionary: dictionary to update
    :param key: key string
    :return: updated dictionary
    """
    if key in dictionary:
        aux = dictionary.get(key)
        dictionary[key] = aux + 1 
    else:
        dictionary[key] = 1 
    return dictionary

    
def add_elem_with_dictionary(dictionary, key, elem, repeat = False):
    """
    This functions adds an element to a passed key and dictionary. \
    With the repet paramenter the user can set if value elements can be redundant in the value list ( repeat = True )\
    or not ( repeat = False ). This is a more efficient function than add_elem_dictionary2().
    :param dictionary: dictionary to update
    :param key: key string
    :param elem: new element string
    :param repeat: False (as default value) or True
    :return: updated dictionary
    """
    if not repeat:
        aux = dictionary.get(key, {}) 
        aux[elem] = 1 
        dictionary[key] = aux 
        return dictionary
    aux = dictionary.get(key, []) 
    aux.append(elem)
    dictionary[key] = aux 
    return dictionary


def add_one_dictionary(dictionary, key):
    """
    This functions increases by one the counter value in the passed key and dictionary. \
    This is a more efficient function than add_one_dictionary2().
    :param dictionary: dictionary to update
    :param key: key string
    :return: updated dictionary
    """
    aux = dictionary.get(key, 0)
    dictionary[key] = aux + 1 
    return dictionary



# def print_nodes(nodes, filename):
#     """This function save nodes into a CSV file."""

    # print output file

    #return

# if __name__ == '__main__':
    #seedList = ['HGNC:17646','HGNC:633']
    #neighbourList = get_neighbours_list(seedList)
    #print(len(neighbourList)) # 353
