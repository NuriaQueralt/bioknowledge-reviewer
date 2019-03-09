# @name: utils.py
# @description: Module for utils
# @version: 1.0
# @date: 21-01-2019
# @author: Núria Queralt Rosinach
# @email: nuriaqr@scripps.edu

# TO DO:

"""Module for utils"""


import datetime
import pandas as pd
#import os, glob, sys
#from gsheets import Sheets
#from biothings_client import get_client
#sys.path.insert(0,'/home/nuria/soft/utils3/lib/')
#import abravo_lib as utils

# VARIABLES
today = datetime.date.today()


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
        df = pd.read_csv('{}'.format(filename), low_memory=False)
    except OSError:
        print('cannot open: ', filename)
        print('Please, provide the correct file path and file name and the file in CSV format.')
        raise
    else:
        return df


def check_format(df, file_type='statements'):
    """This function checks if dataframe contains the expected columns before concatanation."""

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


# def print_nodes(nodes, filename):
#     """This function save nodes into a CSV file."""

    # print output file

    #return

# if __name__ == '__main__':
    #seedList = ['HGNC:17646','HGNC:633']
    #neighbourList = get_neighbours_list(seedList)
    #print(len(neighbourList)) # 353