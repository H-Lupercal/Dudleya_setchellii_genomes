# Query the DU# - Code - Population dictionary to get the code & population of a given sample

# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# December 2025

import os
import argparse
import json

# This may change - update as needed
DICT_PATH = '/WAVE/projects/whittalllab/dudleya/snakemake-testing/samples_config/du#_pop_dictionary.json'
DICT_FILE = open(DICT_PATH, 'r')

# Current format of the dictionary:
'''
[
  {
    "DU#": 10,
    "Code": "UVAC",
    "Population": "Little Uvas Creek"
  },
  ...
'''

def query_dictionary(du_number):
  '''
    input:
      (int) DU# to query the dictionary for
    outputs:
      (str) Code - the shorthand population code for the sample
      (str) Population - the full name of the population the sample is from
  '''
    pop_dict = json.load(DICT_FILE)
    try:
      sample = [s for s in pop_dict if s["DU#"] == du_number][0]
    except:
      print(f"Error: could not find DU# '{du_number}' in the dictionary.")
    
    # print(sample)
    code = sample["Code"]
    population = sample["Population"]
    return(code, population)


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser()
    parser.add_argument('du_number', help='Dudleya sample number (e.g. 173 or 59)', type=int)
    args = parser.parse_args()
    # Function call
    code, population = query_dictionary(args.du_number)
    print(f"Code: {code}")
    print(f"Population: {population}")