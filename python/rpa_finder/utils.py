import numpy as np
import re


def extract_string(str):
    pattern = r'\"([a-zA-Z0-9_\[\]\{\}\^]+)\"'    
    matches = re.findall(pattern, str)
    if (len(matches)==0): 
        return ""
    return matches[0]

def extract_number(str):
    pattern = r'\d+(?=(?:[^"]*"[^"]*")*[^"]*$)'
   
    matches = re.findall(pattern, str)
    if( len(matches)==0 ):
        return 1    
    return int(matches[0])

def delete_all_blanks(str):
    return re.sub(r'\s+', '', str)

#  '2 "A"' -> ["A", 2]
def to_pair(str):
    return [extract_string(str), extract_number(str)]

def from_str_to_list(str):
    str = delete_all_blanks(str)
    str = re.split('\+', str)
    return [to_pair(s) for s in str]

def extract_species(r_list):
    r_list = [item for sublist in r_list for item in sublist]
    r_list = np.reshape(np.array(r_list), (-1,2))
    species = np.unique(r_list[:,0], axis=0)
    species = np.sort(species)
    species = delete_elements(species, '')
    return species

def delete_elements(arr, element):
    return arr[arr != element]

def delete_elements_with_criteria(arr, criteria):
    print('arr')
    print(arr)
    vectorized_criteria = np.vectorize(criteria)
    to_delete = vectorized_criteria(arr)
    return arr[~to_delete]

def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]

def parse_string_to_reaction_list(str):
    lines = str.strip('\n').split('\n')
    lines = [l for l in lines if l]
    ret = [ l.split(',') for l in lines]
    return ret
