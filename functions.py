import re
import numpy as np
import pandas as pd
import os
import json
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')


def get_chemicals_list(filename):
    chemicals = pd.read_csv(filename)
    chemicals_list = list(chemicals[' symbol'].str.strip())
    
    return chemicals_list


def parse_chemical_name(name, chemicals_list):
    temp = ''
    stoi = []
    chem = []
    a = re.findall('[A-Z][^A-Z]*', name)
    
    for ele in a:
        if ele[-1].isdigit():
            stoi.append(float(ele[-1]))
        else:
            stoi.append(1)
            
        chem.append(''.join(i for i in ele if not i.isdigit()))
        
    return chem, stoi


def check_artifacts_name(name, chem_list):
    if '(' and ')' in name:
        return filter_name(name, chem_list)
    else:
        return name


def filter_name(name, chem_list):
    output = ''
    temp = ''
    mult = 1
    i = 0
    j = 0
    count = 0
    k = 0
    for char in name:
        if char == '(':
            i = count
        if char == ')':
            j = count
            if j < len(name) - 1:
                if name[j+1].isdigit():
                    mult = name[j+1]
                    k = j + 1
        count += 1
        
    output = name[:i] + name[k+1:]    
    
    temp = name[i+1:j]    
    
    fix, stoi = parse_chemical_name(temp, chem_list)
    
    for el in fix:
        if int(mult) > 1:
            output += el + mult
        else:
            output += el
    
    return output


def elastic(properties):
        
    for i in range(len(properties)):
        if properties[i]['name'] == 'Stiffness Tensor':
            return True
    return False


def extract_elastic_data():
    output = {}
    
    for filename in os.listdir():
        if '.json' in filename:
            print(filename)
            data = pd.read_json(filename)
            for i in range(len(data)):
                sample = data.iloc[i]

                if elastic(sample.properties):
                    output[sample.chemicalFormula] = sample.properties
    return output


def extract_property(data, prop, units = True):
    output = {}
    
    for key in data.keys():
        sample = data[key]
        for i in range(len(sample)):
            if sample[i]['name'] == prop:
                if units:
                    output[key] = {'scalar' : sample[i]['scalars'], 'units' : sample[i]['units']}
                else:
                    output[key] = {'scalar' : sample[i]['scalars']}
    print('number of samples: ' + str(len(output.keys())))
    
    return output


def check_all_equal_units(bulk, unit):
    unit = 'GPa'
    for key in bulk.keys():
        sample = bulk[key]
        if sample['units'] != unit:
            print(str(key) + ' has different units')
            break
    
    print('OK - all samples have same units')


def plot_histogram(data, prop):
    values = []
    for key in data.keys():
        values.append(data[key]['scalar'])

    fig = plt.figure()
    plt.hist(values, bins='auto', range=[0, 400])  # arguments are passed to np.histogram
#     plt.title("Data Distribution")
    plt.xlabel('Elastic Modulus (GPa)')
#     plt.xlabel(prop)
    plt.ylabel('Frequency')

    return fig


def check_negative_values(data):
    negatives = []
    for key in data.keys():
        if data[key]['scalar'] < 0:
            negatives.append(key)

    print('There is ' + str(len(negatives)) + ' materials with negative property value')
    print(negatives)