import numpy as np
import pandas as pd
import glob, os


"""
return the list of csv files (full directory) inside a given directory
"""
def list_filenames(directory):
    fns = []
    os.chdir(directory)
    for file in glob.glob('*.csv'):
        fns.append(directory+'/'+file)
    return fns

"""
read pandas dataframe csv file and convert it into numpy array
"""
def load_as_nparray(filename):
    df = pd.read_csv(filename)
    data = df.to_numpy()
    return data

size = ['L20', 'L40', 'L60', 'L80', 'L100', 'L120', 'L140']

# Total number of realizations for each size
R ={'L20' : 500,
    'L40' : 400,
    'L60' : 300,
    'L80': 200,
    'L100' : 30,
    'L120': 20,
    'L140': 15
    }


#-----------variables to modify-----------#
raw_dir = '/Users/pcs/data/ABF2D/'
i = 0 # ith size
j = 0 # jth realizations
#-----------------------------------------#


filename = {}
for i in range(len(size)):
    filename[size[i]] = list_filenames(raw_dir + size[i])


data = load_as_nparray(filename[size[i]][j])
print(data)
