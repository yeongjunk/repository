import pandas as pd
import numpy as np
import glob, os

"""
Directory of the fullscan folder
"""
dir = '/Users/pcs/data/ABF3D/fullscan/'


# file name generator
L = [10, 15, 20, 25]
R = [400, 150, 50, 15]
Th = range(1, 27)

# 'i' th size , 'j' theta (start from 1), 'r' realization (start from 1)
i = 0
j = 25
r = 2

fn = 'L%d/L%d_Th%d_R%d.csv'%(L[i], L[i], Th[j], R[i])


df = pd.read_csv(dir + fn)
mask = (df.r == r)
df_r = df.loc[mask,:]
E_r = df_r.E

E_r
