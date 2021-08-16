#!/usr/bin/python
import glob
# Import everything from matplotlib (numpy is accessible via 'np' alias)
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
import os
import mat73


rawpath = '/Users/pcs/data/ABF/rawdata/nu2-tm'
savepath = '/Users/pcs/data/ABF/analysis/nu2-tm' 

#### Reading the file with the averaged mean ROAG ########
filename = '/nu2-tm-sem-det-py-form.mat'

data_dict = mat73.loadmat(rawpath+filename)
#weak_xi = data_dict['data']
