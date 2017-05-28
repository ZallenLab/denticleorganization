
# coding: utf-8

# In[1]:

from __future__ import division

from sys import argv
import glob
import os
import numpy as np
import pandas as pd
import scipy.stats as sps
import time
import string
import csv


# Generate uniform random distributions based on the number of cells given

# In[ ]:




# In[ ]:




# In[4]:

cellstosim = [(1,1284),(2,1140),(3,476),(4,130)]
iterations = 10000

for elem in cellstosim: 
    dent, cells = elem

    positions = np.zeros(((cells*dent),iterations))
    fname = str(dent)+'_montecarlo_positions_replicates.csv'
    
    for it in range(0,iterations): 
        this = np.reshape(np.random.rand(cells, dent),(1,-1))       
        positions[:,it] = this
        

    np.savetxt(fname, positions, delimiter=',')


# In[ ]:




# calculate KS test data, and count how many tests pass for each dentincell number (output in summarydata.csv file)

# In[25]:

def TestPasses(pval, cutoff):
    if pval <= cutoff: 
        return 'different'
    elif pval > cutoff: 
        return 'same'


def IndivStatTest(simdata, filename_out):
# IN: 3D np array, list of strings with length=arr[X,:,:] (array axis 0), name of csv file

    test_ks = sps.ks_2samp(invivo_d, simdata)
    # outputs [ks-statistic, p-value]

    with open(filename_out, 'a') as f:
        csv.writer(f).writerows([[column,  test_ks[0], test_ks[1], TestPasses(test_ks[1], 0.05)]])
        
    return  test_ks[1], TestPasses(test_ks[1], 0.05)
            

dicmap = ['null','A','B','C','D']


invivo_file = 'yw_all_RelativePosition.csv'
dentnumbers = [1,2,3,4]


invivo_data = pd.read_csv(invivo_file)

for dentincell in dentnumbers: 
    # clear out missing data
    invivo = invivo_data[dicmap[dentincell]]
    invivo = invivo.replace(0,np.nan)       # turn zeros into NaNs
    invivo = invivo.dropna(how='all')       # drop any column (axis=0) or row (axis=1) where ALL values are NaN

    invivo_d = invivo/100

    mcname = str(dentincell)+'_montecarlo_positions_replicates.csv'
    sfname = 'summarydata.csv'

    montecarlo = pd.read_csv(mcname,header=None)
    pf = []

    for column in montecarlo: 
        pval, dif = IndivStatTest(montecarlo[column], 'montecarlo_kstests_'+str(dentincell)+'dent.csv')
        
        pf.append(dif)

    pfr = pd.Series(pf)
    with open(sfname,'a') as f:
        f.write(str(dentincell) + ',' + str(pfr[pfr == 'same'].count()) + ',\n')


# In[18]:

pfr = pd.Series(pf)
with open(sfname,'a') as f:
    f.write(str(dentincell) + ',' + str(pfr[pfr == 'same'].count()) + ',\n')


# In[22]:

invivo_data


# In[ ]:

pf


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')

mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)
hist, bins = np.histogram(x, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.show()
# In[ ]:

hist, bins = np.histogram(positions,bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



