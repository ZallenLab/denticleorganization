
# coding: utf-8

# In[188]:

import pandas as pd
import glob
import os
import matplotlib.pyplot as plt

get_ipython().magic('matplotlib inline')


# In[224]:

names = []
for fname in glob.glob('./*fractions.csv'):
    names.append(fname)

name2nummap = {'five':5,'six':6,'seven':7,'eight':8,'nine':9}

aggdata = pd.DataFrame()
for filename in names:
    temp = pd.read_csv(filename,header=0, names=['model','fail_less','pass_greater'])

    temp['number'] = temp.pass_greater*1000/100
    temp['genotype'] = filename.split('/')[-1].split('_')[0]
    temp['denticles'] = int(filename.split('/')[-1].split('_')[1][:2])
    temp['stdev'] = name2nummap[filename.split('/')[-1].split('_')[2]]

    aggdata = aggdata.append(temp,ignore_index=True)

aggdata['acceptable'] = 'no'
aggdata.loc[aggdata.number >= 950, 'acceptable'] = 'YES'

sort_order = ['halfD', 'sixtenthsD', 'twothirdsD', 'sevententhsD', 'threequartersD','eighttenthsD', 'ninetenthsD', 'oneD']
sort_orderIndex = dict(zip(sort_order,range(len(sort_order))))
aggdata['modelrank'] = aggdata['model'].map(sort_orderIndex)

aggdata.set_index(['genotype','stdev','modelrank','model','denticles'],inplace=True)

aggdata.sort_index(level=['genotype','denticles','modelrank','stdev'],inplace=True)

ordereddata = aggdata[['acceptable','number','pass_greater','fail_less']]
ordereddata.to_csv('pretty_statistics_table.csv')


# In[227]:

usort = aggdata.reset_index()


# In[234]:

ax = usort[(usort.acceptable=='YES') & (usort.genotype=='yw')].plot(y='modelrank',x='denticles',kind='scatter',color='blue')
usort[(usort.acceptable=='YES') & (usort.genotype=='Larvae')].plot(y='modelrank',x='denticles',kind='scatter',color='red',ax=ax)


# In[ ]:



