
# coding: utf-8

# ### Variations in alpha over denticle number
# 
# Goal: graphing variations in calculated alpha with denticle number in embryos and larvae

# In[1]:

from IPython.core.display import HTML
HTML("<style>.container { width:90% !important; }</style>")


# In[2]:

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

get_ipython().magic('matplotlib inline')


# In[3]:

sns.set_context('notebook')
sns.set_style('darkgrid')
# sns.set_style('white')
# sns.set_style('ticks')


# In[4]:

numbermap = {'five':5, 'six':6, 'seven':7, 'eight':8, 'nine':9}
modelmap = {'halfD':.5, 'twothirdsD':.666, 'sixtenthsD':.6, 'sevententhsD':.7, 
            'threequartersD':.75, 'eighttenthsD':.8, 'ninetenthsD':.9, 'oneD':1}


# In[7]:

basepath = os.getcwd()
basepath

# In[52]:

iterations = os.getcwd().split('_')[-2][0:-10]


# In[53]:

names = []
for fname in glob.glob('*ks_fractions.csv'):
# for fname in glob.glob('*ks_fractions*.csv'):
    names.append(fname)


fract_data = pd.DataFrame()


for name in names: 
    genotype = name.split('_')[0]
    stdev = numbermap[name.split('_')[3]]

    temp = pd.read_csv(name, header=0, names=['model','less','greater'])
    temp['stage'] = genotype
    temp['stdev'] = stdev
    temp['dentnumber'] = int(name.split('_')[1][0:2])

    fract_data = pd.concat([fract_data, temp], ignore_index=True)

os.chdir(basepath)

fract_data['modelnumber'] = fract_data['model'].map(modelmap)
fract_data.head()


# ### the real graphs

# In[75]:

frac = 99

passes = fract_data[fract_data.greater >= frac]

with sns.axes_style('darkgrid'):
    fig, ax = plt.subplots(figsize=(15,5))
    ax.scatter(passes[passes.stage=='Larvae'].dentnumber, passes[passes.stage=='Larvae'].modelnumber, 
                color='r',marker='o')
    ax.scatter(passes[passes.stage=='yw'].dentnumber, passes[passes.stage=='yw'].modelnumber, 
                color='k',marker='o')
    
    ax.set_xticks(np.arange(0,18))
    ax.set_yticks([0.6, 0.667, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1])
    
    ax.set_xlim(1.5,12.5)
    ax.set_ylim(0.59,1.1)


    ax.set_title('alpha values per denticle number, where at least %i of simulations pass (%s iterations)' % (frac, iterations))
    ax.set_ylabel('model (alpha value)')
    ax.set_xlabel('denticle number')
    
    sns.despine()

fig.savefig(basepath + '/alpha values per denticle number, where at least %i of simulations pass (%s iterations)_grid.svg' % (frac, iterations))


# In[73]:

frac = 99
passes = fract_data[fract_data.greater >= frac]

with sns.axes_style('ticks'):
    fig, ax = plt.subplots(figsize=(15,5))

    ax.scatter(passes[passes.stage=='Larvae'].dentnumber, passes[passes.stage=='Larvae'].modelnumber, 
                color='r',marker='o')
    ax.scatter(passes[passes.stage=='yw'].dentnumber, passes[passes.stage=='yw'].modelnumber, 
                color='k',marker='o')

    ax.set_xticks(np.arange(0,18))
    ax.set_yticks([0.6, 0.667, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1])
    
    ax.set_xlim(1.5,12.5)
    ax.set_ylim(0.59,1.1)
    
    ax.set_title('alpha values per denticle number, where at least %i of simulations pass (%s iterations)' % (frac, iterations))
    ax.set_ylabel('model (alpha value)')
    ax.set_xlabel('denticle number')
    sns.despine()

fig.savefig(basepath +'/alpha values per denticle number, where at least %i of simulations pass (%s iterations).svg' % (frac, iterations))


# ### some checks and validation things

# In[46]:

fract_data.describe()

passes = fract_data[fract_data.greater >= 95]sd, dn = 8, 9

ax = fract_data[fract_data.stdev == sd][fract_data[fract_data.stdev == sd]
                                        .dentnumber == dn].plot('modelnumber','greater',kind='scatter')

ax.set_ylim(0,100)
sns.despine()
ax.set_title(str(dn )+ ' denticles and ' + str(sd) + ' standard deviations')ax = passes[passes.stage=='Larvae'].plot('dentnumber', 'modelnumber', kind='scatter')
ax.set_ylim(0,1.1)
sns.despine()

ax.set_title('Larvae')passes.sort_values(by=['dentnumber','stdev','modelnumber']).head()pd.set_option('display.max_rows',100)

pt = passes.set_index(['stage','model'])
# pt.sort_values(by=['dentnumber','modelnumber'])