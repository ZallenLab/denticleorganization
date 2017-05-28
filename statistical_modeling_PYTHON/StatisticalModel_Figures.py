from __future__ import division

import os
import glob
import numpy as np
import pandas as pd
import scipy.stats as sps
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib

import itertools





# get invivo data
input_file = 'yw_CellbyCell.csv'
invivo_data = pd.read_csv(input_file)

# get genotype from filename
genotype = input_file.split('_')[0]  # split the string at the underscore, returns a list of the parts (no underscores), take the first/before part

# clear out missing invivo_data
invivo_data = invivo_data.replace(0,np.nan)       # turn zeros into NaNs
invivo_data = invivo_data.dropna(how='all')       # drop any column (axis=0) or row (axis=1) where ALL values are NaN
invivo_data = invivo_data.replace(np.nan,0)       # turn NaNs back into zeros, so arithmetic can be normal for 'true-zero' values
invivo_data = invivo_data[invivo_data.dentincell != 1]   # drop columns where the cell has only 1 denticle









# CDF PLOTS
linestodraw = 1000


names = []
for file in glob.glob('*iteration*distance*.csv'):
    names.append(file)

models = ['halfD', 'twothirdsD', 'sevententhsD','oneD']
stdevs = ['six']

pairs = list(itertools.product(models,stdevs))


plt.rc('figure',figsize=(3,2))
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
myfig, ax = plt.subplots(1)



for name in names[:linestodraw]:
    data = pd.read_csv(name,index_col=[0,1,2])
    data.index.names = ['L1','L2','L3']
    data.columns.names = ['cellnumber']

    data = data.drop_duplicates()
    data = data.T

    gt, ite = name.split('_')[:2]

    for elem in itertools.product(models,stdevs):
        a,b = elem
        thisdata = data[a][b].values.flatten()
        thisdata = thisdata[np.nonzero(thisdata)]
        thisdata.sort()
        y = np.cumsum(thisdata)/thisdata.sum()*100

        ax.plot(thisdata,y,linewidth=1.5)

        # set axis ranges
        ax.set_xlim([0,10])
        ax.set_ylim([0,100])


# add invivo data
invd = invivo_data[invivo_data.columns[10:]].values.flatten()
invd = invd[np.nonzero(invd)]
invd.sort()
invd_y = np.cumsum(invd)/invd.sum()*100

ax.plot(invd,invd_y,color='k',linewidth=1.5)
ax.legend(models, loc='right', bbox_to_anchor=(1.5, 0.5))

# make set axis ranges
ax.set_xlim([0,10])
ax.set_ylim([0,100])

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

myfig.savefig(gt +'_'+ ite[9:] + '_' + linestodraw + '_cdf.svg',bbox='tight')
plt.close(myfig)








# SCATTERPLOTS
# TODO
# this is incomplete, but the scatter plots work if you feed them the correct information

names = []
for file in glob.glob('*iteration*distance*.csv'):
    names.append(file)


for name in names:
	data = pd.read_csv(name,index_col=[0,1,2])
	data.index.names = ['L1','L2','L3']
	data.columns.names = ['cellnumber']

	data = data.drop_duplicates()
	data = data.T



sub = data[data['basic info']['dentnumber'] == 2] #['basic info']['celllength']
subi = invivo_data[invivo_data['dentincell'] == 2]


plt.rc('figure',figsize=(5,5))

fig, axes = plt.subplots(3,2) #,sharex=True,sharey=True)

axes[0,0].scatter(sub['basic info']['celllength'],sub['twothirdsD']['five']['a'],color='r')
axes[0,0].scatter(subi['Dvlen'],subi['Intra1-2'],color='k') #,marker='x')
axes[0,0].set_title('five')
axes[0,0].set_xlabel('cell length (µm)')
axes[0,0].set_ylabel('separation (µm)')

axes[0,1].scatter(sub['basic info']['celllength'],sub['twothirdsD']['six']['a'],color='r')
axes[0,1].scatter(subi['Dvlen'],subi['Intra1-2'],color='k') #,marker='x')
axes[0,1].set_title('six')
axes[0,1].set_xlabel('cell length (µm)')
axes[0,1].set_ylabel('separation (µm)')

axes[1,0].scatter(sub['basic info']['celllength'],sub['twothirdsD']['seven']['a'],color='r')
axes[1,0].scatter(subi['Dvlen'],subi['Intra1-2'],color='k') #,marker='x')
axes[1,0].set_title('seven')
axes[1,0].set_xlabel('cell length (µm)')
axes[1,0].set_ylabel('separation (µm)')

axes[1,1].scatter(sub['basic info']['celllength'],sub['twothirdsD']['eight']['a'],color='r')
axes[1,1].scatter(subi['Dvlen'],subi['Intra1-2'],color='k') #,marker='x')
axes[1,1].set_title('eight')
axes[1,1].set_xlabel('cell length (µm)')
axes[1,1].set_ylabel('separation (µm)')

axes[2,0].scatter(sub['basic info']['celllength'],sub['twothirdsD']['nine']['a'],color='r')
axes[2,0].scatter(subi['Dvlen'],subi['Intra1-2'],color='k') #,marker='x')
axes[2,0].set_title('nine')
axes[2,0].set_xlabel('cell length (µm)')
axes[2,0].set_ylabel('separation (µm)')

# plt.scatter(sub['basic info']['celllength'],sub['twothirdsD']['eight']['a'],color='c')
plt.subplots_adjust(left=1,bottom=1,right=2,top=2,wspace=.25,hspace=.35)

plt.savefig('twothirdsD_yw.png',dpi=300,bbox='tight')







# REORGANIZE P VALUE TABLES
os.getcwd()

filenames = []
for file in glob.glob('*pvalue.csv'):
    filenames.append(file)

newframe = pd.DataFrame()

for name in filenames:
    data = pd.read_csv(name,index_col=0)
    # round to six decimal places (make things x*10^-(a zillion) equal to zero, as they are, effectively)
    data = data.round(6)

    genotype, stdev, _, _ = name.split('_')

#     a = data.columns + '_' + stdev
    data.columns = data.columns + '_' + stdev

    newframe = pd.concat([newframe, data], axis=1, join='inner')
    newframe.to_csv(genotype + '_allpvalues_rounded.csv')













# P-VALUE HISTOGRAMS
# TODO
# this is also incomplete, but bits might be useful




# gives an array with one plot for each alpha,sigma combination
filenames = []
for file in glob.glob('*pvalue.csv'):
    filenames.append(file)


for name in filenames:
    data = pd.read_csv(name,index_col=0)
    # round to six decimal places (make things x*10^-(a zillion) equal to zero, as they are, effectively)
    data = data.round(6)

    genotype, stdev, _, _ = name.split('_')

    elemweights = np.ones((10000))/10000

    fig = plt.figure()
    fig = data.hist(range=[0,1],xlabelsize=9,ylabelsize=9,sharey=True,figsize=(20,2),layout=(1,8),bins=20,weights=elemweights)

    plt.savefig(genotype + '_' + stdev+ '_pvaluefigure.png',dpi=300,bbox='tight')






# gives a single overlapped histogram for each sigma
# fig2 = data.plot(kind='hist',range=[0,1],figsize=(3,2),layout=(1,8),bins=20,alpha=0.8)
fig2 = data.plot(kind='hist',range=[0,1],figsize=(3,2),layout=(1,8),bins=20,alpha=0.8)
fig2.set_ylabel('% of simulations')
fig2.set_xlabel('p-value')
fig2.set_title(stdev)
# fig2 = data.plot(kind='hist',subplots=True,range=[0,1],figsize=(20,2),layout=(1,8),bins=20)
# plt.savefig('test.png',dpi=300,bbox='tight')



